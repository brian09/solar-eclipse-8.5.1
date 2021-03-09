#include <iostream>
#include "solar.h"
#include "plinkio.h"
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "gpu-pedigree-data.h"
#include "gpu-selection.h"

using namespace std;
static vector<string> read_id_list(const char * id_list_filename){
	vector<string> output;
	ifstream input_stream(id_list_filename);
	if(!input_stream.is_open()) return output;
	string id;
	while(input_stream >> id){
		output.push_back(id);
	}
	return output;
}
static inline int create_index(const int i, const int j, const int N)
{
   if (i <= j)
      return i * N - (i - 1) * i / 2 + j - i;
   else
      return j * N - (j - 1) * j / 2 + i - j;
}


static void write_pedigree_to_file(string output_filename, float * empirical_pedigree, pio_file_t * plink_file, \
                                    const int n_subjects, const int subset_size,const int * subset_index_map, const bool flatten_indices){
    ofstream output_stream(output_filename.c_str());
    
    pio_sample_t * sample_i;
    pio_sample_t * sample_j;
    output_stream << "IDA,IDB,KIN\n"; 
    if(subset_index_map == 0 ){
        int array_index = 0;
    	for(size_t j = 0; j < n_subjects ;j++){
        	sample_j = fam_get_sample(&plink_file->fam_file, j);
        	if (flatten_indices) 
        	    output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[array_index++] << endl;
        	else
        	    output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[j*n_subjects + j] << endl;
        	for(size_t i = j + 1; i < n_subjects ;i++){
            		sample_i = fam_get_sample(&plink_file->fam_file, i);
			        //const int index = create_index(i,j, n_subjects);
			        if(flatten_indices){
            		    output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  empirical_pedigree[array_index] << endl;
            		    output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[array_index] << endl;
            		    array_index++;
            		}else{
            		    output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  empirical_pedigree[j*n_subjects + i] << endl;
            		    output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[j*n_subjects + i] << endl; 
            		}           		    
        	}
    	}
    }else{
        int array_index = 0;
	    for(size_t j = 0 ; j < subset_size; j++){
		    sample_j = fam_get_sample(&plink_file->fam_file, subset_index_map[j]);
		    if(flatten_indices)
		        output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[array_index++] << endl;
		    else
		       output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[j*subset_size + j] << endl; 
        	for(size_t i = j + 1; i < subset_size ;i++){
            		sample_i = fam_get_sample(&plink_file->fam_file, subset_index_map[i]);
			      //  const int index = create_index(i,j, subset_size);
			        if (flatten_indices){
            		    output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  empirical_pedigree[array_index] << endl;
            		    output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[array_index] << endl;
            		    array_index++;
            		}else{
            		    output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  empirical_pedigree[j*subset_size + i] << endl;
            		    output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  empirical_pedigree[j*subset_size + i] << endl;            		    
            		}

        	}
	    }
    }		
    
    output_stream.close();
}
static string create_empirical_pedigree_correlation(vector<int> gpu_id_list,  int batch_size, int thread_size, const int per_chromosome, pio_file_t * plink_file, \
							const char * frequency_filename, const char * output_filename, const float alpha, int snp_stride, const bool normalize, \
							int * subset_index_map = 0 , const int subset_size = 0){
	if(batch_size % snp_stride != 0 && snp_stride != -1){
		std::cout << "Warning batch size value " << batch_size << " is not divisible by snp stride value " << snp_stride << "\n";
		batch_size = ceil(batch_size/snp_stride)*snp_stride;
		std::cout << "Setting batch size equal to " << batch_size << "\n";
	} 
	
	float * cpu_empirical_pedigree;
//	unsigned * cpu_missing_snp_count;
	const int n_subjects = pio_num_samples(plink_file);
//	cout << "n_subjects : " << n_subjects << " subset size: " << subset_size << endl; 
	unsigned array_size;
	int row_size;
	if(subset_size){
		row_size = subset_size;
		array_size = subset_size*(1+subset_size)/2;
	}else{
		row_size = n_subjects;
		array_size = n_subjects*(1+n_subjects)/2;
	}
	if(gpu_id_list.size() != 1){ 
	    try{		

		    cpu_empirical_pedigree = new float[array_size];
	    }catch(std::bad_alloc& e){
		
		    return string("Failed to allocate host memory for empirical pedigree");
	    }
	    /*try{ 
		    cpu_missing_snp_count = new unsigned[array_size];
	    }catch(std::bad_alloc& e){
		    delete [] cpu_empirical_pedigree;
		
		    return string("Failed to allocate host memory for missing snp count matrix");
	    }*/	
    }

	int pitch =ceil(n_subjects/32.f)*32;
	int subset_pitch = 0;
	if(subset_size){
		subset_pitch = ceil(subset_size/32.f)*32;
		pitch = subset_pitch;
	}
	GPU_Pedigree_Context * gpu_context;

	try{
		gpu_context =  new GPU_Pedigree_Context(plink_file,frequency_filename, batch_size,\
			thread_size,  n_subjects, snp_stride, \
			pitch, subset_size);
		thread_size = gpu_context->thread_size;
		snp_stride = gpu_context->snp_stride;
		batch_size = gpu_context->max_batch_size;
	}catch(GPU_Exception & e){
		if (gpu_id_list.size() != 1) delete [] cpu_empirical_pedigree;
	//	delete [] cpu_missing_snp_count;
	
		return string("Failed to allocate GPU memory for empirical pedigree or missing snp count matrix");
	}catch(std::bad_alloc& e){
		if (gpu_id_list.size() != 1) delete [] cpu_empirical_pedigree;
	//	delete [] cpu_missing_snp_count;
		
		return string("Failed to allocate host buffer for reading the plink files");
	}
		
	GPU_Pedigree_Data ** pedigree_data = new GPU_Pedigree_Data*[gpu_id_list.size()];

	for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
	
		try{
			if(subset_size == 0){
				pedigree_data[gpu_index] = new GPU_Pedigree_Data(gpu_id_list[gpu_index],  alpha, thread_size,\
				                                                 batch_size,plink_file->bed_file.header.num_loci,\
				                                                 n_subjects, pitch, array_size, 0, 0);
			}else{
				pedigree_data[gpu_index] = new GPU_Pedigree_Data(gpu_id_list[gpu_index],  alpha, thread_size,\
				                                                 batch_size, plink_file->bed_file.header.num_loci,\
				                                                 n_subjects, pitch, array_size,subset_index_map,\
				                                                 subset_size);	
			}			
		}catch(GPU_Exception & e){
			cout << e.what() << endl;
			if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
		//	if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;
			
			delete gpu_context;
			string error_message = "Failed to allocate GPU memory for either storing SNP values, allele frequencies, the empirical pedigree, or matrix of missing SNP counts.";
			return error_message;
		}catch(std::bad_alloc& e){
			if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
		//	if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;
		
			delete gpu_context;
			string error_message = "Failed to allocate memory for host array that stores " + to_string(batch_size) + " sets of SNP data.";
			return error_message;
		}
	}
	
	string errmsg;
	if(per_chromosome){
        	struct pio_locus_t * locus;
       	 	locus = bim_get_locus(&plink_file->bim_file, 0);		
	 	unsigned char chromosome = locus->chromosome;						
		int start = 0;
		int end = 0;
		while(start < plink_file->bed_file.header.num_loci){
			cout << "Running GRM creation for chromosome " << chromosome << endl;	
			for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
				try{
					cudaErrorCheck(cudaSetDevice(pedigree_data[gpu_index]->gpu_id));
					cudaErrorCheck(cudaMemset(pedigree_data[gpu_index]->gpu_empirical_pedigree, 0, sizeof(float)*row_size*row_size));
			//		cudaErrorCheck(cudaMemset(pedigree_data[gpu_index]->gpu_missing_snp_count, 0, sizeof(unsigned)*array_size));
				}catch(GPU_Exception & e){
					if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
			//		if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;
					errmsg = "Failure to initialize GPU pedigree arrays to zero on GPU ID " + to_string(pedigree_data[gpu_index]->gpu_id);
					return errmsg;
				}
			}				
			unsigned snp = start;
			while(locus->chromosome == chromosome && snp < plink_file->bed_file.header.num_loci) locus = bim_get_locus(&plink_file->bim_file, snp++);
			end = snp;
			unsigned snp_batch_size = end - start;
			cout << "Chromosome " << chromosome << " contains " << snp_batch_size << " loci\n";
			for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++) pedigree_data[gpu_index]->total_snps = snp_batch_size;
			errmsg = gpu_context->run_pedigree_creation(pedigree_data, gpu_id_list.size(), snp_batch_size);
			if(errmsg.length() != 0){   
				if(gpu_id_list.size() != 1) delete [] cpu_empirical_pedigree;
		//		if(gpu_id_list.size() != 1) delete [] cpu_missing_snp_count;

				return errmsg;
			}
			if(gpu_id_list.size() != 1){
				memset(cpu_empirical_pedigree, 0, sizeof(float)*array_size);
			//	memset(cpu_missing_snp_count, 0, sizeof(unsigned)*array_size);
				errmsg = gpu_context->combine_gpu_results(pedigree_data, cpu_empirical_pedigree, gpu_id_list.size());
				if(errmsg.length() != 0){
					delete [] cpu_empirical_pedigree;
				//	delete [] cpu_missing_snp_count;
					return errmsg;
				}
				int array_index = 0;
			    for(int col = 0 ; col < row_size; col++){
				    for(int row = col; row < row_size; row++){
					    cpu_empirical_pedigree[array_index] /= (snp_batch_size);// - cpu_missing_snp_count[array_index]);
					    array_index++;
				    }
			    }
			}else{
				try{
				 	gpu_context->copy_results_to_cpu(pedigree_data[0]);
				 }catch(GPU_Exception & e){
				 	errmsg = "Failed to copy GPU results to CPU";
				 	if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
				 //	if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;
				 	return errmsg;
				 }
				 int array_index = 0;
				 for(int col = 0; col < row_size; col++){
				    for(int row = col; row < row_size; row++){  
				        gpu_context->temp_cpu_empirical_pedigree[col*row_size + row] /= (snp_batch_size);// - gpu_context->temp_cpu_missing_snp_count[array_index]);
				        array_index++;
				    }
				 }  
				  
				// memcpy(cpu_empirical_pedigree, gpu_context->temp_cpu_empirical_pedigree, sizeof(float)*array_size);
				// memcpy(cpu_missing_snp_count, gpu_context->temp_cpu_missing_snp_count, sizeof(unsigned)*array_size);
			}
							

			if(normalize){
			    
				float norms[row_size];
				if(gpu_id_list.size() != 1){
				    int array_index = 0;
				    for(int index =0 ; index < row_size; index++){
					    norms[index] = sqrt(cpu_empirical_pedigree[array_index]);
					    array_index += row_size - index;
					    
				    }
				    array_index = 0;
				    for(int col = 0; col < row_size; col++){
					    for(int row = col; row < row_size; row++){
						    cpu_empirical_pedigree[array_index]/(norms[row]*norms[col]);
						    array_index++;
					    }
				    }
				}else{
				    for(int index =0 ; index < row_size; index++){
					    norms[index] = sqrt(gpu_context->temp_cpu_empirical_pedigree[index*row_size + index]);
				    }
				    for(int col = 0; col < row_size; col++){
					    for(int row = col; row < row_size; row++){
						    gpu_context->temp_cpu_empirical_pedigree[col*row_size + row]/(norms[row]*norms[col]);
					    }
				    }
				}				    
			}
			char name_buffer[200];
            		sprintf(name_buffer, "%s.chr%u.csv", output_filename, (unsigned int )chromosome);
			string filename = string(name_buffer);
		
			if(gpu_id_list.size() == 1){
			    
			    write_pedigree_to_file(filename, gpu_context->temp_cpu_empirical_pedigree, plink_file,n_subjects, subset_size, subset_index_map, false);
			}else{
			    write_pedigree_to_file(filename, cpu_empirical_pedigree, plink_file,n_subjects, subset_size, subset_index_map, true);
			}
		        if(end != plink_file->bed_file.header.num_loci){
                		chromosome = locus->chromosome;
            		}
            		start = end;
            		cout << "GRM creation completed for chromosome " << chromosome << "\n";
        }
    }else{
            	
         errmsg = gpu_context->run_pedigree_creation(pedigree_data, gpu_id_list.size(),  plink_file->bed_file.header.num_loci);
		if(errmsg.length() != 0){
			if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
	//		if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;

			return errmsg;
		}
		if(gpu_id_list.size() != 1){
			memset(cpu_empirical_pedigree, 0, sizeof(float)*array_size);
		//	memset(cpu_missing_snp_count, 0, sizeof(unsigned)*array_size);
			errmsg = gpu_context->combine_gpu_results(pedigree_data, cpu_empirical_pedigree, gpu_id_list.size());
			if(errmsg.length() != 0){
				delete [] cpu_empirical_pedigree;
			//	delete [] cpu_missing_snp_count;
				return errmsg;
			}
		    int array_index = 0;
			for(int col = 0 ; col < row_size; col++){
		        for(int row = col; row < row_size; row++){
			        cpu_empirical_pedigree[array_index] /= (plink_file->bed_file.header.num_loci);// - cpu_missing_snp_count[array_index]);
			        array_index++;
				}
			}			
		}else{
			try{
				 gpu_context->copy_results_to_cpu(pedigree_data[0]);
			}catch(GPU_Exception & e){
				 errmsg = "Failed to copy GPU results to CPU";
				 if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
			//	 if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;
				 return errmsg;
			}
            int array_index =0;
		    for(int col = 0 ; col < row_size; col++){
			    for(int row = col; row < row_size; row++){
				    gpu_context->temp_cpu_empirical_pedigree[col*row_size + row] /= (plink_file->bed_file.header.num_loci);// - gpu_context->temp_cpu_missing_snp_count[array_index]);
				    array_index++;
			    }
		    }			
		//	memcpy(cpu_empirical_pedigree, gpu_context->temp_cpu_empirical_pedigree, sizeof(float)*array_size);
		//	memcpy(cpu_missing_snp_count, gpu_context->temp_cpu_missing_snp_count, sizeof(unsigned)*array_size);
		}       		
            //	cout << cpu_empirical_pedigree[0] << " " << cpu_missing_snp_count[0] << endl;

		
		if(normalize){
			float norms[row_size];
			if (gpu_id_list.size() != 1){
			    int array_index = 0;
			    for(int index =0 ; index < row_size; index++){
				    norms[index] = sqrt(abs(cpu_empirical_pedigree[array_index]));
				    array_index += row_size - index;
			    }
			    array_index = 0;
			    for(int col = 0; col < row_size; col++){
				    for(int row = col; row < row_size; row++){
					    cpu_empirical_pedigree[array_index] /= norms[row]*norms[col];
					    array_index++;
				    }
			    }
			}else{
		
			    for(int index =0 ; index < row_size; index++){
				    norms[index] = sqrt(abs(gpu_context->temp_cpu_empirical_pedigree[index*row_size + index]));
				  
			    }
			    for(int col = 0; col < row_size; col++){
				    for(int row = col; row < row_size; row++){
					    gpu_context->temp_cpu_empirical_pedigree[col*row_size + row] /= norms[row]*norms[col];
				    }
			    }			
			}
		}
	    if(gpu_id_list.size() == 1){
			    
	        write_pedigree_to_file(string(output_filename), gpu_context->temp_cpu_empirical_pedigree, plink_file,n_subjects, subset_size, subset_index_map, false);
		}else{
			write_pedigree_to_file(string(output_filename), cpu_empirical_pedigree, plink_file,n_subjects, subset_size, subset_index_map, true);
		}		

		cout << "GRM creation completed\n";
	}


	if(gpu_id_list.size() != 1)  delete [] cpu_empirical_pedigree;
//	if(gpu_id_list.size() != 1)  delete [] cpu_missing_snp_count;

	delete [] pedigree_data;
	delete gpu_context;	
	//if(subset_size != 0) delete [] gpu_subset_index_map;
	string empty_string;
	return empty_string;				            		
								
}
static std::string reset_devices(vector<int> gpu_id_list){
	string error_message;
	for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
		try{
			cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
			cudaErrorCheck(cudaDeviceReset());
		}catch(GPU_Exception & e){
			error_message = "Failed to reset device with GPU ID " + to_string(gpu_id_list[gpu_index]);
			return error_message;
		}
	}
	return error_message;
}
extern "C" int gpupedfromsnpsCmd(ClientData clientData, Tcl_Interp *interp,
                              int argc,const char *argv[]){
    

    int per_chromosome = 0;
    const char * plink_filename = 0;
    const char * output_filename = 0;
    const char * frequency_filename = 0;
    float alpha = -1.f;
    const char * id_list_filename = 0;
    bool select_all_gpus = false;
    const char * gpu_list_string = 0;  
    int batch_size = 1000;
    int thread_size = -1;
    int snp_stride = -1;
    bool normalize = false;
    for(unsigned arg = 1; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--i", case_ins) || \
            !StringCmp(argv[arg], "--input", case_ins) \
            || !StringCmp(argv[arg], "-i", case_ins)\
            || !StringCmp(argv[arg], "-input", case_ins)) && arg + 1 < argc){
            plink_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--o", case_ins) || \
                  !StringCmp(argv[arg], "--output", case_ins) \
                  || !StringCmp(argv[arg], "-o", case_ins)\
                  || !StringCmp(argv[arg], "-output", case_ins)) && arg + 1 < argc){
            output_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "-corr", case_ins) || \
                  !StringCmp(argv[arg], "--corr", case_ins)) && arg + 1 < argc){
            alpha = atof(argv[++arg]);
           
        }else if((!StringCmp(argv[arg], "-freq", case_ins) || \
                  !StringCmp(argv[arg], "--freq", case_ins)) && arg + 1 < argc){
            frequency_filename = argv[++arg];
           
        }else if ((!StringCmp(argv[arg], "--gpus", case_ins) || !StringCmp(argv[arg], "-gpus", case_ins)) && arg + 1 < argc){
            
            gpu_list_string = argv[++arg];
            
        }else if (!StringCmp(argv[arg], "--all", case_ins) || !StringCmp(argv[arg], "-all", case_ins)){
            
            select_all_gpus = true;
            
        }else if (!StringCmp(argv[arg], "--normalize", case_ins) || !StringCmp(argv[arg], "-normalize", case_ins)){
            
            normalize = true;
            
        }else if((!StringCmp(argv[arg], "-batch_size", case_ins) || \
                  !StringCmp(argv[arg], "--batch_size", case_ins)) && arg + 1 < argc){
            batch_size = atoi(argv[++arg]);
            if(batch_size < 1){
            	RESULT_LIT("-batch_size must be greater than 0");
            	return TCL_ERROR;
            }
            if(batch_size > 20000){
            	RESULT_LIT("-batch_size must be less than or equal to 20000");
            	return TCL_ERROR;
            }            
        }else if((!StringCmp(argv[arg], "-thread_size", case_ins) || \
                  !StringCmp(argv[arg], "--thread_size", case_ins)) && arg + 1 < argc){
            thread_size = atoi(argv[++arg]);
          /*  if(thread_size != 32 && thread_size != 64 && thread_size != 128 && \
             thread_size != 256 && thread_size != 512 && thread_size != 1024 ){
            	RESULT_LIT("-thread_size must be 32,64,128,256,512, or 1024");
            	return TCL_ERROR;
            } */
            if(thread_size <= 0 || thread_size > 1024 || thread_size % 32 != 0){
            	RESULT_LIT("-thread_size must be a multiple of 32 greater than or equal to 32 and less than or equal to 1024");
            	return TCL_ERROR;
            }           
        }else if((!StringCmp(argv[arg], "-snp_stride", case_ins) || \
                  !StringCmp(argv[arg], "--snp_stride", case_ins)) && arg + 1 < argc){
            snp_stride = atoi(argv[++arg]);
            if(snp_stride < 1 || snp_stride > 10){
            	RESULT_LIT("-snp_stride must be between 1 and 10");
            	return TCL_ERROR;
            }
           
        }else if(!StringCmp(argv[arg], "-per-chromo", case_ins) || \
                 !StringCmp(argv[arg], "--per-chromo", case_ins)){
            per_chromosome = 1;
        }else if((!StringCmp(argv[arg], "-id_list", case_ins) || \
                  !StringCmp(argv[arg], "--id_list", case_ins)) && arg + 1 < argc){
            id_list_filename = argv[++arg];
        }else{
            string error_message = "Invalid argument was entered: " + string(argv[arg]);
            cout << error_message << endl;
           // RESULT_BUF(error_message.c_str());
            return TCL_ERROR;
        }
    }
    vector<int> gpu_id_list;
    try{
     	gpu_id_list = get_gpu_id_list(gpu_list_string, select_all_gpus);
     }catch(GPU_Exception & e){
     	RESULT_BUF(e.what());
     	return TCL_ERROR;
     }catch(...){
     	RESULT_LIT("Unknown error occurred creating GPU list");
     	return TCL_ERROR;
     }
    
    if(gpu_id_list.size() == 0){
        RESULT_LIT("No GPUs were selected");
        return TCL_ERROR;
    }    
    
    if(!plink_filename){
         cout << "No plink file set specified with --i"<< endl;
        return TCL_ERROR;
    }
    
    if(!output_filename){
         cout << "No output filename specified with --o"<< endl;
        return TCL_ERROR;
    }
    if(!frequency_filename){
    	cout << "No frequency filename specified with --freq" << endl;
    	return TCL_ERROR;
    }
    
    pio_file_t * plink_file = new pio_file_t;
    pio_status_t status;
    status = pio_open(plink_file, plink_filename);
    if(status != PIO_OK){
        if(status == P_FAM_IO_ERROR){
            cout << "Error in loading .fam file"<< endl;
            return TCL_ERROR;
        }else if (status == P_BIM_IO_ERROR){
             cout << "Error in loading .bim file"<< endl;
            return TCL_ERROR;
        }else if (status == P_BED_IO_ERROR){
             cout << "Error in loading .bed file"<< endl;
            return TCL_ERROR;
        }else{
             cout << "Error loading plink file"<< endl;;
            return TCL_ERROR;
        }
    }
    
    std::cout << "Number of subjects: " << pio_num_samples(plink_file) << " Number of loci: " << plink_file->bed_file.header.num_loci << std::endl;
    int version = plink_file->bed_file.header.version;
    if (plink_file->bed_file.header.snp_order == BED_UNKNOWN_ORDER){
        pio_close(plink_file);
        cout << "Error in the .bed snp order. Retry creation of file using a different version of plink" << endl;
        delete plink_file;
        return TCL_ERROR;
    }
    ifstream test_stream(frequency_filename);
    if(test_stream.is_open() == false){
    	cout << "Failed to open frequency file: " << frequency_filename << endl;
    	pio_close(plink_file);
    	delete plink_file;
    	return TCL_ERROR;
    }
    string first_line;
    getline(test_stream,first_line);
    const int freq_count = stoi(first_line);
    if(freq_count !=  plink_file->bed_file.header.num_loci){
    	cout << "Number of frequencies in frequency file: " << freq_count << endl;
    	cout << "Number of frequencies in frequency file must be equal to the number of loci in the specfieid plink file\n";
    	pio_close(plink_file);
    	delete plink_file;
    	return TCL_ERROR;
    } 
    test_stream.close();
    /*else if (plink_file->bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
        pio_close(plink_file);
        printf("In order to read efficiently the transpose of specified plink file must be performed\n");
        string transpose_filename = string(plink_filename) + string(".trans");
        string message = "Filename of transposed plink file is " + transpose_filename + "\n";
        printf(message.c_str());
        status = pio_transpose(plink_filename, transpose_filename.c_str());
        if(status != PIO_OK){
            RESULT_LIT("Error in creating transpose");
            return TCL_ERROR;
        }
        
        status = pio_open(plink_file, transpose_filename.c_str());
        
        if(status != PIO_OK){
            printf("Error in opening transposed plink file\n");
            return TCL_ERROR;
        }
    }*/
    vector<string> id_vector; 
    int * subset_index_map = 0;
    int subset_size = 0;
    if(id_list_filename){
	  id_vector = read_id_list(id_list_filename);
	  if(id_vector.size() == 0 ){
		cout << "No IDs read from id list option" << endl;
		pio_close(plink_file);
		delete plink_file;
		return TCL_ERROR;
	  }
	  vector<int> id_index_map;
	  pio_sample_t * sample; 
	  for(int index = 0; index < pio_num_samples(plink_file); index++){
		sample = pio_get_sample(plink_file, index);
		string str_id = string(sample->iid);
		vector<string>::iterator find_iter = find(id_vector.begin(), id_vector.end(), str_id);
		if(find_iter != id_vector.end()){
			id_index_map.push_back(index);
		}
	}
	if(id_index_map.size() == 0){
		cout  << "None of the IDs listed were found within the plink file" << endl;
		pio_close(plink_file);
		delete plink_file;
		
		return TCL_ERROR;
	}
	subset_index_map = new int[id_index_map.size()];
	for(int index = 0 ; index < id_index_map.size(); index++){
		subset_index_map[index] = id_index_map[index];
	}
	subset_size = id_index_map.size();
	id_index_map.clear();
   } 
   id_vector.clear();
   if(subset_size != 0 ) std::cout << "Number of subjects in subset: " << subset_size << std::endl;
   string errmsg;
   if(normalize) cout << "Final values will be normalized so diagonal elements are all one and off diagonal elements are bounded by one and negative one\n";   	
    cout << "Starting GPU GRM creation\n";      
    errmsg = create_empirical_pedigree_correlation(gpu_id_list,  batch_size, thread_size,  per_chromosome,  plink_file, \
					      frequency_filename,output_filename,  alpha, snp_stride, normalize, \
					      subset_index_map, subset_size);
    if(subset_index_map) delete [] subset_index_map;
    pio_close(plink_file);
    delete plink_file;							 
							 
     if(errmsg.length() != 0){
    	cout << errmsg << endl;
    	errmsg = reset_devices(gpu_id_list);
    	if(errmsg.length() != 0) cout << errmsg << endl;
    	return TCL_ERROR;
    }						 
							 
    errmsg = reset_devices(gpu_id_list);
    if(errmsg.length() != 0){
    	cout << errmsg << endl;
    	
    	return TCL_ERROR;
    	
    }
    		
   
						
							 

  //  RESULT_LIT("Empirical pedigree creation is complete");
//else{
   	 return TCL_OK;
   //}
}

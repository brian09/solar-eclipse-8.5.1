//
//  gpu-gwas.cc
//  
//
//  Created by Brian Donohue on 8/21/18.
//
#include <chrono>
#include "gpu-gwas.h"
#include "gpu-gwas-estimator.h"
#include <fstream>
#include "Eigen/Dense"
#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>
using namespace std;

static size_t get_total_adjustable_memory(const int pitch){

    return GPU_GWAS_Estimator_SNP_Data::Adjustable_GPU_Memory_Cost(pitch) + GPU_GWAS_Estimator_Context::Adjustable_GPU_Memory_Cost(pitch) + GPU_GWAS_Estimator_Data::Adjustable_GPU_Memory_Cost();
}

static size_t get_total_static_memory( const int pitch){
  //  if(multi){
   //     return  Multi_Trait_GPU_GWAS_Estimator::Static_GPU_Memory_Cost(n_subjects, pitch);
    //}else{
    	
        return GPU_GWAS_Estimator_Context::Shared_GPU_Memory_Cost(pitch) + GPU_GWAS_Estimator_SNP_Data::Shared_GPU_Memory_Cost(pitch);
    //}
    
}

void determine_GPU_GWAS_sizes(bool screen,vector<int> gpu_id_list, int n_subjects, int n_snps, \
                              int & max_batch_size, int & blockSize,  int & pitch, \
                              const int n_streams, const int defined_blockSize, const int defined_batch_size){

    int blockSize_list[6] = {32, 64, 128, 256, 512, 1024};
    int remainder_list[6];
    int pitch_list[6];
    int scale_list[6];
  //  cout << "n_subjects " << n_subjects << endl;
    for(int index = 0; index < 6; index++){
	scale_list[index] = ceil((double)n_subjects/blockSize_list[index]);
	//cout << "blockSize " << blockSize_list[index] << endl;
	//cout << "scale " << scale_list[index] << endl;

	pitch_list[index] = ceil(n_subjects/32.f)*32;
	//cout << "pitch " << pitch_list[index] << endl;
	remainder_list[index] = scale_list[index]*blockSize_list[index] - n_subjects;
	//cout << "remainder " << remainder_list[index] << endl;
    }
    int selection_index = 0;
    for(int index = 1; index < 6; index++){
	if(remainder_list[selection_index] > remainder_list[index]) selection_index = index;
    }
   // cout << "selection index " << selection_index << endl;
    if(scale_list[selection_index] >= 20 && selection_index < 5){
	int new_selection_index;// = selection_index + 1;
	for(new_selection_index = selection_index + 1; new_selection_index < 6; new_selection_index++){
		if(scale_list[new_selection_index] < 20) {
			selection_index = new_selection_index;
			break;
		}
	}
    }
//cout << "selection index " << selection_index << endl;
    if(defined_blockSize != -1){
	for(int index = 0; index < 6; index++){
		if(blockSize_list[index] == defined_blockSize){
			selection_index = index;
			break;
		}

	}
	if(defined_blockSize != blockSize_list[selection_index])\
	std::cout << "Block size selected is not accepted using block size " << blockSize_list[selection_index] << std::endl;
    }
    blockSize = blockSize_list[selection_index];
    pitch = pitch_list[selection_index];
    
    size_t minimum_gpu_memory = 0;
    
    for(size_t idx = 0; idx < gpu_id_list.size(); idx++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[idx]));
        size_t total_mem;
        size_t free_mem;
        cudaErrorCheck(cudaMemGetInfo(&free_mem,&total_mem));
        if(minimum_gpu_memory == 0) minimum_gpu_memory = free_mem;
        if(minimum_gpu_memory > free_mem) minimum_gpu_memory = free_mem;
        
    }
    
   // const size_t multiple = ceil(n_subjects/32.f);
   // pitch = multiple*32;
    size_t static_memory;
    if(!screen) 
    	static_memory = get_total_static_memory(pitch);
    else
    	static_memory = (pitch*pitch + 2*pitch)*sizeof(float);
    size_t dynamic_memory;
    if(!screen) 
    	dynamic_memory = get_total_adjustable_memory(pitch);
    else
    	dynamic_memory = 4*sizeof(float);
    
    size_t temp_batch_size = floor(double(minimum_gpu_memory - static_memory)/double(n_streams*dynamic_memory)) - 100;
    if(temp_batch_size > MAX_BATCH_SIZE)
        temp_batch_size = MAX_BATCH_SIZE;
    if(defined_batch_size != -1) temp_batch_size = defined_batch_size;
    if(temp_batch_size*n_streams*gpu_id_list.size() > n_snps){
        temp_batch_size = ceil(n_snps/double(n_streams*gpu_id_list.size()));
    }
    
    
    max_batch_size = temp_batch_size;
    

}
extern double gwas_chicdf(double, double);
gwas_data gwas_maximize_newton_raphson_method_null_model(Eigen::VectorXd Y, Eigen::VectorXd mean_column, Eigen::MatrixXd U, const int precision);

GPU_GWAS::~GPU_GWAS(){
  //  output_stream.close();
   // delete [] parameters;
  //  delete [] loglik;

    pio_close(plink_file);
    delete plink_file;
    delete reader;
    plink_ids.clear();
}


chrono::milliseconds GPU_GWAS::calibrate_blockSize(const int set_index, int & best_blockSize, double snp_percentage, int selected_batch_size = -1){
	chrono::milliseconds shortest_time (0);
	if(set_index >= n_sets){ 
		best_blockSize = 0;
		return shortest_time;
	}
	
    Eigen_Data * current_eigen_data_set = reader->get_eigen_data_set(set_index);
    n_subjects = current_eigen_data_set->get_n_subjects();
    int * index_map = new int[plink_ids.size()];
    vector<string> ids = current_eigen_data_set->get_ids();
    for(int index = 0; index < plink_ids.size(); index++){
        string plink_id = plink_ids[index];
        vector<string>::iterator find_iter = find(ids.begin(), ids.end(), plink_id);
        if(find_iter != ids.end()){
            index_map[index] = distance(ids.begin(), find_iter);
        }else{
            index_map[index] = -1;
        }
    }	
    Eigen::MatrixXd eigen_eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>\
    (current_eigen_data_set->get_eigenvectors_transposed(), n_subjects, n_subjects );
    Eigen::VectorXd eigen_mean_column = eigen_eigenvectors_transposed*Eigen::ArrayXd::Ones(n_subjects).matrix();
    Eigen::VectorXd eigen_component_A = eigen_mean_column.cwiseAbs2();
    Eigen::MatrixXd U = Eigen::ArrayXXd::Ones(n_subjects, 2);
    U.col(1) =  Eigen::Map<Eigen::VectorXd>(current_eigen_data_set->get_eigenvalues(), n_subjects );
    best_blockSize = 32;
    for(int size = 32; size <= 1024; size *= 2){
	    determine_GPU_GWAS_sizes(false,gpu_id_list, n_subjects, ceil(snp_percentage*n_snps), \
                             max_batch_size, blockSize, pitch, \
                              n_streams, size, selected_batch_size);	
    	    GPU_GWAS_TYPE * CPU_mean_column_component_A = new GPU_GWAS_TYPE[pitch*2];
    	    memset(CPU_mean_column_component_A,0.f,sizeof(GPU_GWAS_TYPE)*pitch*2);
            for(int row = 0; row < n_subjects; row++){
		CPU_mean_column_component_A[row] = eigen_mean_column(row);
		CPU_mean_column_component_A[pitch + row] = eigen_component_A(row);
   	    }

    	    GPU_GWAS_TYPE * CPU_eigenvectors_transposed = new GPU_GWAS_TYPE[pitch*pitch];
            memset(CPU_eigenvectors_transposed, 0.f, sizeof(GPU_GWAS_TYPE)*pitch*pitch);

            GPU_GWAS_TYPE * CPU_eigenvalues = new GPU_GWAS_TYPE[pitch]; 
            memset(CPU_eigenvalues, 0, sizeof(GPU_GWAS_TYPE)*pitch);
           for(int col = 0; col < n_subjects; col++){

        	CPU_eigenvalues[col] =(GPU_GWAS_TYPE) U(col, 1);
        	 for(int row = 0; row < n_subjects; row++){
            		CPU_eigenvectors_transposed[col*pitch + row] = (GPU_GWAS_TYPE)eigen_eigenvectors_transposed(row, col);
        	 }
   	  }
    	Eigen::VectorXd eigen_Y;
   	 Eigen::VectorXd eigen_component_D;
    	GPU_GWAS_TYPE * CPU_Y_component_D = new GPU_GWAS_TYPE[pitch*2];
   	 memset(CPU_Y_component_D, 0, sizeof(GPU_GWAS_TYPE)*pitch*2);
    	GPU_GWAS_Estimator * estimator;	
    //try{
	estimator = new GPU_GWAS_Estimator(gpu_id_list, plink_file, index_map,CPU_mean_column_component_A,\
                                        CPU_eigenvectors_transposed, CPU_eigenvalues, \
                                          precision, \
                                        max_batch_size,  blockSize, \
                                         ceil(snp_percentage*n_snps), n_subjects, pitch, n_streams, verbose);	
   /* }catch(exception & e){
	string error_message = string(e.what());
	throw GPU_GWAS_Exception(error_message);
    }catch(...){
	string error_message = "Unknown error occurred initializing gpu gwas estimator data";
	throw GPU_GWAS_Exception(error_message);
    }*/
    for(int index = 0; index < 1; index++){
	
              eigen_Y = Eigen::Map<Eigen::VectorXd>( current_eigen_data_set->get_phenotype_column(index), n_subjects);
	      eigen_Y = eigen_eigenvectors_transposed*eigen_Y;
	      eigen_component_D = eigen_Y.cwiseProduct(eigen_mean_column);
        
              double SD;
       	      gwas_data null_result = gwas_maximize_newton_raphson_method_null_model(eigen_Y, eigen_mean_column,U, precision);
       	      for(int row =0 ;row< n_subjects; row++){
     
                   CPU_Y_component_D[row] = (GPU_GWAS_TYPE)eigen_Y(row);
	           CPU_Y_component_D[row + pitch] = eigen_component_D(row);
      	      }


		pio_reset_row(plink_file);
		//try{
			auto start = std::chrono::high_resolution_clock::now();
			//std::string test_filename = "gpu-gwas-calibrate.out";
        	   	estimator->run(current_eigen_data_set->get_trait_name(index), CPU_Y_component_D, "calibrate", null_result);
			auto finish = std::chrono::high_resolution_clock::now();
   
                       chrono::milliseconds run_time = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
			if(size == 32) {
			   shortest_time = run_time;
			   best_blockSize = 32;
			}
			if(run_time < shortest_time){
				shortest_time = run_time;	
				best_blockSize = blockSize;
			}
		/*}catch(exception & e){
		 	string error_message = string(e.what());
			throw GPU_GWAS_Exception(error_message);
    		}catch(...){
			string error_message = "Unknown error occurred during gpu gwas parameter estimation";
			throw GPU_GWAS_Exception(error_message);
    		}*/
        	//write_output(current_eigen_data_set->get_trait_name(index), null_result);
	}
        delete estimator;
   // }
//auto finish = std::chrono::high_resolution_clock::now();
   
     // auto microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
      //  std::cout << "time " << microseconds.count() << "ms\n";
    
    	delete [] CPU_Y_component_D;
    	delete [] CPU_mean_column_component_A;
  //  delete [] CPU_component_A_D;
    	delete [] CPU_eigenvectors_transposed;
   // delete [] CPU_mean_column;
    	delete [] CPU_eigenvalues;
    }
    delete [] index_map;
    return shortest_time;	    
		
}
int GPU_GWAS::process_next_trait_set(){
    if(set_index >= n_sets) {
        return 1;
    }
auto start = std::chrono::high_resolution_clock::now();
    Eigen_Data * current_eigen_data_set = reader->get_eigen_data_set(set_index);
    n_subjects = current_eigen_data_set->get_n_subjects();
    determine_GPU_GWAS_sizes(false,gpu_id_list, n_subjects, n_snps, \
                             max_batch_size, blockSize, pitch, \
                              n_streams, defined_blockSize, defined_batch_size);
    

    
    int * index_map = new int[plink_ids.size()];
    vector<string> ids = current_eigen_data_set->get_ids();
    for(int index = 0; index < plink_ids.size(); index++){
        string plink_id = plink_ids[index];
        vector<string>::iterator find_iter = find(ids.begin(), ids.end(), plink_id);
        if(find_iter != ids.end()){
            index_map[index] = distance(ids.begin(), find_iter);
        }else{
            index_map[index] = -1;
        }
    }
    Eigen::MatrixXd eigen_eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>\
    (current_eigen_data_set->get_eigenvectors_transposed(), n_subjects, n_subjects );
    
    Eigen::VectorXd eigen_mean_column = eigen_eigenvectors_transposed*Eigen::ArrayXd::Ones(n_subjects).matrix();
    Eigen::VectorXd eigen_component_A = eigen_mean_column.cwiseAbs2();
    GPU_GWAS_TYPE * CPU_mean_column_component_A = new GPU_GWAS_TYPE[pitch*2];
    memset(CPU_mean_column_component_A,0.f,sizeof(GPU_GWAS_TYPE)*pitch*2);
    for(int row = 0; row < n_subjects; row++){
	CPU_mean_column_component_A[row] = eigen_mean_column(row);
	CPU_mean_column_component_A[pitch + row] = eigen_component_A(row);
    }
    Eigen::MatrixXd U = Eigen::ArrayXXd::Ones(n_subjects, 2);
    U.col(1) =  Eigen::Map<Eigen::VectorXd>(current_eigen_data_set->get_eigenvalues(), n_subjects );
    GPU_GWAS_TYPE * CPU_eigenvectors_transposed = new GPU_GWAS_TYPE[pitch*pitch];
    memset(CPU_eigenvectors_transposed, 0.f, sizeof(GPU_GWAS_TYPE)*pitch*pitch);
   // double * CPU_mean_column = new double[n_subjects];
  //  double * CPU_component_A_D = new double[pitch*2];
  //  memset(CPU_component_A_D, 0.f , sizeof(double)*pitch*2);
    GPU_GWAS_TYPE * CPU_eigenvalues = new GPU_GWAS_TYPE[pitch]; 
    memset(CPU_eigenvalues, 0, sizeof(GPU_GWAS_TYPE)*pitch);
    for(int col = 0; col < n_subjects; col++){
     //   CPU_mean_column[col] = (double)eigen_mean_column(col);
    //    CPU_component_A_D[col] = (double)eigen_component_A(col);
        CPU_eigenvalues[col] =(GPU_GWAS_TYPE) U(col, 1);
        for(int row = 0; row < n_subjects; row++){
            CPU_eigenvectors_transposed[col*pitch + row] = (GPU_GWAS_TYPE)eigen_eigenvectors_transposed(row, col);
        }
    }
    Eigen::VectorXd eigen_Y;
    Eigen::VectorXd eigen_component_D;
    GPU_GWAS_TYPE * CPU_Y_component_D = new GPU_GWAS_TYPE[pitch*2];
    memset(CPU_Y_component_D, 0, sizeof(GPU_GWAS_TYPE)*pitch*2);
    GPU_GWAS_Estimator * estimator;	
   // try{
	estimator = new GPU_GWAS_Estimator(gpu_id_list, plink_file, index_map,CPU_mean_column_component_A,\
                                        CPU_eigenvectors_transposed, CPU_eigenvalues, \
                                         precision, \
                                        max_batch_size,  blockSize, \
                                         n_snps, n_subjects, pitch, n_streams, verbose);	
   /* }catch(exception & e){
	string error_message = string(e.what());
	throw GPU_GWAS_Exception(error_message);
    }catch(...){
	string error_message = "Unknown error occurred initializing gpu gwas estimator data";
	throw GPU_GWAS_Exception(error_message);
    }*/
    for(int index = 0; index < current_eigen_data_set->get_n_phenotypes(); index++){
	
   /* if(run_multi_trait){
        const unsigned current_n_phenotypes = current_eigen_data_set->get_n_phenotypes();
        double null_loglik;
        double SD;
        Multi_Trait_GPU_GWAS_Estimator * estimator =  new Multi_Trait_GPU_GWAS_Estimator(gpu_id_list,CPU_eigenvalues,\
                                       CPU_mean_column, \
                                       h2r, loglik, beta, \
                                       beta_se,  variance,\
                                       precision, max_batch_size,  \
                                       dataset_thread_count, blockSize,\
                                       scale, n_snps, n_subjects, pitch);
        estimator->Process_SNP_Data(plink_file, CPU_eigenvectors_transposed, index_map);
       
        for(size_t index = 0 ; index < current_n_phenotypes;index++){
            eigen_Y = eigen_eigenvectors_transposed*Eigen::Map<Eigen::VectorXd>(current_eigen_data_set->get_phenotype_column(index), n_subjects);
            eigen_component_D = eigen_Y.cwiseProduct(eigen_mean_column);
            gwas_data null_result = compute_null_model_MLE(eigen_Y, eigen_mean_column, U, precision);
            
            for(size_t row =0 ;row< n_subjects; row++){
                CPU_component_A_D[pitch + row] = (double)eigen_component_D(row);
                CPU_Y[row] = (double)eigen_Y(row);
            }
            
            estimator->run(CPU_Y,  CPU_component_A_D);
            write_output(current_eigen_data_set->get_trait_name(index) , null_result);
            
        }
        
        delete estimator;
        
    }else{*/
              eigen_Y = Eigen::Map<Eigen::VectorXd>(current_eigen_data_set->get_phenotype_column(index), n_subjects);
	      eigen_Y = eigen_eigenvectors_transposed*eigen_Y;
	      eigen_component_D = eigen_Y.cwiseProduct(eigen_mean_column);
        //eigen_component_D = eigen_Y.cwiseProduct(eigen_mean_column);
              double SD;
       	      gwas_data null_result = gwas_maximize_newton_raphson_method_null_model(eigen_Y, eigen_mean_column,U, precision);
       	      for(int row =0 ;row< n_subjects; row++){
        //    CPU_component_A_D[pitch+ row] = (double)eigen_component_D(row);
                   CPU_Y_component_D[row] = (GPU_GWAS_TYPE)eigen_Y(row);
	           CPU_Y_component_D[row + pitch] = eigen_component_D(row);
      	      }
      // cout << "n_streams " << n_streams << endl;
     //  cout << "pitch " << pitch << endl;
      // cout << "blockSize " << blockSize << endl;
    //   cout << "gpu list " << gpu_id_list.size() << " " << gpu_id_list[0] << endl;
      // cout << "max batch size " << max_batch_size << endl;
		string str_output_filename = current_eigen_data_set->get_trait_name(index) + "-gpu-gwas.out";
		pio_reset_row(plink_file);
		//try{
        	   	estimator->run(current_eigen_data_set->get_trait_name(index), CPU_Y_component_D, str_output_filename.c_str(), null_result);
		/*}catch(exception & e){
		 	string error_message = string(e.what());
			throw GPU_GWAS_Exception(error_message);
    		}catch(...){
			string error_message = "Unknown error occurred during gpu gwas parameter estimation";
			throw GPU_GWAS_Exception(error_message);
    		}*/
        	//write_output(current_eigen_data_set->get_trait_name(index), null_result);
	}
        delete estimator;
   // }
auto finish = std::chrono::high_resolution_clock::now();
   
      auto microseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
        std::cout << "time " << microseconds.count() << "ms\n";
    set_index++;
    delete [] CPU_Y_component_D;
    delete [] CPU_mean_column_component_A;
  //  delete [] CPU_component_A_D;
    delete [] CPU_eigenvectors_transposed;
   // delete [] CPU_mean_column;
    delete [] CPU_eigenvalues;
    return 0;
}

GPU_GWAS::GPU_GWAS(vector<int> _gpu_id_list, vector<string> trait_list, const char * phenotype_filename,const char * evd_data_filename,\
                    const char * plink_filename, int _n_streams, const int _precision, \
		   const int _blockSize, const int _defined_batch_size, const bool _verbose){
    verbose = _verbose;
    gpu_id_list = _gpu_id_list;

    precision = _precision;
    //output_filename = _output_filename;
    plink_file = new pio_file_t;
    defined_blockSize = _blockSize;
    n_streams = _n_streams;
    defined_batch_size = _defined_batch_size;
    pio_status_t plink_status = pio_open(plink_file, plink_filename);
    if(plink_status != PIO_OK){
        string error_message = "Error loading plink file";
        throw GPU_GWAS_Exception(error_message);
    }
    n_snps = plink_file->bed_file.header.num_loci;
    pio_sample_t * sample;
    for(size_t sample_index = 0 ;sample_index < plink_file->bed_file.header.num_samples; sample_index++){
        sample = pio_get_sample(plink_file, sample_index);
        plink_ids.push_back(string(sample->iid));
    }
   
    //if(trait_list.size() == 1){
        run_multi_trait = false;
    //}else{
    //    run_multi_trait = true;
  //  }
    
   // initialize_output_stream(); 
    if(verbose){
    	cout << "Reading phenotype file for traits\n";
    }
    if(evd_data_filename){
	//try{
	 	reader = new Solar_Trait_Reader(phenotype_filename,evd_data_filename, trait_list);
	/*}catch(exception & e){
		string error_message = string(e.what());
		throw GPU_GWAS_Exception(error_message);
	}catch(...){
		string error_message = "Unknown error occurred reading phenotype or pedigree data";
		throw GPU_GWAS_Exception(error_message);
	}*/
     }else{
	//try{
		reader = new Solar_Trait_Reader(phenotype_filename, trait_list, plink_ids);
	/*}catch(exception & e){
		string error_message = string(e.what());
		throw GPU_GWAS_Exception(error_message);
	}catch(...){
		string error_message = "Unknown error occurred reading phenotype or pedigree data";
		throw GPU_GWAS_Exception(error_message);
	}*/
    }
    n_sets = reader->get_n_sets();
    if(n_sets >= 1){
   	set_index = 0;
   	//loglik = new GPU_GWAS_TYPE[n_snps];
    	
    }

}

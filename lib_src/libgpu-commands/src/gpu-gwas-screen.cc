#include "gpu-gwas-screen-variables.h"
#include "gpu-gwas-estimator-context.h"
#include "solar-trait-reader.h"
#include "gpu-exception.h"
#include "Eigen/Dense"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "plinkio.h"
using namespace std;
static unsigned  N_SNPS_COMPUTED;
static unsigned LAST_AMOUNT;
extern double gwas_chicdf(double, double);
void gpu_gwas_screen_kernels(GPU_Data * gpu_data, GPU_GWAS_Screen_Data * gpu_screen_data, GPU_GWAS_Estimator_SNP_Data * snp_data);
Eigen::VectorXd calculate_theta(Eigen::VectorXd Y, Eigen::MatrixXd aux);
static void initialize_output_stream(ofstream & output_stream,string trait_name){
	string output_filename = trait_name + "-gpu-gwas-screen.out";
	output_stream.open(output_filename.c_str());
	output_stream << "SNP,p-value,beta,beta_se\n";
}
static void write_data_to_file(ofstream & output_stream, pio_file_t * plink_file, float * chi_squared, float * beta_data, const int start_index, const int batch_size){
	pio_locus_t * locus;
	
	for(int index = 0; index < batch_size; index++){
		
		locus = pio_get_locus(plink_file, start_index + index);
        	const double chi = chi_squared[index];
       		if(chi > 0 && chi == chi){
            		const double pvalue = gwas_chicdf(chi, 1.0);
 			output_stream  << locus->name << "," << pvalue << "," << beta_data[index*2] << "," << beta_data[index*2 + 1] << endl;
        	}else{
            		output_stream << locus->name << "," << 1 << "," << 0 << "," << 0 << endl;
        	}	
	}	

}
static void Thread_Launch(const int gpu_index, const int stream_index, GPU_Data * gpu_data, GPU_GWAS_Screen_Data * screen_data, GPU_GWAS_Estimator_SNP_Data * snp_data, ofstream & output_stream, const bool verbose){
    cudaErrorCheck(cudaSetDevice(gpu_data->gpu_id));
    	const int n_snps= pio_num_loci(snp_data->plink_file);
	int snp_batch_size;
	int last_snp_index = 0;
	int last_snp_batch_size = 0;
	int current_start;
	do{
		snp_batch_size = snp_data->prepare_snp_data();
		if(snp_batch_size != 0){	
        		snp_data->copy_SNP_data_to_gpu();
        		current_start = snp_data->get_next_start();

        		gpu_gwas_screen_kernels(gpu_data, screen_data, snp_data);
        		if(last_snp_batch_size != 0){
		#pragma omp critical
				{	
			 		write_data_to_file(output_stream, snp_data->plink_file, screen_data->cpu_chi_squared, screen_data->cpu_beta_data, last_snp_index, last_snp_batch_size);	
				}
			}
        		screen_data->copy_results_to_cpu(snp_batch_size, gpu_data->stream);
			last_snp_index = current_start;
			last_snp_batch_size = snp_batch_size;
			cudaErrorCheck(cudaStreamSynchronize(gpu_data->stream));
		}
		
          	if(verbose){
		#pragma omp critical
		    { 
			N_SNPS_COMPUTED += snp_batch_size;
			if(N_SNPS_COMPUTED - LAST_AMOUNT == snp_batch_size){	
				std::cout << "SNPs Computed: " << N_SNPS_COMPUTED << " SNPs Left: " << n_snps - N_SNPS_COMPUTED << " Percent Complete " << 100*N_SNPS_COMPUTED/n_snps << "% \r";
				std::cout.flush();
		 		LAST_AMOUNT = N_SNPS_COMPUTED;
			}
	   	  }
		}
	}while(snp_batch_size != 0);
        if(last_snp_batch_size != 0){
	#pragma omp critical
		{
			write_data_to_file(output_stream, snp_data->plink_file, screen_data->cpu_chi_squared, screen_data->cpu_beta_data, last_snp_index, last_snp_batch_size);		
		}
	}	
  	 std::cout << std::endl;	
	

}

extern void determine_GPU_GWAS_sizes(bool multi,vector<int> gpu_id_list, int n_subjects, int n_snps, \
                              int & max_batch_size, int & blockSize,  int & pitch, \
                              const int n_streams, const int defined_blockSize, const int defined_batch_size);
const char * gpu_gwas_screen(vector<int> gpu_ids, vector<string> trait_list, const char*  phenotype_filename, const char* plink_filename, const char * evd_data_filename, int defined_max_batch_size, int defined_thread_size,\
				 int n_streams, const bool verbose){
	pio_file_t * plink_file = new pio_file_t;
	
	pio_open(plink_file, plink_filename);
	const int n_snps = plink_file->bed_file.header.num_loci;
	cout << "n_snps: " << n_snps << endl;	
	vector<string> plink_ids;
	pio_sample_t * sample; 
	for(unsigned id = 0; id < plink_file->bed_file.header.num_samples; id++){
		sample = pio_get_sample(plink_file, id);
		string iid = string(sample->iid);
		plink_ids.push_back(iid);
	}

	Solar_Trait_Reader * file_reader;
	if(!evd_data_filename){
		//try{
	 		file_reader = new Solar_Trait_Reader(phenotype_filename, trait_list,  plink_ids);
	 	//}catch(Solar_Trait_Reader_Exception & e){
	 		//return e.what();
	 	//}
	 }else{
		//try{
	 		file_reader = new Solar_Trait_Reader(phenotype_filename,evd_data_filename, trait_list);
	 	//}catch(Solar_Trait_Reader_Exception & e){
	 		//return e.what();
	 	//}	 	
	 }
	for(unsigned set = 0; set < file_reader->get_n_sets(); set++){
		Eigen_Data * eigen_data = file_reader->get_eigen_data_set(set);
		int n_subjects = eigen_data->get_n_subjects();
		int blockSize;
		int pitch;
		int max_batch_size;
		
		determine_GPU_GWAS_sizes(true,gpu_ids, n_subjects,  n_snps, \
                              max_batch_size,  blockSize,  pitch, \
                               n_streams, defined_thread_size, defined_max_batch_size);	
                
                 cout << "n_subjects: " << n_subjects << " pitch " << pitch <<  "max_batch_size  " << max_batch_size  << "blockSize " << blockSize << endl;        
          	int * index_map = new int[plink_ids.size()];
    		vector<string> ids = eigen_data->get_ids();
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
    (eigen_data->get_eigenvectors_transposed(), n_subjects, n_subjects );
    		float * cpu_eigenvectors = new float[pitch*pitch];
    		memset(cpu_eigenvectors, 0, sizeof(float)*pitch*pitch);
    		for(int col = 0; col < n_subjects; col++){
    			for(int row= 0;row < n_subjects; row++){
    				cpu_eigenvectors[col*pitch + row] = eigen_eigenvectors_transposed(row, col);
    			}
    		}
		Eigen::VectorXd eigenvalues =  Eigen::Map<Eigen::VectorXd>(eigen_data->get_eigenvalues(), n_subjects );
		Eigen::MatrixXd U = Eigen::MatrixXd::Ones(eigenvalues.rows(), 2);
		U.col(1) = eigenvalues;
		float ** gpu_eigenvectors = new float*[gpu_ids.size()];
		for(int gpu = 0; gpu < gpu_ids.size(); gpu++){
			cudaErrorCheck(cudaSetDevice(gpu_ids[gpu]));
			cudaErrorCheck(cudaMalloc((void**)&gpu_eigenvectors[gpu], sizeof(float)*pitch*pitch));
			cudaErrorCheck(cudaMemcpy(gpu_eigenvectors[gpu], cpu_eigenvectors, sizeof(float)*pitch*pitch, cudaMemcpyHostToDevice)); 
		}
		for(int trait = 0 ; trait < eigen_data->get_n_phenotypes(); trait++){
			pio_reset_row(plink_file);
			string trait_name = eigen_data->get_trait_name(trait);
			Eigen::VectorXd trait_vector = Eigen::Map<Eigen::VectorXd>(eigen_data->get_phenotype_column(trait), n_subjects);
			Eigen::VectorXd Y = eigen_eigenvectors_transposed*(trait_vector.array() - trait_vector.mean()).matrix();
			Eigen::VectorXd theta = calculate_theta(Y, U);
			if(theta(0) != theta(0) || theta(1) != theta(1)){
				cout << "Cannot run GWAS Screen on Trait " << eigen_data->get_trait_name(trait) << " since variance components could not be calculated\n";
				continue;
			}
			Eigen::VectorXd eigen_Sigma = (U*theta).cwiseInverse();
			float * Sigma_Y = new float[pitch];
			float * Sigma = new float [pitch];
			memset(Sigma_Y, 0.f, sizeof(float)*pitch);
			memset(Sigma, 0.f, sizeof(float)*pitch);
			for(int row = 0; row < n_subjects; row++){
				Sigma_Y[row] = Y(row)*eigen_Sigma(row);
				Sigma[row] = eigen_Sigma(row);
			}
			float ** gpu_Sigma_Y = new float*[gpu_ids.size()];
			float ** gpu_Sigma = new float*[gpu_ids.size()];
			
			GPU_Data ** gpu_data = new GPU_Data*[gpu_ids.size()*n_streams];
			GPU_GWAS_Estimator_SNP_Data ** snp_data = new GPU_GWAS_Estimator_SNP_Data*[gpu_ids.size()*n_streams];
			GPU_GWAS_Screen_Data ** screen_data = new GPU_GWAS_Screen_Data*[gpu_ids.size()*n_streams];
			for(int gpu = 0; gpu < gpu_ids.size(); gpu++){
				cudaErrorCheck(cudaSetDevice(gpu_ids[gpu]));
				cudaErrorCheck(cudaMalloc((void**)&gpu_Sigma_Y[gpu], sizeof(float)*pitch));
				cudaErrorCheck(cudaMalloc((void**)&gpu_Sigma[gpu], sizeof(float)*pitch));
				
				cudaErrorCheck(cudaMemcpy(gpu_Sigma_Y[gpu], Sigma_Y, sizeof(float)*pitch,cudaMemcpyHostToDevice)); 
				cudaErrorCheck(cudaMemcpy(gpu_Sigma[gpu], Sigma, sizeof(float)*pitch,cudaMemcpyHostToDevice));
				 
				for(int stream = 0; stream < n_streams; stream++){
					gpu_data[gpu*n_streams + stream] = new GPU_Data(gpu_ids[gpu], blockSize, max_batch_size);
					snp_data[gpu*n_streams + stream] = new GPU_GWAS_Estimator_SNP_Data(gpu_eigenvectors[gpu], plink_file,gpu_data[gpu*n_streams + stream], \
                                blockSize, n_subjects, n_snps, max_batch_size, pitch,  index_map);
                                	screen_data[gpu*n_streams + stream] = new GPU_GWAS_Screen_Data(gpu_Sigma_Y[gpu], gpu_Sigma[gpu],   n_subjects , max_batch_size , pitch);
                                }
                         }
                        
                         ofstream output_stream;
                         N_SNPS_COMPUTED = 0;
                         initialize_output_stream(output_stream, trait_name);
    			if(gpu_ids.size() == 1 && n_streams == 1){
       		 		Thread_Launch(0, 0, gpu_data[0], screen_data[0], snp_data[0], output_stream, verbose);
  			 }else{
        			omp_set_num_threads(gpu_ids.size()*n_streams);
			#pragma omp parallel
        			{
            				const int total_index  = omp_get_thread_num();
	    				const int gpu_index = floor((double)total_index/n_streams);
	    				const int stream_index = total_index % n_streams;
	    				//try{
            					Thread_Launch(gpu_index,stream_index, gpu_data[gpu_index*n_streams + stream_index], screen_data[gpu_index*n_streams + stream_index], snp_data[gpu_index*n_streams + stream_index], output_stream,verbose);
            				//}catch(std::exception & e){
            				//	ex = e;
            				//	exception_flagged = true;
            				//}
       				 }
       			}
       			output_stream.close();
       			for(int gpu = 0; gpu < gpu_ids.size(); gpu++){
       				cudaErrorCheck(cudaSetDevice(gpu_ids[gpu]));
				
       				for(int stream = 0; stream < n_streams; stream++){
       					delete gpu_data[gpu*n_streams + stream];
       					delete snp_data[gpu*n_streams + stream];
       					delete screen_data[gpu*n_streams + stream];
       				}
       			       	cudaErrorCheck(cudaFree(gpu_Sigma_Y[gpu]));
       				cudaErrorCheck(cudaFree(gpu_Sigma[gpu]));	
       			}
       			delete [] gpu_data;
       			delete [] snp_data;
       			delete [] screen_data;
       			delete [] gpu_Sigma_Y;
       			delete [] gpu_Sigma;
       			delete [] Sigma_Y;
       			delete [] Sigma;
       				                   
                         
			
		}
		for(int gpu = 0; gpu < gpu_ids.size(); gpu++){
			cudaErrorCheck(cudaSetDevice(gpu_ids[gpu]));
			cudaErrorCheck(cudaFree(gpu_eigenvectors[gpu]));
		}
		delete [] gpu_eigenvectors;		
		delete [] cpu_eigenvectors;
		delete [] index_map;
	}
	delete file_reader;
	pio_close(plink_file);
	delete plink_file;

	return 0;
}


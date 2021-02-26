//
//  gpu-gwas-estimator.cc
//  
//
//  Created by Brian Donohue on 8/10/18.
//

#include <stdio.h>
#include "gpu-gwas-estimator.h"
#include <iostream>
#include <omp.h>
#include <exception>
#include <string>
#include <iomanip>
//static int total_snps_computed;
static unsigned  N_SNPS_COMPUTED;
static unsigned LAST_AMOUNT;
extern double gwas_chicdf(double, double);
void GPU_GWAS_Estimator::write_data_to_file(GPU_GWAS_TYPE * parameters, int snp_index, int snp_batch_size){
	pio_locus_t * locus;
	for(int index = 0; index < snp_batch_size; index++){
		locus = pio_get_locus(plink_file, snp_index + index);
        	double chi2 = 2.0*(parameters[index*5 + 4] - null_result.loglik);
        	double pvalue;
       		if(chi2 >= 0 && parameters[5*index] == parameters[5*index] && parameters[5*index] != 99){
            		pvalue = gwas_chicdf(chi2, 1.0);
 		output_stream  << locus->name << "," << parameters[index*5] << "," << parameters[index*5 + 4] << ",";
        	output_stream <<  sqrt(parameters[index*5 +3]) << "," << parameters[index*5 + 1] << "," << parameters[index*5 + 2] << "," << chi2 << "," << pvalue << ",Success\n";
       		 }else if(parameters[5*index] == 99){
            	output_stream << locus->name << "," << null_result.h2r << "," << null_result.loglik << "," << null_result.SD << "," << null_result.beta << "," << null_result.SE << "," << null_result.chi << "," << null_result.pvalue << ",Iteration Limit Reached\n";
        	}else{
            	output_stream << locus->name << "," << null_result.h2r << "," << null_result.loglik << "," << null_result.SD << "," << null_result.beta << "," << null_result.SE << "," << null_result.chi << "," << null_result.pvalue << ",Failure\n";
        	}	
	}		

}
void GPU_GWAS_Estimator::initialize_output_stream(const char * filename){
	 output_stream.open(filename);
	 output_stream << "SNP,h2r,loglik,SD,beta_snp,beta_snp_se,chi2,p-value,Status\n";
}
void GPU_GWAS_Estimator::close_output_stream(){
	output_stream.close();
}
void GPU_GWAS_Estimator::Thread_Launch(const int gpu_index, const int stream_index){
    cudaErrorCheck(cudaSetDevice(gpu_data[gpu_index*n_streams + stream_index]->gpu_id));
	int snp_batch_size;
	int last_snp_index = 0;
	int last_snp_batch_size = 0;
	int output_str_largest_length = 0;
	int current_start;
	do{
		snp_batch_size = snp_data[gpu_index*n_streams + stream_index]->prepare_snp_data();
		if(snp_batch_size != 0){	
        		snp_data[gpu_index*n_streams + stream_index]->copy_SNP_data_to_gpu();
        		current_start = snp_data[gpu_index*n_streams + stream_index]->get_next_start();
        		Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer(gpu_data[gpu_index*n_streams + stream_index], snp_data[gpu_index*n_streams + stream_index]);
       			Run_GPU_GWAS_Function_Pointer(gpu_data[gpu_index*n_streams + stream_index], context[gpu_index*n_streams + stream_index], results[gpu_index*n_streams\
					 + stream_index], snp_data[gpu_index*n_streams + stream_index]->SNP_data);
        		if(last_snp_batch_size != 0 && open_output_stream){
		#pragma omp critical
				{	
			 		write_data_to_file(results[gpu_index*n_streams + stream_index]->cpu_parameters, last_snp_index, last_snp_batch_size);	
				}
			}
        		results[gpu_index*n_streams + stream_index]->copy_results_to_cpu(snp_batch_size);
			last_snp_index = current_start;
			last_snp_batch_size = snp_batch_size;
			cudaErrorCheck(cudaStreamSynchronize(gpu_data[gpu_index*n_streams + stream_index]->stream));
		}
		
          	if(verbose){
		#pragma omp critical
		    { 
			N_SNPS_COMPUTED += snp_batch_size;
			if(!gpu_index && !stream_index){
				std::string output_str =  "Trait: " + trait_name + " SNPs Left: " + std::to_string(n_snps - N_SNPS_COMPUTED) + " Percent Complete " + std::to_string(floor(100.f*N_SNPS_COMPUTED/n_snps)) + "%\r";
				if(output_str_largest_length == 0) output_str_largest_length = output_str.length();
				if(output_str_largest_length <= output_str.length()){
					output_str_largest_length = output_str.length();
					std::cout << output_str << std::flush;
				}else if (output_str_largest_length > output_str.length()){
					output_str =  "Trait: " + trait_name + " SNPs Left: " + std::to_string(n_snps - N_SNPS_COMPUTED) + " Percent Complete " + std::to_string(floor(100.f*N_SNPS_COMPUTED/n_snps)) + "%";
					std::cout << output_str << std::setw(output_str_largest_length - output_str.length()) << "\r" << std::flush;	
				}
				//std::cout << "Trait: " << trait_name << " SNPs Computed: " << N_SNPS_COMPUTED << " SNPs Left: " << n_snps - N_SNPS_COMPUTED << " Percent Complete " << floor(100.0*N_SNPS_COMPUTED/n_snps) << "% \r";
				//std::cout.flush();
				//std::cout << output_str << std::flush;
		 		LAST_AMOUNT = N_SNPS_COMPUTED;
			}
	   	  }
		}
	}while(snp_batch_size != 0);
        if(last_snp_batch_size != 0 && open_output_stream){
	#pragma omp critical
		{
			write_data_to_file(results[gpu_index*n_streams + stream_index]->cpu_parameters, last_snp_index, last_snp_batch_size);		
		}
	}	

	if(verbose && !gpu_index && !stream_index){
 
		const int str_length = 59 + 2*std::to_string(n_snps).length() + 5;
		for (int i = 0 ; i < str_length; i++){
			std::cout <<  " ";
		}
		std::cout << " \n";
	}
}

void GPU_GWAS_Estimator::run(std::string _trait_name, const GPU_GWAS_TYPE * const CPU_Y_component_D, const char * filename, gwas_data _null_result){
	trait_name = _trait_name;
	open_output_stream = true;
	if(!strcmp(filename,"calibrate")) open_output_stream = false;
		
	null_result = _null_result;
	
    for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
	cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
	cudaErrorCheck(cudaMemcpy(GPU_Y_component_D_mean_column_component_A[gpu_index], CPU_Y_component_D, sizeof(GPU_GWAS_TYPE)*pitch*2, cudaMemcpyHostToDevice));
    }
    
    N_SNPS_COMPUTED = 0;
    LAST_AMOUNT = 0;
    if(open_output_stream) initialize_output_stream(filename);
    //total_snps_computed = 0;
    std::exception ex;
    bool exception_flagged = false;
    if(gpu_id_list.size() == 1 && n_streams == 1){
        Thread_Launch(0, 0);
    }else{
    	omp_set_dynamic(0);
        omp_set_num_threads(gpu_id_list.size()*n_streams);
#pragma omp parallel
        {
            const int total_index  = omp_get_thread_num();
	    const int gpu_index = floor((double)total_index/n_streams);
	    const int stream_index = total_index % n_streams;
	    try{
            	Thread_Launch(gpu_index,stream_index);
            }catch(std::exception & e){
            	ex = e;
            	exception_flagged = true;
            }
        }
        if(exception_flagged){
        	throw ex;
        }
    }
    close_output_stream();
    
}

GPU_GWAS_Estimator::~GPU_GWAS_Estimator(){
 
    for(int gpu_index =0 ; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
	for(int stream_index = 0; stream_index < n_streams; stream_index++){
      		 delete context[gpu_index*n_streams + stream_index];
       		 delete results[gpu_index*n_streams + stream_index];
       		 delete gpu_data[gpu_index*n_streams + stream_index];
        	 delete snp_data[gpu_index*n_streams + stream_index];
	}
    }
    free_gpu_memory();
    
    delete [] gpu_data;
    delete [] context;
    delete [] results;
    delete [] snp_data;
    
    delete [] GPU_eigenvectors_transposed;
    delete [] GPU_Y_component_D_mean_column_component_A;
    delete [] GPU_eigenvalues;
   /* for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
	cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
	cudaErrorCheck(cudaDeviceReset());
    } */  
}

void GPU_GWAS_Estimator::free_gpu_memory(){
    for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        
        cudaErrorCheck(cudaFree(GPU_eigenvectors_transposed[gpu_index]));
        
        cudaErrorCheck(cudaFree(GPU_Y_component_D_mean_column_component_A[gpu_index]));
        
       // cudaErrorCheck(cudaFree(GPU_mean_column[gpu_index]));
        
      //  cudaErrorCheck(cudaFree(GPU_raw_components_A_D[gpu_index]));
        
       
        
        cudaErrorCheck(cudaFree(GPU_eigenvalues[gpu_index]));
        
    }
}

void GPU_GWAS_Estimator::allocate_gpu_memory(){
    for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_eigenvectors_transposed[gpu_index], sizeof(GPU_GWAS_TYPE)*pitch*pitch));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_Y_component_D_mean_column_component_A[gpu_index], sizeof(GPU_GWAS_TYPE)*pitch*4));
        
        //cudaErrorCheck(cudaMalloc((void**)&GPU_mean_column[gpu_index], sizeof(double)*pitch));
        
      //  cudaErrorCheck(cudaMalloc((void**)&GPU_raw_components_A_D[gpu_index], sizeof(double)*pitch*2));
        
        
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_eigenvalues[gpu_index], sizeof(GPU_GWAS_TYPE)*pitch));
        
    }
}

GPU_GWAS_Estimator::GPU_GWAS_Estimator(std::vector<int> _gpu_id_list, pio_file_t * const _plink_file, const int * const index_map,const GPU_GWAS_TYPE * const CPU_mean_column_component_A,\
                                                                const GPU_GWAS_TYPE * const CPU_eigenvectors_transposed, const GPU_GWAS_TYPE * const  CPU_eigenvalues,\
							        int  _precision, \
                                                               int _max_batch_size, int _blockSize, int total_snps, int _n_subjects,\
						               int _pitch,int _n_streams, const bool _verbose){
    gpu_id_list = _gpu_id_list;
    verbose = _verbose;    
    plink_file = _plink_file;
   // loglik = _loglik;
 
    n_streams = _n_streams;
    
    precision = _precision;
    n_snps = total_snps;
    max_batch_size  = _max_batch_size;
    n_subjects = _n_subjects;
    pitch = _pitch;
    
    blockSize = _blockSize;

    
    snp_data = new GPU_GWAS_Estimator_SNP_Data*[gpu_id_list.size()*n_streams];
    context = new GPU_GWAS_Estimator_Context*[gpu_id_list.size()*n_streams];
    results = new GPU_GWAS_Estimator_Data*[gpu_id_list.size()*n_streams];
    gpu_data = new GPU_Data*[gpu_id_list.size()*n_streams];
    
    GPU_Y_component_D_mean_column_component_A = new GPU_GWAS_TYPE*[gpu_id_list.size()];
   // GPU_mean_column = new double*[gpu_id_list.size()];
  //  GPU_raw_components_A_D = new double*[gpu_id_list.size()];
    GPU_eigenvalues = new GPU_GWAS_TYPE*[gpu_id_list.size()];
    GPU_eigenvectors_transposed = new GPU_GWAS_TYPE*[gpu_id_list.size()];

    allocate_gpu_memory();
   
    for(int gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
       // cudaErrorCheck(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
      //  cudaErrorCheck(cudaSetDeviceFlags(cudaDeviceScheduleYield));
	//cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
       // cudaErrorCheck(cudaMemcpy(GPU_Y[gpu_index], CPU_Y, sizeof(GPU_GWAS_TYPE)*pitch, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpy(GPU_eigenvectors_transposed[gpu_index], CPU_eigenvectors_transposed, sizeof(GPU_GWAS_TYPE)*pitch*pitch, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpy(GPU_eigenvalues[gpu_index], CPU_eigenvalues, sizeof(GPU_GWAS_TYPE)*pitch, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpy(GPU_Y_component_D_mean_column_component_A[gpu_index] + pitch*2, CPU_mean_column_component_A, sizeof(GPU_GWAS_TYPE)*pitch*2, cudaMemcpyHostToDevice));
	for(int stream_index = 0; stream_index < n_streams; stream_index++){
        	gpu_data[gpu_index*n_streams + stream_index] = new GPU_Data(gpu_id_list[gpu_index], blockSize, max_batch_size);
        

        
       
        
        
      //  cudaErrorCheck(cudaMemcpy(GPU_mean_column[gpu_index], CPU_mean_column, sizeof(GPU_GWAS_TYPE)*n_subjects, cudaMemcpyHostToDevice));
      //  cudaErrorCheck(cudaMemcpy(GPU_raw_components_A_D[gpu_index], CPU_component_A_D, sizeof(double)*pitch*2, cudaMemcpyHostToDevice));
        
        	context[gpu_index*n_streams + stream_index] = new GPU_GWAS_Estimator_Context(gpu_data[gpu_index*n_streams + stream_index], GPU_Y_component_D_mean_column_component_A[gpu_index], GPU_eigenvalues[gpu_index],\
                                                          n_subjects, max_batch_size,\
                                                          blockSize, pitch, precision);
        
        	snp_data[gpu_index*n_streams + stream_index] = new GPU_GWAS_Estimator_SNP_Data(GPU_eigenvectors_transposed[gpu_index], plink_file,gpu_data[gpu_index*n_streams + stream_index], \
                                                            blockSize, n_subjects, n_snps, max_batch_size, pitch,index_map);
        	results[gpu_index*n_streams + stream_index] = new GPU_GWAS_Estimator_Data(gpu_data[gpu_index*n_streams + stream_index],max_batch_size);
	}
    }
    
    
    
    
    set_function_pointers(blockSize);
    
}
/*

Multi_Trait_GPU_GWAS_Estimator::Multi_Trait_GPU_GWAS_Estimator(std::vector<int> _gpu_id_list, \
                               const GPU_GWAS_TYPE * const  CPU_eigenvalues,GPU_GWAS_TYPE * _h2r, GPU_GWAS_TYPE * _loglik, GPU_GWAS_TYPE * _beta, \
                               GPU_GWAS_TYPE * _beta_se,  GPU_GWAS_TYPE * _variance, int _precision, unsigned _max_batch_size,\
			       unsigned _blockSize, unsigned total_snps, unsigned _n_subjects, unsigned _pitch, unsigned _n_streams){
    gpu_id_list = _gpu_id_list;
    
    h2r = _h2r;
    loglik = _loglik;
    variance = _variance;
    beta = _beta;
    beta_se = _beta_se;
    
    
    precision = _precision;
    n_snps = total_snps;
    max_batch_size  = _max_batch_size;
    n_subjects = _n_subjects;
    pitch = _pitch;
    n_streams = _streams;
    blockSize = _blockSize;
  //  scale = _scale;
    
    context = new GPU_GWAS_Estimator_Context*[gpu_id_list.size()*n_streamc];
    results = new GPU_GWAS_Estimator_Data*[gpu_id_list.size()*n_streams];
    gpu_data = new GPU_Data*[gpu_id_list.size()*n_streams];
    
    GPU_Y = new GPU_GWAS_TYPE*[gpu_id_list.size()];
   // GPU_mean_column = new double*[gpu_id_list.size()];
   // GPU_raw_components_A_D = new double*[gpu_id_list.size()];
    GPU_eigenvalues = new GPU_GWAS_TYPE*[gpu_id_list.size()];
    GPU_SNP_Data = new GPU_GWAS_TYPE*[gpu_id_list.size()*n_streams];
    
    CPU_SNP_Data = new GPU_GWAS_TYPE[pitch*n_snps];
    
    CPU_Pinned_SNP_Data = new GPU_GWAS_TYPE*[gpu_id_list.size()];
    allocate_gpu_memory();
    
    for(size_t gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        gpu_data[gpu_index] = new GPU_Data(gpu_id_list[gpu_index], blockSize, max_batch_size);
        cudaErrorCheck(cudaMallocHost((void**)&CPU_Pinned_SNP_Data[gpu_index],sizeof(GPU_GWAS_TYPE)*pitch*max_batch_size,cudaHostAllocWriteCombined));
        
        cudaErrorCheck(cudaMemcpy(GPU_eigenvalues[gpu_index], CPU_eigenvalues, sizeof(GPU_GWAS_TYPE)*pitch, cudaMemcpyHostToDevice));
       // cudaErrorCheck(cudaMemcpy(GPU_mean_column[gpu_index], CPU_mean_column, sizeof(GPU_GWAS_TYPE)*n_subjects, cudaMemcpyHostToDevice));
        
        context[gpu_index] = new GPU_GWAS_Estimator_Context(gpu_data[gpu_index], GPU_Y[gpu_index], GPU_eigenvalues[gpu_index],\
                                                          n_subjects, max_batch_size, blockSize,  pitch);
        
        
        results[gpu_index] = new GPU_GWAS_Estimator_Data(gpu_data[gpu_index], h2r, loglik, beta, beta_se, variance, max_batch_size);
    }
    
    
    
    
    set_function_pointer(blockSize);
    
}

Multi_Trait_GPU_GWAS_Estimator::~Multi_Trait_GPU_GWAS_Estimator(){
    delete [] CPU_SNP_Data;

    for(size_t gpu_index =0 ; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        cudaErrorCheck(cudaFreeHost(CPU_Pinned_SNP_Data[gpu_index]));
        delete context[gpu_index];
        delete results[gpu_index];
        delete gpu_data[gpu_index];
    }
    free_gpu_memory();
    delete [] gpu_data;
    delete [] CPU_Pinned_SNP_Data;
    
    delete [] context;
    delete [] results;
    
    delete [] GPU_SNP_Data;
    delete [] GPU_Y;
   // delete [] GPU_mean_column;
    delete [] GPU_eigenvalues;
   // delete [] GPU_raw_components_A_D;
     for(unsigned gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        cudaErrorCheck(cudaDeviceReset());
    }
}

void Multi_Trait_GPU_GWAS_Estimator::Thread_Launch(const int gpu_index){
    cudaErrorCheck(cudaSetDevice(gpu_data[gpu_index]->gpu_id));
    unsigned snp_batch_size = 0;
    unsigned snp_index;
#pragma omp critical
    {
        if(n_snps_computed < n_snps){
            snp_index = n_snps_computed;
            if(n_snps_computed + max_batch_size <= n_snps){
                snp_batch_size = max_batch_size;
                n_snps_computed += snp_batch_size;
            }else{
                snp_batch_size = n_snps - n_snps_computed;
                n_snps_computed = n_snps;
            }
        }else{
	    snp_batch_size = 0;
	}
    }
    while(snp_batch_size != 0){
        results[gpu_index]->pin_cpu_results(snp_index, snp_batch_size);
        memcpy(CPU_Pinned_SNP_Data[gpu_index], CPU_SNP_Data + snp_index*pitch, sizeof(GPU_GWAS_TYPE)*pitch*snp_batch_size);
        cudaErrorCheck(cudaMemcpyAsync(GPU_SNP_Data[gpu_index], CPU_Pinned_SNP_Data[gpu_index], \
                        sizeof(GPU_GWAS_TYPE)*pitch*snp_batch_size, cudaMemcpyHostToDevice, gpu_data[gpu_index]->gpu_stream));
        
        Run_GPU_GWAS_Function_Pointer(gpu_data[gpu_index], context[gpu_index], results[gpu_index], GPU_SNP_Data[gpu_index]);
        
        results[gpu_index]->copy_results_to_cpu(snp_batch_size);
#pragma omp critical
        {
            if(n_snps_computed < n_snps){
                snp_index = n_snps_computed;
                if(n_snps_computed + max_batch_size <= n_snps){
                    snp_batch_size = max_batch_size;
                    n_snps_computed += snp_batch_size;
                }else{
                    snp_batch_size = n_snps - n_snps_computed;
                    n_snps_computed = n_snps;
                }
            }else{
		snp_batch_size = 0;
	    }
        }
        cudaErrorCheck(cudaStreamSynchronize(gpu_data[gpu_index]->gpu_stream));
    }
    if(results[gpu_index]->is_memory_pinned) results[gpu_index]->unpin_cpu_results();
    
}
void Multi_Trait_GPU_GWAS_Estimator::run(const double * const CPU_Y, const double * const CPU_component_A_D){
    
    for(size_t gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        cudaErrorCheck(cudaMemcpy(GPU_Y[gpu_index], CPU_Y, sizeof(double)*n_subjects, cudaMemcpyHostToDevice));
        cudaErrorCheck(cudaMemcpy(GPU_raw_components_A_D[gpu_index], CPU_component_A_D, sizeof(double)*pitch*2, cudaMemcpyHostToDevice));
    }
    
    n_snps_computed = 0;
    
    if(gpu_id_list.size() == 1){
        Thread_Launch(0);
    }else{
        omp_set_num_threads(gpu_id_list.size());
#pragma omp parallel
        {
            const int gpu_index  = omp_get_thread_num();
            Thread_Launch(gpu_index);
        }
    }
    
    
}
void Multi_Trait_GPU_GWAS_Estimator::free_gpu_memory(){
    for(size_t gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        
        cudaErrorCheck(cudaFree(GPU_SNP_Data[gpu_index]));
        
        cudaErrorCheck(cudaFree(GPU_Y[gpu_index]));
        
        cudaErrorCheck(cudaFree(GPU_mean_column[gpu_index]));
        
        cudaErrorCheck(cudaFree(GPU_raw_components_A_D[gpu_index]));
        
    
        
        cudaErrorCheck(cudaFree(GPU_eigenvalues[gpu_index]));
        
    }
}
void Multi_Trait_GPU_GWAS_Estimator::allocate_gpu_memory(){
    for(size_t gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_SNP_Data[gpu_index], sizeof(double)*pitch*max_batch_size));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_Y[gpu_index], sizeof(double)*n_subjects));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_mean_column[gpu_index], sizeof(double)*n_subjects));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_raw_components_A_D[gpu_index], sizeof(double)*pitch*2));
        
        
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_eigenvalues[gpu_index], sizeof(double)*n_subjects));
        
    }
}

void Multi_Trait_GPU_GWAS_Estimator::Process_SNP_Data_Thread_Launch(const size_t gpu_index, GPU_GWAS_Estimator_SNP_Data * const snp_data){
    cudaErrorCheck(cudaSetDevice(gpu_data[gpu_index]->gpu_id));
    size_t snp_batch_size;
#pragma omp critical
    {
        snp_batch_size = snp_data->prepare_snp_data();
    }
    size_t current_start;
    while(snp_batch_size != 0){
        snp_data->copy_SNP_data_to_gpu();
        current_start = snp_data->get_start();
        Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer(snp_data);
        cudaErrorCheck(cudaMemcpyAsync(CPU_SNP_Data + current_start*pitch, snp_data->SNP_data, sizeof(double)*pitch*snp_batch_size, cudaMemcpyDeviceToHost, gpu_data[gpu_index]->gpu_stream));
        
#pragma omp critical
        {
            snp_batch_size = snp_data->prepare_snp_data();
        }
        cudaErrorCheck(cudaStreamSynchronize(gpu_data[gpu_index]->gpu_stream));
    }
}
void Multi_Trait_GPU_GWAS_Estimator::Process_SNP_Data(pio_file_t * plink_file, const double * const CPU_eigenvectors_transposed, const int * const index_map){
    double ** GPU_eigenvectors_transposed = new double*[gpu_id_list.size()];
    GPU_GWAS_Estimator_SNP_Data ** snp_data = new GPU_GWAS_Estimator_SNP_Data*[gpu_id_list.size()];
    for(size_t gpu_index = 0 ; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        
        cudaErrorCheck(cudaMalloc((void**)&GPU_eigenvectors_transposed[gpu_index], sizeof(double)*pitch*pitch));
        cudaErrorCheck(cudaMemcpy(GPU_eigenvectors_transposed[gpu_index], CPU_eigenvectors_transposed, sizeof(double)*pitch*pitch, cudaMemcpyHostToDevice));
        snp_data[gpu_index] = new GPU_GWAS_Estimator_SNP_Data(GPU_eigenvectors_transposed[gpu_index], plink_file,gpu_data[gpu_index], \
                                                            blockSize, scale, n_subjects, max_batch_size, pitch,index_map);
    }
    
    if(gpu_id_list.size() == 1){
        Process_SNP_Data_Thread_Launch(0,snp_data[0]);
    }else{
        omp_set_num_threads(gpu_id_list.size());
#pragma omp parallel
        {
            const int gpu_index  = omp_get_thread_num();
            Process_SNP_Data_Thread_Launch(gpu_index, snp_data[gpu_index]);
        }
    }
    
    for(size_t gpu_index = 0 ; gpu_index < gpu_id_list.size(); gpu_index++){
        cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
        cudaErrorCheck(cudaFree(GPU_eigenvectors_transposed[gpu_index]));
        delete snp_data[gpu_index];
    }
    delete [] GPU_eigenvectors_transposed;
    delete [] snp_data;
}*/
/*
Multi_Trait_GPU_GWAS_Estimator::Multi_Trait_GPU_GWAS_Estimator(std::vector<int> _gpu_id_list,const float * const CPU_eigenvalues,\
                               const float * const CPU_mean_column, const float * const CPU_component_A,\
                               float * _h2r, float * _loglik, float2 * _beta, \
                               float2 * _beta_se,  float * _variance,\
                               size_t _precision, size_t _max_batch_size,  size_t \
                               _dataset_thread_count, size_t _blockSize,\
                               size_t _scale, size_t _total_snps, size_t _n_subjects, size_t _pitch){
    gpu_id_list = _gpu_id_list;
    h2r = _h2r;
    loglik = _loglik;
    beta = _beta;
    variance = _variance;
    beta_se = _beta_se;
    blockSize = _blockSize;
    scale = _scale;
    dataset_thread_count = _dataset_thread_count;
    total_snps = _total_snps;
    pitch = _pitch;
    n_subjects = _n_subjects;
    max_batch_size = _max_batch_size;
    eigenvalues = _eigenvalues;
    precision = _precision;
    gpu_data = new GPU_Data*[gpu_id_list.size()];
    context = new GPU_GWAS_Estimator_Context*[gpu_id_list.size()];
    results = new GPU_GWAS_Estimator_Data*[gpu_id_list.size()];
    GPU_Y = new float*[gpu_id_list.size()];
    GPU_mean_column = new float*[gpu_id_list.size()];
    GPU_component_A = new float*[gpu_id_list.size()];
    GPU_component_D = new float*[gpu_id_list.size()];
    GPU_eigenvalues = new float*[gpu_id_list.size()];
    allocate_gpu_memory();
    CPU_SNP_Data = new float[pitch*total_snps];
    
    for(size_t gpu_index = 0; gpu_index < gpu_id_list.size(); gpu_index++){
        
        cudaErrorCheck(cudaSetDevice(gpu_data_list[gpu_index]));
        
        cudaErrorCheck(cudaMemcpy(GPU_eigenvalues[gpu_index], CPU_eigenvalues, sizeof(float)*n_subjects, cudaMemcpyHostToDevice));
        
        cudaErrorCheck(cudaMemcpy(GPU_mean_column[gpu_index], CPU_mean_column, sizeof(float)*n_subjects, cudaMemcpyHostToDevice));
        
        cudaErrorCheck(cudaMemcpy(GPU_component_A[gpu_index], CPU_component_A, sizeof(float)*n_subjects, cudaMemcpyHostToDevice));
        
        context[gpu_index] = new GPU_GWAS_Estimator_Context(gpu_data_list[gpu_index], GPU_Y[gpu_index],GPU_mean_column[gpu_index], GPU_eigenvalues[gpu_index],\
                                                          GPU_component_A[gpu_index],GPU_component_D[gpu_index], \
                                                          n_subjects, max_batch_size,\
                                                          blockSize, scale, pitch);
        
        results[gpu_index] = new GPU_GWAS_Estimator_Data(gpu_data_list[gpu_index], h2r, loglik, beta, beta_se, variance, max_batch_size);
    }
    
    
    set_function_pointer_blockSize();
    
}*/


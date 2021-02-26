//
//  gpu-gwas-estimator-context.cc
//  
//
//  Created by Brian Donohue on 8/18/18.
//

#include "gpu-gwas-estimator-context.h"

GPU_GWAS_Estimator_Data::GPU_GWAS_Estimator_Data(GPU_Data *  const _gpu_data, \
                                                 const int _max_batch_size){
    max_batch_size = _max_batch_size;
    gpu_data = _gpu_data;
    cudaErrorCheck(cudaMallocHost((void**)&cpu_parameters, sizeof(GPU_GWAS_TYPE)*5*max_batch_size));
  //  cpu_parameters = _cpu_parameters;
    //cpu_loglik = _cpu_loglik;
   // is_memory_pinned = false;
  /*  cudaErrorCheck(cudaMallocHost((void**)&cpu_pinned_h2r, sizeof(float)*max_batch_size));
    cudaErrorCheck(cudaMallocHost((void**)&cpu_pinned_loglik, sizeof(float)*max_batch_size));
    cudaErrorCheck(cudaMallocHost((void**)&cpu_pinned_beta, sizeof(float2)*max_batch_size));
    cudaErrorCheck(cudaMallocHost((void**)&cpu_pinned_beta_se, sizeof(float2)*max_batch_size));
    cudaErrorCheck(cudaMallocHost((void**)&cpu_pinned_variance, sizeof(float)*max_batch_size));*/
    allocate_gpu_memory();
}


GPU_GWAS_Estimator_Data::~GPU_GWAS_Estimator_Data(){
  //  unpin_cpu_results();
   // cpu_pinned_parameters = 0;
   // cpu_pinned_loglik = 0;
    cudaErrorCheck(cudaFreeHost(cpu_parameters));
   // cpu_parameters = 0;
    //cpu_loglik =0;

  /*  cudaErrorCheck(cudaFreeHost(cpu_pinned_h2r));
    cudaErrorCheck(cudaFreeHost(cpu_pinned_loglik));
    cudaErrorCheck(cudaFreeHost(cpu_pinned_beta));
    cudaErrorCheck(cudaFreeHost(cpu_pinned_beta_se));
    cudaErrorCheck(cudaFreeHost(cpu_pinned_variance));*/
    free_gpu_memory();
}


void GPU_GWAS_Estimator_Data::copy_results_to_cpu(const int batch_size){
   // cudaErrorCheck(cudaMemcpyAsync(cpu_pinned_parameters, parameters, sizeof(GPU_GWAS_TYPE)*batch_size*4, cudaMemcpyDeviceToHost, gpu_data->stream));
    cudaErrorCheck(cudaMemcpyAsync(cpu_parameters, parameters, sizeof(GPU_GWAS_TYPE)*batch_size*5, cudaMemcpyDeviceToHost, gpu_data->stream));
    //cudaErrorCheck(cudaMemcpyAsync(cpu_pinned_loglik, loglik, sizeof(GPU_GWAS_TYPE)*batch_size, cudaMemcpyDeviceToHost, gpu_data->stream));
   // cudaErrorCheck(cudaMemcpyAsync(cpu_loglik + index, loglik, sizeof(GPU_GWAS_TYPE)*batch_size, cudaMemcpyDeviceToHost, gpu_data->stream));	
	
}
/*
void GPU_GWAS_Estimator_Data::copy_results_to_unpinned_cpu(const size_t start, const size_t batch_size){
    memcpy(cpu_h2r + start, cpu_pinned_h2r, sizeof(float)*batch_size);
    memcpy(cpu_loglik + start, cpu_pinned_loglik, sizeof(float)*batch_size);
    memcpy(cpu_beta + start, cpu_pinned_beta, sizeof(float2)*batch_size);
    memcpy(cpu_beta_se + start, cpu_pinned_beta_se, sizeof(float2)*batch_size);
    memcpy(cpu_variance + start, cpu_pinned_variance, sizeof(float)*batch_size);
    
}*/

void GPU_GWAS_Estimator_Data::free_gpu_memory(){
    cudaErrorCheck(cudaFree(parameters));
   // cudaErrorCheck(cudaFree(loglik));

    
}

void GPU_GWAS_Estimator_Data::allocate_gpu_memory(){
    
    cudaErrorCheck(cudaMalloc((void**)&parameters, sizeof(GPU_GWAS_TYPE)*max_batch_size*5));
   // cudaErrorCheck(cudaMalloc((void**)&loglik, sizeof(GPU_GWAS_TYPE)*max_batch_size));
    
}



GPU_GWAS_Estimator_SNP_Data::GPU_GWAS_Estimator_SNP_Data(const GPU_GWAS_TYPE  * const _eigenvectors, pio_file_t * const file,GPU_Data * const _gpu_data, \
                                                         const int _blockSize, const int _n_subjects,const int _total_snps, const int batch_size,const int _pitch, \
                                                         const int * const _index_map){
    eigenvectors_transposed = _eigenvectors;
    gpu_data = _gpu_data;
    plink_file = file;
    blockSize = _blockSize;
    total_snps = _total_snps;
   // scale = _scale;
    max_batch_size = batch_size;
    pitch = _pitch;
    n_subjects = _n_subjects;
    index_map = _index_map;
    n_snps = batch_size;
    snp_buffer = new snp_t[plink_file->bed_file.header.num_samples*n_snps];
// cudaErrorCheck(cudaMallocHost((void**)&CPU_SNP_data,sizeof(GPU_GWAS_TYPE)*pitch*batch_size,cudaHostAllocWriteCombined));
    CPU_SNP_data = new GPU_GWAS_TYPE[pitch*batch_size];
    memset(CPU_SNP_data, 0, sizeof(GPU_GWAS_TYPE)*pitch*max_batch_size);
    start = 0;
    next_start =0;
    next_n_snps =0;
    SNP_data = 0;
    temp_SNP_data = 0;
    allocate_gpu_memory();
}

GPU_GWAS_Estimator_SNP_Data::~GPU_GWAS_Estimator_SNP_Data(){
    if(snp_buffer) delete [] snp_buffer;
    snp_buffer = 0;
    eigenvectors_transposed = 0;
    index_map = 0;
    free_gpu_memory();
   // cudaErrorCheck(cudaFreeHost(CPU_SNP_data));
   delete [] CPU_SNP_data;
}

int GPU_GWAS_Estimator_SNP_Data::prepare_snp_data(){
    bool skip = true;
    int end;
    
    #pragma omp critical
    {
    	next_start = plink_file->bed_file.cur_row;
   	if(next_start < total_snps){
		skip = false;
		next_n_snps =  max_batch_size;
		end = next_start + max_batch_size;
		if(end > total_snps){
			next_n_snps = total_snps - next_start;
		}
		for(int snp =0; snp < next_n_snps; snp++){
			pio_next_row(plink_file, snp_buffer + snp*plink_file->bed_file.header.num_samples);	
		}
     	}else{

		next_n_snps = 0;
    	 }
     }
    if(!skip){
	for(int snp = 0; snp < next_n_snps; snp++){
		//double mean = 0.0;
		//int n = n_subjects;
		GPU_GWAS_TYPE * CPU_SNP_data_ptr = CPU_SNP_data + snp*pitch;
		snp_t * snp_buffer_ptr = snp_buffer + snp*plink_file->bed_file.header.num_samples;
		for(int row =0 ; row < plink_file->bed_file.header.num_samples;row++){
                	const int index = index_map[row];
               		 if(index != -1){
                    	 	 const snp_t value = snp_buffer_ptr[row];

                    		 CPU_SNP_data_ptr[index] = value;
			  }
                }
		/*mean /= n;
		for(int row = 0; row < n_subjects; row++){
			if(CPU_SNP_data_ptr[row] != 3){
				CPU_SNP_data_ptr[row] -= mean;
			}else{
				CPU_SNP_data_ptr[row] = 0;
			}

		}*/
        }
    }
    /*if(next_start < total_snps){
        
        skip = false;
        end = next_start + max_batch_size;
        next_n_snps = max_batch_size;
        if(end > total_snps){
            end = total_snps;
            next_n_snps = total_snps  - next_start;
            //  eigen_CPU_SNP_data = eigen_CPU_SNP_data.block(0, 0, n_subjects, n_snps);
        }
        
        for(int snp =0 ;snp < next_n_snps; snp++){
            pio_next_row(plink_file, snp_buffer);
            //   locus = pio_get_locus(plink_file, snp +next_start);
            // free(SNP_names[snp]);
            //    SNP_names[snp] = strdup((const char*) locus->name);
            GPU_GWAS_TYPE * CPU_SNP_data_ptr = CPU_SNP_data + snp*pitch;
            for(int row =0 ; row < plink_file->bed_file.header.num_samples;row++){
                const int index = index_map[row];
                if(index != -1){
                    const snp_t value = snp_buffer[row];
                    CPU_SNP_data_ptr[index] = value;
                }
            }
        }
    }
    //  }
    if(skip) return 0;
    
    */
    return next_n_snps;
    
}

void GPU_GWAS_Estimator_SNP_Data::copy_SNP_data_to_gpu(){
    n_snps = next_n_snps;
    gpu_data->n_data_sets = n_snps;
    start = next_start;
    cudaErrorCheck(cudaMemcpyAsync(temp_SNP_data, CPU_SNP_data, sizeof(GPU_GWAS_TYPE)*pitch*n_snps, cudaMemcpyHostToDevice, gpu_data->stream));
    
}



void GPU_GWAS_Estimator_SNP_Data::allocate_gpu_memory(){
    
    cudaErrorCheck(cudaMalloc((void**)&SNP_data,sizeof(GPU_GWAS_TYPE)*pitch*max_batch_size));
    cudaErrorCheck(cudaMalloc((void**)&temp_SNP_data,sizeof(GPU_GWAS_TYPE)*pitch*max_batch_size));
}

void GPU_GWAS_Estimator_SNP_Data::free_gpu_memory(){
    
    cudaErrorCheck(cudaFree(SNP_data));
    cudaErrorCheck(cudaFree(temp_SNP_data));
}

void GPU_GWAS_Estimator_Context::allocate_gpu_memory(){
    
    cudaErrorCheck(cudaMalloc((void**)&omega, sizeof(GPU_GWAS_TYPE)*pitch*max_batch_size));
    cudaErrorCheck(cudaMalloc((void**)&raw_beta_components_B_C_E, sizeof(GPU_GWAS_TYPE)*pitch*max_batch_size*3));
 
  //  cudaErrorCheck(cudaMalloc((void**)&components_A_D, sizeof(double)*12*2*max_batch_size));


}

void GPU_GWAS_Estimator_Context::free_gpu_memory(){
    cudaErrorCheck(cudaFree(omega));
    cudaErrorCheck(cudaFree(raw_beta_components_B_C_E));
   
}

GPU_GWAS_Estimator_Context::GPU_GWAS_Estimator_Context(GPU_Data * const _gpu_data, GPU_GWAS_TYPE * const _Y_component_D_mean_column_component_A,  \
                           GPU_GWAS_TYPE  * const _eigenvalues, const int _n_subjects, const int _max_batch_size,\
			   const int _blockSize, const int _pitch, const int _precision){
    blockSize = _blockSize;
    Y_component_D_mean_column_component_A = _Y_component_D_mean_column_component_A;
    eigenvalues = _eigenvalues;
    precision = _precision;
    pitch = _pitch;
    max_batch_size = _max_batch_size;
    n_subjects = _n_subjects;
    raw_beta_components_B_C_E = 0;
    
    omega =0;

    gpu_data = _gpu_data;
    allocate_gpu_memory();
}

GPU_GWAS_Estimator_Context::~GPU_GWAS_Estimator_Context(){
    free_gpu_memory();
}


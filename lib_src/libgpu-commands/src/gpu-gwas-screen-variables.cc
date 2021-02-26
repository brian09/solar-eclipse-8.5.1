#include "gpu-gwas-screen-variables.h"
#include "gpu-exception.h"
GPU_GWAS_Screen_Data::GPU_GWAS_Screen_Data(float * _Sigma_Y, float * _Sigma,   const int _n_subjects , const int _max_batch_size , const int _pitch){
	Sigma_Y = _Sigma_Y;
	Sigma = _Sigma;
	max_batch_size  = _max_batch_size;
	pitch  = _pitch;
	n_subjects = _n_subjects;
	cudaErrorCheck(cudaMallocHost((void**)&cpu_chi_squared, sizeof(float)*max_batch_size));
	cudaErrorCheck(cudaMallocHost((void**)&cpu_beta_data, sizeof(float)*max_batch_size*2));
	test_statistic_denom = beta_data = chi_squared = 0;
	allocate_gpu_memory();
	//test_statistic_denom = beta_data = chi_squared = 0;
}
void GPU_GWAS_Screen_Data::allocate_gpu_memory(){ 
	if(!test_statistic_denom){
		cudaErrorCheck(cudaMalloc((void**)&test_statistic_denom, sizeof(float)*max_batch_size));
	}
	if(!chi_squared){
		cudaErrorCheck(cudaMalloc((void**)&chi_squared, sizeof(float)*max_batch_size));
	}
	if(!beta_data){
		cudaErrorCheck(cudaMalloc((void**)&beta_data, sizeof(float)*max_batch_size*2));
	}		
}

void GPU_GWAS_Screen_Data::free_gpu_memory(){
	cudaErrorCheck(cudaFreeHost(cpu_chi_squared));
	cudaErrorCheck(cudaFreeHost(cpu_beta_data));
	if(test_statistic_denom){
		cudaErrorCheck(cudaFree(test_statistic_denom));
	}
	if(chi_squared){
		cudaErrorCheck(cudaFree(chi_squared));
	}
	if(beta_data){
		cudaErrorCheck(cudaFree(beta_data));
	}		
}

GPU_GWAS_Screen_Data::~GPU_GWAS_Screen_Data(){
	free_gpu_memory();
	Sigma_Y  = 0;
	Sigma = 0;
}

void GPU_GWAS_Screen_Data::copy_results_to_cpu(const int batch_size, cudaStream_t stream){
	cudaErrorCheck(cudaMemcpyAsync(cpu_chi_squared, chi_squared, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream));
	cudaErrorCheck(cudaMemcpyAsync(cpu_beta_data , beta_data, sizeof(float)*batch_size*2, cudaMemcpyDeviceToHost, stream));
} 

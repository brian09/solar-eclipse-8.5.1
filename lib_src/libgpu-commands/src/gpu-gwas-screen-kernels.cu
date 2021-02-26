#include "gpu-exception.h"
#include "gpu-gwas-screen-variables.h"
#include "gpu-gwas-estimator-context.h"
static __device__ void add_reduce_shared_memory(float * shared_mem){

	for (int stride=blockDim.x/2;stride>32;stride>>=1){
		if (threadIdx.x<stride) shared_mem[threadIdx.x]+=shared_mem[threadIdx.x+stride];
		__syncthreads();
	}
	if(threadIdx.x < 32){
		volatile float * volatile_mem = shared_mem;
		volatile_mem[threadIdx.x] += volatile_mem[threadIdx.x + 32];
		volatile_mem[threadIdx.x] += volatile_mem[threadIdx.x + 16];
		volatile_mem[threadIdx.x] += volatile_mem[threadIdx.x + 8];
		volatile_mem[threadIdx.x] += volatile_mem[threadIdx.x + 4];
		volatile_mem[threadIdx.x] += volatile_mem[threadIdx.x + 2];
		volatile_mem[threadIdx.x] += volatile_mem[threadIdx.x + 1];
	}

	
}
/*
 Shared GPU Variables
....eigenvectors
....Sigma_Y
....Sigma
*/
/*
Local GPU Variables
....snp_data
....temp_snp_data
....sum
....n_subjects_missing
....chi_squared
....test_statistic_denom
*/

static __global__ void demean_snp_data(float *  const snp_data,  const int pitch, const int n_subjects){
	extern __shared__ float shared_sum [];
	__shared__ unsigned int missing_subjects;
	__shared__ float shared_mean;
	if(!threadIdx.x) missing_subjects  = 0;
	__syncthreads();
	int row_index;
	float l_sum = 0.f;
	for(row_index = threadIdx.x; row_index < pitch ; row_index += blockDim.x){
		const float snp_value =  snp_data[blockIdx.x*pitch + row_index];
		if(snp_value == 3) atomicInc(&missing_subjects, n_subjects);
		if(snp_value != 3) l_sum += snp_value;
	}
	shared_sum[threadIdx.x] = l_sum;
	__syncthreads();	
	add_reduce_shared_memory(shared_sum);
	if(!threadIdx.x) shared_mean = shared_sum[0]/(n_subjects - missing_subjects);
	__syncthreads();
	
	for(row_index = threadIdx.x ; row_index < n_subjects; row_index += blockDim.x){
		const float snp_value  =  snp_data[blockIdx.x*pitch + row_index];
		if(snp_value != 3){
			 snp_data[blockIdx.x*pitch + row_index] -= shared_mean;
		}else{
			snp_data[blockIdx.x*pitch + row_index]  = 0;
		}
	}
}


static __global__ void square_snp_data(float * const snp_data, const int pitch){

	for(int row_index = threadIdx.x; row_index < pitch; row_index+= blockDim.x){
		const double snp_value = snp_data[blockIdx.x*pitch + row_index];
		
		snp_data[blockIdx.x*pitch + row_index] *= snp_value;
	} 

}

static __global__ void calculate_chi_beta_and_beta_se(float * chi_squared, float * beta_data,  const float * test_statistic_denom, const int max_batch_size){
	const int index = threadIdx.x + blockDim.x*blockIdx.x;
	if(index >= max_batch_size) return;
	const float numerator = chi_squared[index];
	const float denominator = test_statistic_denom[index];
	beta_data[index*2] = numerator/denominator;
	beta_data[index*2 + 1] = rsqrtf(denominator);
	float l_chi_squared = numerator*numerator/denominator;
	if(l_chi_squared != l_chi_squared || l_chi_squared < 0) l_chi_squared = 0;
	
	chi_squared[index] = l_chi_squared;
}
 
void gpu_gwas_screen_kernels(GPU_Data * gpu_data, GPU_GWAS_Screen_Data * gpu_screen_data, GPU_GWAS_Estimator_SNP_Data * snp_data){
	cudaStream_t stream = gpu_data->stream;
	//cublasHandle_t handle = gpu_data->handle;
	demean_snp_data<<<gpu_data->gridSize(), gpu_data->blockSize(), sizeof(float)*gpu_data->thread_count, stream>>>(snp_data->temp_SNP_data,gpu_screen_data->pitch, gpu_screen_data->n_subjects);
	const float alpha = 1.f;
	const float beta = 0.f;
	cudaErrorCheck(cudaMemsetAsync(gpu_screen_data->chi_squared, 0, sizeof(float)*gpu_data->n_data_sets, stream));
	cudaErrorCheck(cudaMemsetAsync(gpu_screen_data->test_statistic_denom, 0, sizeof(float)*gpu_data->n_data_sets, stream))
	cublasErrorCheck(cublasSgemm(gpu_data->cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,gpu_screen_data->pitch,\
			 gpu_data->n_data_sets, gpu_screen_data->pitch, &alpha, \
			 snp_data->eigenvectors_transposed, gpu_screen_data->pitch,\
			 snp_data->temp_SNP_data, gpu_screen_data->pitch,&beta, \ 
			 snp_data->SNP_data, gpu_screen_data->pitch));


	cublasErrorCheck(cublasSgemv(gpu_data->cublas_handle, CUBLAS_OP_T,gpu_screen_data->pitch, gpu_data->n_data_sets, &alpha, \
			snp_data->SNP_data, gpu_screen_data->pitch, gpu_screen_data->Sigma_Y, 1, \
			&beta, gpu_screen_data->chi_squared, 1));	

	
 	square_snp_data<<<gpu_data->gridSize(), gpu_data->blockSize(), 0, stream>>>(snp_data->SNP_data, gpu_screen_data->pitch);

	cublasErrorCheck(cublasSgemv(gpu_data->cublas_handle, CUBLAS_OP_T,gpu_screen_data->pitch, gpu_data->n_data_sets, &alpha, \
			snp_data->SNP_data, gpu_screen_data->pitch, gpu_screen_data->Sigma, 1, \
			&beta, gpu_screen_data->test_statistic_denom, 1));	
	dim3 results_gridSize(gpu_data->n_data_sets/64.f, 1, 1);
	dim3 results_blockSize(64, 1, 1);
	calculate_chi_beta_and_beta_se<<<results_gridSize,results_blockSize, 0, stream>>>(gpu_screen_data->chi_squared, gpu_screen_data->beta_data,  gpu_screen_data->test_statistic_denom, gpu_data->n_data_sets);	
}


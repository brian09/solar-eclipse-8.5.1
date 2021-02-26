#include "gpu-exception.h"
#include "gpu-fphi-variables.h"
static const unsigned FULL_MASK = 0xffffffff;

__inline__ __device__ 
static float reduce_warp(float value, const unsigned mask = FULL_MASK){
	
	value+=__shfl_down_sync(mask, value, 16,32);
	
	value+=__shfl_down_sync(mask, value, 8,32);
	
	value+=__shfl_down_sync(mask, value, 4,32);
	
	value+=__shfl_down_sync(mask, value, 2,32);
       
	value+=__shfl_down_sync(mask, value, 1,32);
	return value;
}
template<int blockSize>
__inline__ __device__ 
static float warp_reduction(float * shared_data, float input){
	unsigned mask = FULL_MASK; 
	input = reduce_warp(input, mask);
	if(blockSize != 32){
		const int lane = threadIdx.x % 32;
		//const int wid = threadIdx.x/32;
		//static __shared__ GPU_GWAS_TYPE shared_data[blockSize/32];
		if(lane == 0) shared_data[threadIdx.x/32] = input;
		__syncthreads();
		input = (threadIdx.x < blockSize / 32) ? shared_data[lane] : 0;
		if(threadIdx.x/32 == 0) input = reduce_warp(input);
	}
	return input;
}

template<int blockSize>
static __global__ void demean_columns(float * trait_matrix, const int pitch, const int n_subjects){

	__shared__ float shared_data[blockSize/32];
	__shared__ float mean;
	int row_index;
	float sum = 0.f;
	for(row_index = threadIdx.x ; row_index < pitch ;row_index += blockSize){
		sum += trait_matrix[blockIdx.x*pitch + row_index];
	}
	
	sum = warp_reduction<blockSize>(shared_data, sum);
	if(!threadIdx.x) mean = sum/n_subjects;
	__syncthreads();
	
	for(row_index = threadIdx.x ; row_index < n_subjects ;row_index += blockSize){
		trait_matrix[blockIdx.x*pitch + row_index] -= mean;
	}	
}

template<int blockSize>
static __global__ void calculate_theta_and_sigma(float2 * __restrict__ theta, float * __restrict__ sigma, const float * __restrict__ trait_matrix, \
		 const float * __restrict__ eigenvalues,  const float ZTZI_0, const float ZTZI_1, const float ZTZI_2, const float ZTZI_3, const int pitch, const int n_subjects){

   int row_index;
   __shared__ float shared_data[blockSize/32];
   float sum_0 = 0.f;
   float sum_1 = 0.f;
   float sum_sigma = 0.f;
   float2 result;
   for(row_index = threadIdx.x; row_index < n_subjects; row_index += blockSize){
   	float F = trait_matrix[blockIdx.x*pitch + row_index];
   	const float eigenvalue = eigenvalues[row_index];
   	F *= F;
   	sum_0 += F*(ZTZI_0 + ZTZI_2*eigenvalue);
   	sum_1 += F*(ZTZI_1 + ZTZI_3*eigenvalue);
   	sum_sigma += F;
   }
   sum_sigma = warp_reduction<blockSize>(shared_data,sum_sigma);
   if(!threadIdx.x) sigma[blockIdx.x] = sum_sigma;
   __syncthreads();
   sum_0 = warp_reduction<blockSize>(shared_data, sum_0);
   if(!threadIdx.x) result.x = sum_0;
   __syncthreads();
   sum_1 = warp_reduction<blockSize>(shared_data, sum_1);
   if(!threadIdx.x){
   	result.y = sum_1;
   	theta[blockIdx.x] = result;
   }
   
 }
 
/*
static __global__ void calculate_weights(float * __restrict__ weights, const float * __restrict__ eigenvalues,
		const float2 * __restrict__ theta, const int pitch, const int n_subjects){
   __shared__ float2 shared_theta;
   if(!threadIdx.x) shared_theta = theta[blockIdx.x];
   int row_index ;
   for(row_index = threadIdx.x; row_index < n_subjects; row_index += blockDim.x){
   	  const float eigenvalue = eigenvalues[row_index];
   	  float weight  = shared_theta.x + shared_theta.y*eigenvalue;
   	  if(weight != 0.f) 
   	  	weights[blockIdx.x*pitch + row_index] = 1.f/(weight*weight);
   	  else
   	  	weights[blockIdx.x*pitch + row_index] = 0.f;
   }
   
   if(row_index >= n_subjects && row_index < pitch) weights[blockIdx.x*pitch + row_index] = 0.f;
   	  
}*/

 

template <int blockSize>
static __global__ void calculate_sigma_a_and_sigma_e(float2 * __restrict__ sigma_e_and_a, float * __restrict__ h2,\
						const float2 * __restrict__ theta, const float * __restrict__ trait_matrix, const float * __restrict__ eigenvalues,\
						const int pitch, const int n_subjects){
						
  __shared__ float shared_data[blockSize/32];
  __shared__ float2 shared_theta;
  if(!threadIdx.x) shared_theta = theta[blockIdx.x];
  __syncthreads();
 
  int row_index;
  float sum_A = 0.f;
  float sum_B = 0.f;
  float sum_C = 0.f;
  float sum_D = 0.f;
  float sum_E = 0.f;

  for(row_index = threadIdx.x ; row_index < n_subjects; row_index+=blockSize){
 	const float lambda = eigenvalues[row_index];
        float weight = lambda*shared_theta.y + shared_theta.x;
        if(weight != 0.f)
        	weight = 1.f/(weight*weight);
        else
        	weight = 0.f;
  	float F = trait_matrix[blockIdx.x*pitch + row_index];
  	F *= F;
  	
   	sum_A += weight;
  	sum_D += weight*F;
  	//sum_sigma += F; 	
  	sum_B += weight*lambda;
  	sum_C += weight*lambda*lambda;
  	sum_E += F*weight*lambda;
  	
  }
    sum_A = warp_reduction<blockSize>(shared_data, sum_A);
  //if(!threadIdx.x) A[blockIdx.x] = sum_A;
  __syncthreads();
  
   sum_D = warp_reduction<blockSize>(shared_data, sum_D);
 // if(!threadIdx.x) D[blockIdx.x] = sum_D;
  __syncthreads();
  
  // sum_sigma = warp_reduction<blockSize>(shared_data, sum_sigma);
 // if(!threadIdx.x) Sigma[blockIdx.x] = sum_sigma;
 // __syncthreads();
  
  sum_B = warp_reduction<blockSize>(shared_data, sum_B);
 // if(!threadIdx.x) B[blockIdx.x] = sum_B;
  __syncthreads();
  
   sum_C = warp_reduction<blockSize>(shared_data, sum_C);
 // if(!threadIdx.x) C[blockIdx.x] = sum_C;
  __syncthreads();
  
   sum_E = warp_reduction<blockSize>(shared_data, sum_E);
   if(!threadIdx.x){
   	const float denom = sum_A*sum_C - sum_B*sum_B;
   	float2 l_sigma_e_and_a;
   	l_sigma_e_and_a.x = (sum_C*sum_D - sum_B*sum_E)/denom;
   	l_sigma_e_and_a.y = (sum_A*sum_E - sum_B*sum_D)/denom;
   	sigma_e_and_a[blockIdx.x] = l_sigma_e_and_a;
   }
   
 // if(!threadIdx.x) E[blockIdx.x] = sum_E;

}


template <size_t blockSize>
static __global__ void calculate_score(float * __restrict__ score, const float * __restrict__ eigenvalues, const float * __restrict__ trait_matrix,
										const float * __restrict__ Sigma,const int pitch, const int n_subjects){
  __shared__ float shared_data[blockSize/32];
  __shared__ float mean;
  if(!threadIdx.x) mean = Sigma[blockIdx.x]/n_subjects;
  __syncthreads();
  int row_index;
  
  float sum = 0.f;

  for(row_index = threadIdx.x ; row_index < n_subjects; row_index+=blockSize){
  	float F = trait_matrix[blockIdx.x*pitch + row_index];
  	F *= F;
  	const float lambda = eigenvalues[row_index];
  	sum += lambda*((F/mean) - 1.f);

  	
  }
  sum = warp_reduction<blockSize>(shared_data, sum);
  if(!threadIdx.x) score[blockIdx.x] = sum;

}
template<int blockSize>
static __global__ void calculate_h2_and_SE(bool * __restrict__ boolean_score, float * __restrict__ h2,  float *  __restrict__ indicator,
		float * __restrict__ SE, const float * __restrict__ eigenvalues, const float2 * __restrict__ sigma_e_and_a,\
		 const int pitch, const int n_subjects){
	__shared__ float shared_data[blockSize/32];
	__shared__ float2 shared_sigma_e_and_a;
	if(!threadIdx.x){
		shared_sigma_e_and_a = sigma_e_and_a[blockIdx.x];
	}
	__syncthreads();
	float sum_A = 0.f;
	float sum_B = 0.f;
	float sum_C = 0.f;
	int row_index;
	for(row_index = threadIdx.x; row_index < n_subjects; row_index += blockSize){
		const float eigenvalue = eigenvalues[row_index];
		float weight = shared_sigma_e_and_a.x + shared_sigma_e_and_a.y*eigenvalue;
		weight = 1.f/(weight*weight);
		sum_A += weight;
		sum_B += weight*eigenvalue;
		sum_C  += weight*eigenvalue*eigenvalue;
	}
	
	 sum_A = warp_reduction<blockSize>(shared_data, sum_A);
	 __syncthreads();
	 
	 	
	sum_B = warp_reduction<blockSize>(shared_data, sum_B);
	 __syncthreads();
	sum_C = warp_reduction<blockSize>(shared_data, sum_C);
	
	
	if(!threadIdx.x){
		const float score = indicator[blockIdx.x];
		const float denom = shared_sigma_e_and_a.x + shared_sigma_e_and_a.y;
		if(score < 0.f || denom == 0.f || shared_sigma_e_and_a.x != shared_sigma_e_and_a.x || shared_sigma_e_and_a.y != shared_sigma_e_and_a.y){
			h2[blockIdx.x] = nan("");
			boolean_score[blockIdx.x] = false;
			return;
		}
		
		h2[blockIdx.x] = shared_sigma_e_and_a.y/denom;
		const float G = shared_sigma_e_and_a.y/(denom*denom);
		const float E = shared_sigma_e_and_a.x/(denom*denom);
		const float det = sum_A*sum_C - sum_B*sum_B;
		const float var = 2.f*(G*G*sum_C + 2*G*E*sum_B + E*E*sum_A)/det;
		SE[blockIdx.x] = sqrt(var);
		boolean_score[blockIdx.x] = true;
	}
}
template<int blockSize>
static __global__ void calculate_chi_squared(const bool * __restrict__ boolean_score, float * __restrict__ chi_squared, const float2 * __restrict__ sigma_e_and_a, \
						const float * __restrict__ sigma,const float * __restrict__ eigenvalues, const int pitch, const int n_subjects){
 if(boolean_score[blockIdx.x] == false){
   if(!threadIdx.x) chi_squared[blockIdx.x] = -1.f;
   return;
  }
  __shared__ float shared_data[blockSize/32];
 							
   __shared__ float2 shared_sigma_e_and_a;   
   if(!threadIdx.x) shared_sigma_e_and_a = sigma_e_and_a[blockIdx.x];
   __syncthreads();
  
   float sum = 0.f;
   int row_index;
   for(row_index = threadIdx.x; row_index < n_subjects; row_index += blockSize){
   	sum += logf(shared_sigma_e_and_a.y*eigenvalues[row_index] + shared_sigma_e_and_a.x);
   }

   sum = warp_reduction<blockSize>(shared_data, sum); 
   
   if(!threadIdx.x) chi_squared[blockIdx.x] = (n_subjects*logf(sigma[blockIdx.x]/n_subjects) + n_subjects) - (sum + n_subjects);	
}
template<int blockSize>
void run_fphi_gpu_functions(GPU_Data * const gpu_data, GPU_FPHI_Shared_Variables * const shared_variables, GPU_FPHI_Stream_Variables * const stream_variables, GPU_FPHI_Results * const results){
	demean_columns<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(stream_variables->trait_matrix, shared_variables->pitch, shared_variables->n_subjects);
	const float alpha = 1.f;
	const float beta = 0.f;
	cublasErrorCheck(cublasSgemm(gpu_data->cublas_handle, CUBLAS_OP_N,CUBLAS_OP_N,\
		shared_variables->pitch, gpu_data->n_data_sets,shared_variables->pitch, &alpha,\
		shared_variables->eigenvectors_transposed,shared_variables->pitch,\
				stream_variables->trait_matrix, shared_variables->pitch,\
				&beta, stream_variables->temp_trait_matrix, shared_variables->pitch));		
	if(shared_variables->use_covariates){
		cublasErrorCheck(cublasSgemm(gpu_data->cublas_handle, CUBLAS_OP_N,CUBLAS_OP_N,\
			shared_variables->pitch, gpu_data->n_data_sets,shared_variables->pitch, &alpha,\
			shared_variables->hat_matrix,shared_variables->pitch,\
					stream_variables->temp_trait_matrix, shared_variables->pitch,\
					&beta, stream_variables->trait_matrix, shared_variables->pitch));	
	}else{
		cudaErrorCheck(cudaMemcpyAsync(stream_variables->trait_matrix, stream_variables->temp_trait_matrix, sizeof(float)*gpu_data->n_data_sets*shared_variables->pitch, cudaMemcpyDeviceToDevice, gpu_data->stream));
	}
	
	
	calculate_theta_and_sigma<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(stream_variables->theta, stream_variables->sigma, stream_variables->trait_matrix, \
		 shared_variables->eigenvalues,  shared_variables->ZTZI_0, shared_variables->ZTZI_1,shared_variables->ZTZI_2, \
		 shared_variables->ZTZI_3, shared_variables->pitch, shared_variables->n_subjects);


	calculate_sigma_a_and_sigma_e<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(stream_variables->sigma_e_and_a, results->h2r,\
						stream_variables->theta, stream_variables->trait_matrix, shared_variables->eigenvalues,\
						shared_variables->pitch, shared_variables->n_subjects);
						
						
	calculate_score<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(results->score, shared_variables->eigenvalues, stream_variables->trait_matrix,\
										stream_variables->sigma,shared_variables->pitch, shared_variables->n_subjects);	
										
	calculate_h2_and_SE<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(stream_variables->boolean_score, results->h2r,  results->score,
		results->SE, shared_variables->eigenvalues, stream_variables->sigma_e_and_a,\
		shared_variables->pitch,shared_variables->n_subjects);	

	calculate_chi_squared<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(stream_variables->boolean_score, results->chi_squared, stream_variables->sigma_e_and_a, \
						stream_variables->sigma, shared_variables->eigenvalues, shared_variables->pitch, shared_variables->n_subjects);				
						
}
void GPU_FPHI_Estimator::set_function_pointer(){
    switch(blockSize)
    {
        case 32:
            Run_GPU_FPHI_Function_Pointer = &run_fphi_gpu_functions<32>;
           
            break;
        case 64:
            Run_GPU_FPHI_Function_Pointer = &run_fphi_gpu_functions<64>;
          
            break;
        case 128:
            Run_GPU_FPHI_Function_Pointer = &run_fphi_gpu_functions<128>;
            
            break;
        case 256:
            Run_GPU_FPHI_Function_Pointer = &run_fphi_gpu_functions<256>;
           
            break;
            
        case 512:
            Run_GPU_FPHI_Function_Pointer = &run_fphi_gpu_functions<512>;
            
            break;

	case 1024:
            Run_GPU_FPHI_Function_Pointer = &run_fphi_gpu_functions<1024>;
        
            break;
        
       
    }
    	
	

}


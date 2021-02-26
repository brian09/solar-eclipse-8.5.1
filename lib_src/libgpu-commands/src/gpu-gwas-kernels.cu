#include "gpu-gwas-estimator.h"
#include "gpu-gwas-settings.h"

static const unsigned FULL_MASK = 0xffffffff;

__inline__ __device__ 
static GPU_GWAS_TYPE reduce_warp(GPU_GWAS_TYPE value, const unsigned mask = FULL_MASK){
	
	value+=__shfl_down_sync(mask, value, 16,32);
	
	value+=__shfl_down_sync(mask, value, 8,32);
	
	value+=__shfl_down_sync(mask, value, 4,32);
	
	value+=__shfl_down_sync(mask, value, 2,32);
       
	value+=__shfl_down_sync(mask, value, 1,32);
	return value;
}
template<int blockSize>
__inline__ __device__ 
static GPU_GWAS_TYPE warp_reduction(GPU_GWAS_TYPE * shared_data, GPU_GWAS_TYPE input){
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

template<int  blockSize>
__device__ static void   reduction (GPU_GWAS_TYPE  shared_data[blockSize]){
		if(blockSize >= 1024){
			if(threadIdx.x < 512){
				shared_data[threadIdx.x] += shared_data[threadIdx.x + 512];
			}
			__syncthreads();
		}
    
		if(blockSize >= 512){
			if(threadIdx.x < 256){
				shared_data[threadIdx.x] += shared_data[threadIdx.x + 256];
			}
			__syncthreads();
		}
		if(blockSize >= 256){
			if(threadIdx.x < 128){
				shared_data[threadIdx.x] += shared_data[threadIdx.x + 128];
			}
			__syncthreads();
		}
		if(blockSize >= 128){
			if(threadIdx.x < 64){
				shared_data[threadIdx.x] += shared_data[threadIdx.x + 64];
			}
			__syncthreads();
		}
		if(blockSize >= 64){
			if(threadIdx.x < 32){
				volatile GPU_GWAS_TYPE  * volatile_data  = shared_data;
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 32];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 16];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 8];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 4];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 2];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 1];
			}
		}else{
			if(threadIdx.x < 16){
				volatile GPU_GWAS_TYPE  * volatile_data  = shared_data;
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 16];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 8];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 4];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 2];
				volatile_data[threadIdx.x] += volatile_data[threadIdx.x + 1];			
			}
		
		
		}	
	
    
}

__device__ static __inline__ GPU_GWAS_TYPE calculate_constraint(const GPU_GWAS_TYPE x){
	//const GPU_GWAS_TYPE x_squared = x*x;
	return 1.f/(1.f + expf(-x));
	//return exp(x)/(1.0+exp(x));
}

__device__ static __inline__  GPU_GWAS_TYPE calculate_dconstraint(const GPU_GWAS_TYPE x){
	const GPU_GWAS_TYPE e_x = expf(-x);
	const GPU_GWAS_TYPE denom = e_x + 1.f;
	return e_x/(denom*denom);// //pow((1.0 + x*x), -2);
	//return exp(x)*pow(1 + exp(x), -2);
}

__device__ static __inline__  GPU_GWAS_TYPE calculate_ddconstraint(const GPU_GWAS_TYPE x){
	const GPU_GWAS_TYPE e_x = expf(x);
	const GPU_GWAS_TYPE denom = e_x + 1.f;
	return -(e_x-1.f)*e_x/(denom*denom*denom);//pow((x_squared + 1.0), -3);
	//return -exp(x)*(exp(x) - 1.0)*pow(exp(x) + 1.0, -3);
}
/*
gradient = dconstraint*0.5*(variance*(lambda - 1)*(sigma*(residual_squared - 1)))

*/
__device__ static __inline__  GPU_GWAS_TYPE calculate_dloglik(const GPU_GWAS_TYPE var, const GPU_GWAS_TYPE lambda, const GPU_GWAS_TYPE residual_squared, const GPU_GWAS_TYPE omega){
	return 0.5f*(lambda - 1.f)*(omega*(-1.f + residual_squared/var));
}
__device__ static __inline__ GPU_GWAS_TYPE calculate_dloglik_part_one(const GPU_GWAS_TYPE lambda, const GPU_GWAS_TYPE omega, const GPU_GWAS_TYPE residual_squared){
	return 0.5f*(lambda - 1.f)*omega*residual_squared;
}
__device__ static __inline__ GPU_GWAS_TYPE calculate_dloglik_part_two(const GPU_GWAS_TYPE lambda, const GPU_GWAS_TYPE omega){
	return 0.5f*(lambda - 1.f)*-omega;
}
__device__ static __inline__ GPU_GWAS_TYPE calculate_ddloglik(const GPU_GWAS_TYPE var, const GPU_GWAS_TYPE lambda, const GPU_GWAS_TYPE residual_squared, const GPU_GWAS_TYPE omega){
	const GPU_GWAS_TYPE lambda_minus_one_omega = (lambda - 1.f)*omega;
	return 0.5f*lambda_minus_one_omega*lambda_minus_one_omega*(1.f - 2.f*residual_squared/var);
}
__device__ static __inline__ GPU_GWAS_TYPE calculate_ddloglik_part_one(const GPU_GWAS_TYPE lambda_minus_one_omega, const GPU_GWAS_TYPE residual_squared){
	return 0.5f*lambda_minus_one_omega*lambda_minus_one_omega*-2.f*residual_squared;
}
__device__ static __inline__ GPU_GWAS_TYPE calculate_ddloglik_part_two(const GPU_GWAS_TYPE lambda_minus_one_omega){
	return 0.5f*lambda_minus_one_omega*lambda_minus_one_omega;
}


__device__ static __inline__ GPU_GWAS_TYPE calculate_ddloglik_with_constraint(const GPU_GWAS_TYPE t,const GPU_GWAS_TYPE dloglik, const GPU_GWAS_TYPE ddloglik){
	const GPU_GWAS_TYPE dconstraint = calculate_dconstraint(t);
	return dconstraint*dconstraint*ddloglik + calculate_ddconstraint(t)*dloglik;
}


__device__ static __inline__ GPU_GWAS_TYPE calculate_dloglik_with_constraint(const GPU_GWAS_TYPE t, const GPU_GWAS_TYPE dloglik){

	return calculate_dconstraint(t)*dloglik;
}

template< int blockSize>
__global__ static void initialize_parameters(GPU_GWAS_TYPE * const raw_beta_components_B_C_E, const GPU_GWAS_TYPE * const Y_component_D_mean_column_component_A,\
					     const GPU_GWAS_TYPE * const  snp_data, const int n_subjects, const int pitch){
	
	int row_index;// = threadIdx.x + blockIdx.x*blockSize;
	const int end = (int)floorf((float)n_subjects/blockSize)*blockSize;
#pragma unroll
	for(row_index = threadIdx.x; row_index < end; row_index += blockSize){
		const GPU_GWAS_TYPE snp_value = snp_data[blockIdx.x*pitch + row_index];
		raw_beta_components_B_C_E[blockIdx.x*3*pitch + row_index] = Y_component_D_mean_column_component_A[pitch*2 + row_index]*snp_value;
		raw_beta_components_B_C_E[blockIdx.x*3*pitch + pitch + row_index] = snp_value*snp_value;
		raw_beta_components_B_C_E[blockIdx.x*3*pitch + pitch*2+ row_index] = Y_component_D_mean_column_component_A[row_index]*snp_value;
		;
	}
	if(row_index < pitch){
		GPU_GWAS_TYPE snp_value = 0.f;
		if(row_index < n_subjects) snp_value = snp_data[blockIdx.x*pitch + row_index];
		raw_beta_components_B_C_E[blockIdx.x*3*pitch + row_index] = Y_component_D_mean_column_component_A[pitch*2+row_index]*snp_value;
		raw_beta_components_B_C_E[blockIdx.x*3*pitch + pitch + row_index] = snp_value*snp_value;
		raw_beta_components_B_C_E[blockIdx.x*3*pitch + pitch*2+ row_index] = Y_component_D_mean_column_component_A[row_index]*snp_value;	
	//raw_beta_components[blockIdx.x*2*pitch + row_index] = Y[row_index]*snp_value;
	//raw_beta_components[blockIdx.x*2*pitch + row_index + pitch] = snp_value*snp_value;

	}

}

template< int blockSize>
__global__ static void calculate_loglik(const GPU_GWAS_TYPE * const omega,  GPU_GWAS_TYPE * const parameters,  const int n_subjects, const int pitch){
	//__shared__ GPU_GWAS_TYPE shared_h2r;
	__shared__ GPU_GWAS_TYPE shared_data[blockSize/32];
	//if(threadIdx.x == 0) shared_h2r = parameters[blockIdx.x*4];
	//__syncthreads();
	//const GPU_GWAS_TYPE local_h2r = parameters[blockIdx.x*4];
	//if(local_h2r != local_h2r) return;
	int row_index;
	GPU_GWAS_TYPE sum = 0;
	const int end = (int)floorf((float)n_subjects/blockSize)*blockSize;
#pragma unroll
	for(row_index = threadIdx.x; row_index < end; row_index += blockSize){
		sum -= logf(omega[blockIdx.x*pitch + row_index]);
	}
	//if(row_index < n_subjects) sum -= logf(omega[blockIdx.x*pitch + row_index]);
	
	 sum = (row_index < n_subjects) ? sum - logf(omega[blockIdx.x*pitch + row_index]) : sum;
	//__shared__ GPU_GWAS_TYPE shared_data[blockSize];
	//shared_data[threadIdx.x] = sum;
	//__syncthreads();
	//if(blockSize == pitch)
	//	sum = warp_reduction<blockSize>(shared_data, sum, n_subjects);
	//else
		sum = warp_reduction<blockSize>(shared_data,sum);
	//reduction< blockSize>(shared_data);
	if(threadIdx.x == 0) parameters[blockIdx.x*5 + 4] = -0.5f*(sum + n_subjects*(1.f + logf(parameters[blockIdx.x*5 + 3])));
}
template<int blockSize>
__global__ static void maximize_h2r(GPU_GWAS_TYPE * __restrict__ const parameters, GPU_GWAS_TYPE * __restrict__ const omega,\
					  GPU_GWAS_TYPE * __restrict__ const raw_beta_components_B_C_E,  \
					  const GPU_GWAS_TYPE * __restrict__ const Y_component_D_mean_column_component_A, const GPU_GWAS_TYPE * __restrict__ const eigenvalues, const GPU_GWAS_TYPE * __restrict__ const snp_data, const int n_subjects,\
					 const int pitch, const float precision){
	//GPU_GWAS_TYPE * l_omega = omega + blockIdx.x*pitch;
	//GPU_GWAS_TYPE * l_raw_beta_components_B_C_E = raw_beta_components_B_C_E + blockIdx.x*pitch*3;
	//const GPU_GWAS_TYPE * l_snp_data = snp_data + blockIdx.x*pitch;
	//__shared__ GPU_GWAS_TYPE parameter_t;
 	const int end = (int)floorf((float)n_subjects/blockSize)*blockSize;
 	GPU_GWAS_TYPE value[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
 	int row_index;
	__shared__ GPU_GWAS_TYPE next_h2r;
	GPU_GWAS_TYPE parameter_t = 0.f;
	__shared__ bool continue_loop;
	__shared__ GPU_GWAS_TYPE shared_parameters[5];
	#ifdef USE_WARPS
	__shared__ GPU_GWAS_TYPE shared_data[blockSize/32];
	#endif 
	#ifdef USE_SHARED_MEMORY
	__shared__ GPU_GWAS_TYPE shared_data[blockSize];
	#endif
	int iteration_count = 0;
	if(!threadIdx.x) next_h2r = 0.5f;
	
	if(!threadIdx.x) continue_loop = true;
	__syncthreads();					 
					 
	do{
		if(!threadIdx.x) shared_parameters[0] = next_h2r;
	#pragma unroll
		for(row_index = threadIdx.x; row_index < end; row_index += blockSize){
			omega[blockIdx.x*pitch + row_index] = 1.f/((eigenvalues[row_index]-1.f)*next_h2r + 1.f);
		}
		if(row_index < pitch){
			GPU_GWAS_TYPE local_omega = 0.f;
			if(row_index < n_subjects) local_omega = 1.f/((eigenvalues[row_index]-1.f)*next_h2r + 1.f);
			omega[blockIdx.x*pitch + row_index] = local_omega;
		}	   
	   
		//__syncthreads();

		value[0] = value[1] = value[2] = value[3] = value[4] = 0.f;
	//GPU_GWAS_TYPE A = 0.f, B = 0.f, C = 0.f, D = 0.f, E = 0.f;
	//GPU_GWAS_TYPE sum = 0.f;
	#pragma unroll
		for(row_index = threadIdx.x; row_index < pitch; row_index += blockSize){
			const GPU_GWAS_TYPE local_omega = omega[blockIdx.x*pitch + row_index];
			value[3] += Y_component_D_mean_column_component_A[pitch + row_index]*local_omega;
			value[0] += Y_component_D_mean_column_component_A[pitch*3 + row_index]*local_omega;
			value[1] += raw_beta_components_B_C_E[blockIdx.x*3*pitch+row_index]*local_omega;
			value[2] += raw_beta_components_B_C_E[blockIdx.x*3*pitch+pitch + row_index]*local_omega;
			value[4] += raw_beta_components_B_C_E[blockIdx.x*3*pitch+ pitch*2 + row_index]*local_omega;
			
		}
	#ifdef USE_WARPS
		value[0] = warp_reduction<blockSize>(shared_data,value[0]);
		__syncthreads();
		value[1] = warp_reduction<blockSize>(shared_data,value[1]);
		__syncthreads();
		value[2] = warp_reduction<blockSize>(shared_data,value[2]);
		__syncthreads();
		value[3] = warp_reduction<blockSize>(shared_data,value[3]);
		__syncthreads();
		value[4] = warp_reduction<blockSize>(shared_data,value[4]);
	#endif
	
	#ifdef USE_SHARED_MEMORY
		shared_data[threadIdx.x] = value[0];
		__syncthreads();
		reduction< blockSize>(shared_data);
		if(!threadIdx.x) value[0] = shared_data[0];
		shared_data[threadIdx.x] = value[1];
		__syncthreads();
		reduction<blockSize>(shared_data);
		if(!threadIdx.x) value[1] = shared_data[0];
		shared_data[threadIdx.x] = value[2];
		__syncthreads();
		reduction<blockSize>(shared_data);
		if(!threadIdx.x) value[2] = shared_data[0];
		shared_data[threadIdx.x] = value[3];
		__syncthreads();
		reduction<blockSize>(shared_data);
		if(!threadIdx.x) value[2] = shared_data[0];
		shared_data[threadIdx.x] = value[3];
		__syncthreads();
		reduction<blockSize>(shared_data);
		if(!threadIdx.x) value[3] = shared_data[0];
		shared_data[threadIdx.x] = value[4];
		__syncthreads();
		reduction<blockSize>(shared_data);
		if(!threadIdx.x) value[4] = shared_data[0];
	#endif 
									
		if(!threadIdx.x){
			const GPU_GWAS_TYPE denom = value[0]*value[2] - value[1]*value[1];
			shared_parameters[4] = (value[2]*value[3] - value[1]*value[4])/denom;
		//parameters[blockIdx.x*4 + 1] 
			shared_parameters[1] = (value[0]*value[4] - value[1]*value[3])/denom;	
		//parameters[blockIdx.x*4 + 2] 
			shared_parameters[2] = sqrtf(value[0]/denom);
		}		


		__syncthreads();

		value[0] = value[1] = value[2] = 0.f;

		#pragma unroll
		for(row_index = threadIdx.x; row_index < pitch; row_index += blockSize){
			GPU_GWAS_TYPE local_residual_squared =  Y_component_D_mean_column_component_A[row_index];
			local_residual_squared -= shared_parameters[1]*snp_data[blockIdx.x*pitch + row_index];
			local_residual_squared -= Y_component_D_mean_column_component_A[pitch*2 + row_index]*shared_parameters[4];

			const GPU_GWAS_TYPE local_omega =omega[blockIdx.x*pitch + row_index];
			const GPU_GWAS_TYPE lambda_minus_one_omega = (eigenvalues[row_index] - 1.f)*local_omega;
			local_residual_squared *= local_residual_squared*local_omega;
		//residual_squared[blockIdx.x*pitch + row_index] = local_residual_squared;
			value[0] += local_residual_squared;
			value[1] +=  lambda_minus_one_omega*local_residual_squared;
			value[2] += lambda_minus_one_omega*lambda_minus_one_omega*-2.f*local_residual_squared;
		}
		#ifdef USE_WARPS
			value[0] = warp_reduction<blockSize>(shared_data, value[0]);
		#endif
	
		#ifdef USE_SHARED_MEMORY
			shared_data[threadIdx.x] = value[0];
			__syncthreads();
			reduction<blockSize>(shared_data);
			if(!threadIdx.x) value[0] = shared_data[0];
		#endif
		if(!threadIdx.x) shared_parameters[3] = value[0]/n_subjects;
		__syncthreads();
		value[1] /= shared_parameters[3];
		value[2] /= shared_parameters[3];
		#pragma unroll
		for(row_index = threadIdx.x ; row_index < pitch; row_index += blockSize){
			const GPU_GWAS_TYPE lambda_minus_one_omega = (eigenvalues[row_index] - 1.f)*omega[blockIdx.x*pitch + row_index];
			value[1] -= lambda_minus_one_omega;
			value[2] += lambda_minus_one_omega*lambda_minus_one_omega;
		}
		#ifdef USE_WARPS
			value[1] = warp_reduction<blockSize>(shared_data, value[1]);
		#endif
	
		#ifdef USE_SHARED_MEMORY
			shared_data[threadIdx.x] = value[1];
			__syncthreads();
			reduction<blockSize>(shared_data);
			if(!threadIdx.x) value[1] = shared_data[0];
		#endif	
	
		#ifdef USE_WARPS
			value[2] = warp_reduction<blockSize>(shared_data, value[2]);
		#endif
	
		#ifdef USE_SHARED_MEMORY
			shared_data[threadIdx.x] = value[2];
			__syncthreads();
			reduction<blockSize>(shared_data);
			if(!threadIdx.x)value[2] = shared_data[0];
		#endif			


		if(!threadIdx.x){
			const GPU_GWAS_TYPE dconstraint = calculate_dconstraint(parameter_t);
		
			parameter_t -= dconstraint*value[1]/( calculate_ddconstraint(parameter_t)*value[1] + dconstraint*dconstraint*value[2]);
			next_h2r = calculate_constraint(parameter_t);
			if(parameter_t != parameter_t || fabs(next_h2r - shared_parameters[0]) < precision) continue_loop = false;
		}
          

         	//if(!threadIdx.x) next_h2r = calculate_constraint(parameter_t);
	     	__syncthreads();

	}while(continue_loop && ++iteration_count < 50);
	if(next_h2r == next_h2r && iteration_count < 50){
		if(!threadIdx.x){
			parameters[blockIdx.x*5] = shared_parameters[0];
			parameters[blockIdx.x*5 + 1] = shared_parameters[1];
			parameters[blockIdx.x*5 + 2] = sqrtf(shared_parameters[3])*shared_parameters[2];
			parameters[blockIdx.x*5 + 3] = shared_parameters[3];
		}
	
		

		/*
		value[0] = 0.f;
		
		#pragma unroll
		for(row_index = threadIdx.x; row_index < end; row_index += blockSize){
			value[0] -= logf(omega[blockIdx.x*pitch + row_index]);
		}
		value[0] = (row_index < n_subjects) ? value[0] - logf(omega[blockIdx.x*pitch + row_index]) : value[0];
		#ifdef USE_WARPS
			value[0] = warp_reduction<blockSize>(shared_data,value[0]);
		#endif
		
		#ifdef USE_SHARED_MEMORY
			shared_data[threadIdx.x] = value[0]
			__syncthreads();
			reduction<blockSize>(shared_data);
			if(!threadIdx.x)value[0]= shared_data[0];
		#endif
	
	
		if(!threadIdx.x) parameters[blockIdx.x*5 + 4] = -0.5f*(value[0] + n_subjects*(1.f + logf(shared_parameters[3])));*/		
	}else if(iteration_count < 50){
		 if(!threadIdx.x) parameters[blockIdx.x*5] = nanf("");
		 return;
	}else{
		 if(!threadIdx.x) parameters[blockIdx.x*5] = 99;
		return;
	}
	value[0] = 0.f;
		
	#pragma unroll
	for(row_index = threadIdx.x; row_index < end; row_index += blockSize){
		value[0] -= logf(omega[blockIdx.x*pitch + row_index]);
	}
	value[0] = (row_index < n_subjects) ? value[0] - logf(omega[blockIdx.x*pitch + row_index]) : value[0];
	#ifdef USE_WARPS
		value[0] = warp_reduction<blockSize>(shared_data,value[0]);
	#endif
		
	#ifdef USE_SHARED_MEMORY
		shared_data[threadIdx.x] = value[0];
		__syncthreads();
		reduction<blockSize>(shared_data);
		if(!threadIdx.x)value[0]= shared_data[0];
	#endif
	
	
	if(!threadIdx.x) parameters[blockIdx.x*5 + 4] = -0.5f*(value[0] + n_subjects*(1.f + logf(shared_parameters[3])));						 
}


template< int blockSize>
static void gpu_gwas_call_main_kernel(GPU_Data * const gpu_data, GPU_GWAS_Estimator_Context * const context, GPU_GWAS_Estimator_Data * const results,const GPU_GWAS_TYPE * const snp_data){
	initialize_parameters< blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(context->raw_beta_components_B_C_E, \ 
						  context->Y_component_D_mean_column_component_A, snp_data, context->n_subjects, context->pitch);
	
	maximize_h2r<blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(results->parameters,  context->omega, \
					 context->raw_beta_components_B_C_E, \
				        context->Y_component_D_mean_column_component_A, context->eigenvalues, snp_data, context->n_subjects, context->pitch, powf(10, -context->precision));
	//calculate_loglik< blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(context->omega, results->parameters,  context->n_subjects, context->pitch);
}

template<int blockSize>
static __global__ void calculate_SNP_means(GPU_GWAS_TYPE * __restrict__ const snp_data, const int n_subjects, const int pitch){
	#ifdef USE_WARPS
	__shared__ GPU_GWAS_TYPE shared_data[blockSize/32];
	#endif
	#ifdef USE_SHARED_MEMORY
	__shared__ GPU_GWAS_TYPE shared_data[blockSize];
	#endif
	__shared__ unsigned missing_subject_count;
	__shared__ GPU_GWAS_TYPE shared_mean;
	if(threadIdx.x == 0) missing_subject_count = 0;
	__syncthreads();
	//float * l_SNP_values = &SNP_values[blockIdx.x*blockSize*scale];
	int row_index;
	GPU_GWAS_TYPE sum = 0;
	GPU_GWAS_TYPE value = 0;
#pragma unroll
	for(row_index = threadIdx.x; row_index < pitch; row_index += blockSize){
		value = snp_data[blockIdx.x*pitch + row_index];
		if(value != 3) sum += value;
		if(value == 3) atomicInc(&missing_subject_count, n_subjects);
	}
	//shared_data[threadIdx.x] = sum;
	//__syncthreads();
	//if(missing_subject_count == 0) return;
	#ifdef USE_WARPS
	sum = warp_reduction< blockSize>(shared_data, sum);
	#endif
	#ifdef USE_SHARED_MEMORY
	shared_data[threadIdx.x] = sum;
	__syncthreads();
	reduction<blockSize>(shared_data);
	sum = shared_data[0];
	#endif
       // sum = warp_reduction<blockSize>(sum);
	if(threadIdx.x == 0) shared_mean = sum/(n_subjects - missing_subject_count);
	__syncthreads();
	const int end = (int)floorf((float)n_subjects/blockSize)*blockSize;
#pragma unroll
	for(row_index = threadIdx.x; row_index < end; row_index += blockSize){
		value = snp_data[blockIdx.x*pitch + row_index];
		if(value != 3) snp_data[blockIdx.x*pitch + row_index] -= shared_mean;
		if(value == 3)  snp_data[blockIdx.x*pitch + row_index] = 0;
	}
	if(row_index < pitch) value = snp_data[blockIdx.x*pitch + row_index];
	if(row_index < n_subjects && value != 3) snp_data[blockIdx.x*pitch + row_index] -= shared_mean;
	if(row_index < n_subjects && value == 3) snp_data[blockIdx.x*pitch + row_index] = 0;
	
}

template< int blockSize>
void demean_and_multiply_snps_wrapper(GPU_Data * const gpu_data, GPU_GWAS_Estimator_SNP_Data * const snp_data){
	

        calculate_SNP_means< blockSize><<<gpu_data->gridSize(), gpu_data->blockSize(), 0, gpu_data->stream>>>(snp_data->temp_SNP_data, snp_data->n_subjects,snp_data->pitch);
	//cudaErrorCheck(cudaStreamSynchronize(gpu_data->stream));
	//cudaErrorCheck(cudaDeviceSynchronize());
	GPU_GWAS_TYPE alpha = 1;
	GPU_GWAS_TYPE beta = 0;
	//if(typeid(alpha) == typeid(double())){
	//  cublasErrorCheck(cublasDgemm(gpu_data->cublas_handle, CUBLAS_OP_N,CUBLAS_OP_N,\
					snp_data->pitch, gpu_data->n_data_sets,snp_data->pitch,(double*) &alpha,\
					(const double * const)snp_data->eigenvectors_transposed,snp_data->pitch,\
					(double * )snp_data->temp_SNP_data, snp_data->pitch,\
					(double *)&beta, snp_data->SNP_data, snp_data->pitch));	
	//cudaErrorCheck(cudaStreamSynchronize(gpu_data->stream));
	//}else{
	   cublasErrorCheck(cublasSgemm(gpu_data->cublas_handle, CUBLAS_OP_N,CUBLAS_OP_N,\
					snp_data->pitch, gpu_data->n_data_sets,snp_data->pitch,(float*) &alpha,\
					(const float * const)snp_data->eigenvectors_transposed,snp_data->pitch,\
					(float * )snp_data->temp_SNP_data, snp_data->pitch,\
					(float *)&beta, snp_data->SNP_data, snp_data->pitch));	
	//cudaErrorCheck(cudaDeviceSynchronize());
	//}

}

void GPU_GWAS_Estimator::set_function_pointers(const int _blockSize){
    switch(_blockSize)
    {
        case 32:
            Run_GPU_GWAS_Function_Pointer = &gpu_gwas_call_main_kernel<32>;
            Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer = &demean_and_multiply_snps_wrapper<32>;
            break;
        case 64:
            Run_GPU_GWAS_Function_Pointer = &gpu_gwas_call_main_kernel<64>;
            Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer = &demean_and_multiply_snps_wrapper<64>;
            break;
        case 128:
            Run_GPU_GWAS_Function_Pointer = &gpu_gwas_call_main_kernel<128>;
            Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer = &demean_and_multiply_snps_wrapper<128>;
            break;
        case 256:
            Run_GPU_GWAS_Function_Pointer = &gpu_gwas_call_main_kernel<256>;
            Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer = &demean_and_multiply_snps_wrapper<256>;
            break;
            
        case 512:
            Run_GPU_GWAS_Function_Pointer = &gpu_gwas_call_main_kernel<512>;
            Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer = &demean_and_multiply_snps_wrapper<512>;
            break;

	case 1024:
            Run_GPU_GWAS_Function_Pointer = &gpu_gwas_call_main_kernel<1024>;
            Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer = &demean_and_multiply_snps_wrapper<1024>;
            break;
        
       
    }
    
    
}


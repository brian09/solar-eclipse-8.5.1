#include "gpu-pedigree-data.h"
#include "cuda_runtime.h"
#include <cmath>

static const unsigned FULL_MASK = 0xffffffff;
/*__device__ __inline__
static  int create_index(const int i, const int j, const int N)
{
   if (i <= j)
      return i * N - (i - 1) * i / 2 + j - i;
   else
      return j * N - (j - 1) * j / 2 + i - j;
}*/
__host__ __device__ __inline__ 
static void calculate_indices(const long k, const int n, int & i, int & j){
	const long v = 2*n+1;
	j = floor( ( v - sqrt( (double)(v*v - 8*k) ) ) / 2.0 ) ;
	i = k  + j*(j-1)/2  + j - n*j;

}


template<int blockSize>
static __global__ void combine_results(float * const empirical_pedigree, unsigned * const missing_snp_count, const  float * const other_empirical_pedigree, \
					const unsigned * const other_missing_snp_count, const int array_size){
	
	for(int index = threadIdx.x; index < array_size ; index +=blockSize){
		empirical_pedigree[index] += other_empirical_pedigree[index];
		missing_snp_count[index] += other_missing_snp_count[index];
	}
}



template<int snp_stride>
static __global__ void calculate_empirical_pedigree(float * const empirical_pedigree, unsigned * const missing_snp_count, const float * const allele_frequencies, const snp_t * const snp_data,\
						const float alpha, const int n_subjects, const int pitch, const int array_size, const int batch_size){
						//const int2 * const index_map
	if(blockIdx.x*snp_stride >= batch_size) return;
	switch(snp_stride){
		case 1:
		{
		__shared__ float freq;
		__shared__ float variance;
		if(!threadIdx.x) {
			freq = allele_frequencies[blockIdx.x];
			const float abs_alpha = fabsf(alpha);
			variance = (abs_alpha != 1.f) ? pow(2.f*freq*(1.f-freq),abs_alpha) : 2.f*freq*(1.f-freq);
			if(alpha < 0.f) variance = 1.f/variance;
			freq *= 2.f;
		}
		
		__syncthreads();
		if(freq != 2.f && freq != 0.f && freq != -2.f){
			int col_index;
			int row_index = n_subjects;
			snp_t col_value;
			for(int index = threadIdx.x ; index < array_size; index += blockDim.x){
				if(row_index + blockDim.x >= n_subjects){
					calculate_indices(index, n_subjects, row_index, col_index);
					col_value = snp_data[blockIdx.x*pitch + col_index];
				}else{
					row_index += blockDim.x;
				}
				const snp_t row_value = snp_data[blockIdx.x*pitch + row_index];
	
	                      /*
				const int2 indices = index_map[index];
				if(col_index != indices.y){
					col_index = indices.y;
					col_value = snp_data[blockIdx.x*pitch + col_index];
				}
				const int row_value = snp_data[blockIdx.x*pitch + indices.x];*/
				if(row_value != 3 && col_value != 3){
					atomicAdd(&empirical_pedigree[index], (row_value-freq)*(col_value-freq)*variance);
				}else{
					atomicAdd(&missing_snp_count[index], 1);
				}
			}
		}else if (freq == 2.f || freq == 0.f){
			int col_index;
			int row_index = n_subjects;
			snp_t col_value;		
			for(int index = threadIdx.x ; index < array_size; index += blockDim.x){
				if(row_index + blockDim.x >= n_subjects){
					calculate_indices(index, n_subjects, row_index, col_index);
					col_value = snp_data[blockIdx.x*pitch + col_index];
				}else{
					row_index += blockDim.x;
				}
				const snp_t row_value = snp_data[blockIdx.x*pitch + row_index];			
			/*	const int2 indices = index_map[index];
				if(col_index != indices.y){
					col_index = indices.y;
					col_value = snp_data[blockIdx.x*pitch + col_index];
				}
				const int row_value = snp_data[blockIdx.x*pitch + indices.x];*/
				if(row_value == 3 || col_value == 3){
					atomicAdd(&missing_snp_count[index], 1);
				}
			}
		}else{
			for(int index = threadIdx.x ; index < array_size; index += blockDim.x){
				atomicAdd(&missing_snp_count[index], 1);
			}
		}		
			
		break;
		}
		default:
		{
		__shared__ float freq[snp_stride];
		__shared__ float variance[snp_stride];
		//if(blockIdx.x + 1 != gridDim.x || batch_size % snp_stride == 0){
		if(!threadIdx.x) {
			const float abs_alpha = fabsf(alpha);
			#pragma unroll snp_stride
			for(int snp_index = 0; snp_index < snp_stride; snp_index++){
				freq[snp_index] = allele_frequencies[blockIdx.x*snp_stride + snp_index];
				variance[snp_index] = (abs_alpha != 1.f ) ? powf(2.f*freq[snp_index]*(1.f-freq[snp_index]),abs_alpha) : 2.f*freq[snp_index]*(1.f-freq[snp_index]);
				if(alpha < 0.f ) variance[snp_index] = 1.f/variance[snp_index];
				freq[snp_index] *= 2.f;
			}
		}
		__syncthreads();
		int col_index = -1;
		snp_t col_values[snp_stride];
		int row_index = n_subjects;
				
		for(int index = threadIdx.x ; index < array_size; index += blockDim.x){
			int skip_col_read = 1;
			if(row_index + blockDim.x >= n_subjects){
				calculate_indices(index, n_subjects, row_index, col_index);
				//col_value = snp_data[blockIdx.x*pitch + col_index];
				skip_col_read = 0;
			}else{
			        row_index += blockDim.x;
			}
					
			/*
			const int2 indices = index_map[index];
			int skip_col_read = 1;	
			if(col_index != indices.y){
				col_index = indices.y;
				skip_col_read = 0;
			}*/			
			float sum = 0.f;
			int missing_count = 0;
			#pragma unroll snp_stride
			for(int snp_index = 0; snp_index < snp_stride; snp_index++){
				const float local_freq = freq[snp_index];
				if(local_freq == -2.f) {
					++missing_count;
					continue;
				}
				const float local_variance = variance[snp_index];
				//const int2 indices = index_map[index];

				if(!skip_col_read) col_values[snp_index] =  snp_data[(blockIdx.x*snp_stride + snp_index)*pitch + col_index];			
				//const int col_value = snp_data[(blockIdx.x*snp_stride + snp_index)*pitch + indices.y];
				const snp_t row_value = snp_data[(blockIdx.x*snp_stride + snp_index)*pitch + row_index];
				if(row_value != 3 && col_values[snp_index] != 3){
					if(local_freq != 2.f && local_freq != 0.f){ 
						sum += (col_values[snp_index]-local_freq)*(row_value-local_freq)*local_variance;
					}
						
				}else{
					++missing_count;
				}
			
			}
			if(sum != 0.f) atomicAdd(&empirical_pedigree[index], sum);
				//if(sum != 0.f) atomicAdd(&empirical_pedigree[index], 0.5f);
			if(missing_count != 0) atomicAdd(&missing_snp_count[index], missing_count);
		} 			
			
		break;
		}
		
	}		  
		
}
static int blockSize;
static int minGridSize;
template<int snp_stride>
void call_gpu_kernels(GPU_Pedigree_Data * const pedigree_data){
    const int batch_size = pedigree_data->current_batch_size;
	
	        if(batch_size % snp_stride == 0 && batch_size >= snp_stride){
		        if(pedigree_data->subset_size == 0){
			        calculate_empirical_pedigree<snp_stride><<<ceil(batch_size/snp_stride), blockSize, 0, pedigree_data->stream>>>\
				    (pedigree_data->gpu_empirical_pedigree,pedigree_data->gpu_missing_snp_count, pedigree_data->gpu_frequencies,pedigree_data->gpu_snp_data,\		
				    pedigree_data->alpha, pedigree_data->n_subjects, pedigree_data->pitch,\
				    pedigree_data->array_size, batch_size);
							
		        }else{

			        calculate_empirical_pedigree<snp_stride><<<ceil(batch_size/snp_stride), blockSize, 0, pedigree_data->stream>>>\
				    (pedigree_data->gpu_empirical_pedigree,pedigree_data->gpu_missing_snp_count, pedigree_data->gpu_frequencies,pedigree_data->gpu_snp_data,\		
				    pedigree_data->alpha, pedigree_data->subset_size, pedigree_data->pitch,\
				    pedigree_data->array_size, batch_size);
		        }						

	        }else{
		        if(pedigree_data->subset_size == 0){

			        calculate_empirical_pedigree<snp_stride><<<batch_size, blockSize, 0, pedigree_data->stream>>>\
				    (pedigree_data->gpu_empirical_pedigree,pedigree_data->gpu_missing_snp_count, pedigree_data->gpu_frequencies,pedigree_data->gpu_snp_data,\		
				    pedigree_data->alpha, pedigree_data->n_subjects, pedigree_data->pitch,\
				    pedigree_data->array_size, batch_size);				
		        }else{

			        calculate_empirical_pedigree<snp_stride><<<batch_size, blockSize, 0, pedigree_data->stream>>>\
				    (pedigree_data->gpu_empirical_pedigree,pedigree_data->gpu_missing_snp_count, pedigree_data->gpu_frequencies,pedigree_data->gpu_snp_data,\		
				    pedigree_data->alpha, pedigree_data->subset_size, pedigree_data->pitch,\
				    pedigree_data->array_size, batch_size);
		        }	

	        }	
          
	         
}
void GPU_Pedigree_Context::set_blockSize_gpu_kernel_pointer(){
	if(snp_stride != -1){
	
		if(max_batch_size % snp_stride != 0){
			max_batch_size = ceil(max_batch_size/snp_stride)*snp_stride;
		}
		
		switch(snp_stride){
			case 1:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<1>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<1>;
			break;
			case 2:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<2>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels< 2>;
			break;
			case 3:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<3>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<3>;
			break;
			case 4:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<4>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<4>;
			break;
			case 5:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<5>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<5>;
			break;
			case 6:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<6>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels< 6>;
			break;
			case 7:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<7>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<7>;
			break;
			case 8:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<8>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<8>;
			break;
			case 9:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<9>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<9>;
			break;
			case 10:
			if(thread_size == -1){
				cudaOccupancyMaxPotentialBlockSize ( &minGridSize, &blockSize, calculate_empirical_pedigree<10>);
				thread_size = blockSize;
			//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
			}else{
				blockSize = thread_size;
			}
			Call_GPU_Kernels = &call_gpu_kernels<10>;
			break;	
		}
	}else{
		int best_snp_stride;
		int best_blockSize = -1;
		int best_gridSize = -1;
		if(thread_size != -1) best_blockSize = thread_size;
		for(int current_snp_stride = 10; current_snp_stride >= 1; current_snp_stride--){
			
			int current_blockSize;
			int current_gridSize;
			switch(current_snp_stride){
				case 1:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<1>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<1>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 1;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<1>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 1;
					 	}
					 }
					 break;
				case 2:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<2>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<2>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 2;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<2>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 2;
					 	}
					 }
					 break;			
				case 3:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<3>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<3>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 3;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<3>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 3;
					 	}
					 }
					 break;
				case 4:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<4>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<4>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 4;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<4>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 4;
					 	}
					 }
					 break;
				case 5:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<5>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<5>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 5;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<5>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 5;
					 	}
					 }
					 break;
				case 6:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<6>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<6>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 6;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<6>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 6;
					 	}
					 }
					 break;		
				case 7:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<7>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<7>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 7;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<7>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 7;
					 	}
					 }
					 break;	
				case 8:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<8>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<8>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 8;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<8>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 8;
					 	}
					 }
					 break;	
				case 9:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<9>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<9>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 9;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<9>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 9;
					 	}
					 }
					 break;	
				case 10:
					if(thread_size == -1){
						cudaOccupancyMaxPotentialBlockSize ( &current_gridSize, &current_blockSize, calculate_empirical_pedigree<10>);
						cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<10>, current_blockSize, 0);
					//	if(blockSize > array_size) blockSize = ceil(array_size/32.f)*32;
					 	if(current_gridSize > best_gridSize || best_gridSize == -1){
					 		best_gridSize = current_gridSize;
					 		best_blockSize = current_blockSize;
					 		best_snp_stride = 10;
					 	}						
					}else{						
					 	cudaOccupancyMaxActiveBlocksPerMultiprocessor ( &current_gridSize, calculate_empirical_pedigree<10>, thread_size, 0);
					 	if(current_gridSize > best_gridSize){
					 		best_gridSize = current_gridSize;
					 		best_snp_stride = 10;
					 	}
					 }
					 break;
			}
		}
		
		thread_size = best_blockSize;
		blockSize = thread_size;
		snp_stride = best_snp_stride;
		
		if(max_batch_size % snp_stride != 0){
			max_batch_size = ceil(max_batch_size/snp_stride)*snp_stride;
		}
		
								 					 					 					 			 					 			switch(snp_stride){
			case 1:
				Call_GPU_Kernels = &call_gpu_kernels<1>;
				break;
			case 2:
				Call_GPU_Kernels = &call_gpu_kernels<2>;
				break;
			case 3:
				Call_GPU_Kernels = &call_gpu_kernels<3>;
				break;
			case 4:
				Call_GPU_Kernels = &call_gpu_kernels<4>;
				break;
			case 5:
				Call_GPU_Kernels = &call_gpu_kernels<5>;
				break;
			case 6:
				Call_GPU_Kernels = &call_gpu_kernels<6>;
				break;
			case 7:
				Call_GPU_Kernels = &call_gpu_kernels<7>;
				break;
			case 8:
				Call_GPU_Kernels = &call_gpu_kernels<8>;
				break;
			case 9:
				Call_GPU_Kernels = &call_gpu_kernels<9>;
				break;
			case 10:
				Call_GPU_Kernels = &call_gpu_kernels<10>;
				break;
		}		
	}					 						
}


template<int blockSize>
void call_combine_results_kernels(GPU_Pedigree_Data ** pedigree_data, const int n_gpus){
	for(int gpu_index = 1 ; gpu_index < n_gpus; gpu_index++){
		combine_results<blockSize><<<1, blockSize, 0, pedigree_data[0]->stream>>>(pedigree_data[0]->gpu_empirical_pedigree, pedigree_data[0]->gpu_missing_snp_count,  pedigree_data[gpu_index]->gpu_empirical_pedigree, \
					pedigree_data[gpu_index]->gpu_missing_snp_count, pedigree_data[0]->array_size);
	}

}


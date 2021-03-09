#include "gpu-pedigree-data.h"
#include "cuda_runtime.h"
#include <cmath>
#include "cublas_v2.h"
#define METHOD 3
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
static void calculate_indices(const long k, const long n, int & i, int & j){
	const long v = 2*n+1;
	j = floor( ( v - sqrt( (double)(v*v - 8*k) ) ) / 2.0 ) ;
	i = k - n*j + j*(j-1)/2  + j;

}



static __global__ void demean_snp_data(float * __restrict__ snp_data,  const float * __restrict__ const allele_frequencies, const int n_subjects, const int pitch){
    const float mean = 2.f*allele_frequencies[blockIdx.x];
    if(mean != -2.f && mean != 0.f && mean != 2.f){
        const float sd = sqrtf(mean*(1.f-mean/2.f));
        for(int index = threadIdx.x; index < n_subjects; index += blockDim.x){
            const float value = snp_data[blockIdx.x*pitch + index];
            if(value != 3.f){
                snp_data[blockIdx.x*pitch + index] = (value - mean)/sd;
            }else{
                snp_data[blockIdx.x*pitch + index] = 0.f;
            }
        }
    }else if (mean != -2.f){
        for(int index = threadIdx.x ; index < n_subjects; index += blockDim.x){
              const float value = snp_data[blockIdx.x*pitch + index];
              snp_data[blockIdx.x*pitch + index] = 0.f;
        }
    }else{
        for(int index = threadIdx.x ; index < n_subjects; index += blockDim.x){
              snp_data[blockIdx.x*pitch + index] = 0.f;
        }
    }
}
static __global__ void count_snps(unsigned * __restrict__ const snp_count, const bool * __restrict__ const missing_snps, const unsigned n_snps, const int array_size, const int n_subjects, const int pitch){
    int row_index, col_index;
    bool col_value;
    row_index = n_subjects;
  //  const int warp_index = threadIdx.x % 32;
  //  const int warp_id = threadIdx.x / 32;
    for(int index = threadIdx.x; index < array_size; index += blockDim.x){
        if(row_index + blockDim.x >= n_subjects){
            calculate_indices(index, n_subjects, row_index, col_index);
            col_value = missing_snps[blockIdx.x*pitch + col_index];
        }else{
            row_index += blockDim.x;
        }
        if(col_value){
            atomicInc(&snp_count[index], n_snps);
            continue;
        }
        if(missing_snps[blockIdx.x*pitch + row_index]) atomicInc(&snp_count[index], n_snps);
    }
} 
             
    
void call_gpu_kernels(GPU_Pedigree_Data * const pedigree_data){
    const float alpha = 1.f;
    const float beta = 1.f;
    const int n_rows = (pedigree_data->subset_size == 0) ? pedigree_data->n_subjects : pedigree_data->subset_size;
 //   cudaErrorCheck(cudaMemsetAsync(pedigree_data->missing_snps, false, sizeof(bool)*pedigree_data->pitch\
    *pedigree_data->current_batch_size, pedigree_data->stream));
    
    demean_snp_data<<<pedigree_data->current_batch_size, pedigree_data->blockSize, 0, pedigree_data->stream>>>\
    (pedigree_data->gpu_snp_data,  pedigree_data->gpu_frequencies,\ 
    n_rows, pedigree_data->pitch);  
   /*      
    count_snps<<<pedigree_data->current_batch_size, pedigree_data->blockSize, 0, pedigree_data->stream>>>\
    (pedigree_data->gpu_missing_snp_count, pedigree_data->missing_snps, pedigree_data->total_snps,\
    pedigree_data->array_size, n_rows, pedigree_data->pitch);
    */
    cublasErrorCheck(cublasSsyrk(pedigree_data->handle,CUBLAS_FILL_MODE_LOWER, CUBLAS_OP_N,\
                                 n_rows, pedigree_data->current_batch_size,&alpha,\
                                 pedigree_data->gpu_snp_data,pedigree_data->pitch, &beta,\
                                 pedigree_data->gpu_empirical_pedigree,n_rows));
}            
void GPU_Pedigree_Context::set_blockSize_gpu_kernel_pointer(){
    Call_GPU_Kernels = &call_gpu_kernels;
}
/*
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
}*/



/*										
void GPU_Pedigree_Data::set_blockSize_gpu_kernel_pointer(){
	switch(thread_size){
		case 32:
		Call_GPU_Kernels = &call_gpu_kernels<32>;
		break;
		case 64:
		Call_GPU_Kernels = &call_gpu_kernels<64>;
		break;
		case 128:
		Call_GPU_Kernels = &call_gpu_kernels<128>;
		break;
		case 256:
		Call_GPU_Kernels = &call_gpu_kernels<256>;
		break;
		case 512:
		Call_GPU_Kernels = &call_gpu_kernels<512>;
		break;
		case 1024:
		Call_GPU_Kernels = &call_gpu_kernels<1024>;
		break;	
	}			
}*/

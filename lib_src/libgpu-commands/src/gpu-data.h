#ifndef gpu_data_h
#define gpu_data_h
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "gpu-exception.h"
typedef struct GPU_Data{
    const int gpu_id;
    cudaStream_t  stream;
    cublasHandle_t cublas_handle;
    const int thread_count;
    int n_data_sets;
    
    inline dim3 blockSize(){
        //dim3 return_val(thread_count, 1, 1);
        return dim3(thread_count, 1, 1);
    }
    inline dim3 gridSize(){
        //dim3 return_val(n_data_sets, 1, 1);
        return dim3(n_data_sets, 1, 1);
    }
   /* inline dim3 setup_blockSize(){
        return dim3 (thread_count, 1 , 1);
    }
    
    inline dim3 setup_gridSize(){
        
        return dim3(n_data_sets, 1, 1);
    }
    
    inline dim3 dataset_gridSize(){
        //dim3 return_val(ceil(float(n_data_sets)/dataset_thread_count), 1, 1);
        return dim3(ceil(float(n_data_sets)/dataset_thread_count), 1, 1);
    }
    
    inline dim3 dataset_blockSize(){
        //dim3 return_val(ceil(float(n_data_sets)/dataset_thread_count), 1, 1);
        return dim3(dataset_thread_count, 1, 1);
    }
    */
    GPU_Data(const int _gpu_id, const int _thread_count,  \
             const int _n_data_sets):\
    thread_count(_thread_count),gpu_id(_gpu_id){
        //   gpu_id = _gpu_id;
        //  thread_count = _thread_count;
        n_data_sets = _n_data_sets;
        cudaErrorCheck(cudaSetDevice(gpu_id));
        cudaErrorCheck(cudaStreamCreate(&stream));
        cublasErrorCheck(cublasCreate(&cublas_handle));
        cublasErrorCheck(cublasSetStream(cublas_handle, stream));
    }
    ~GPU_Data(){
        cudaErrorCheck(cudaSetDevice(gpu_id));
        cudaErrorCheck(cudaStreamDestroy(stream));
        cublasErrorCheck(cublasDestroy(cublas_handle))
    }
}GPU_Data;
#endif 

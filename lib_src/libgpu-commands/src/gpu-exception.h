//
//  gpu-exception.h
//  
//
//  Created by Brian Donohue on 7/17/18.
//
//

#ifndef gpu_excpetion_h
#define gpu_excpetion_h
#include <string.h>
#include <stdio.h>
#include "cuda_runtime.h"
#include <stdlib.h>
#include "cublas_v2.h"
#include <string>
#include <iostream>
#include <exception>
class GPU_Exception  : virtual public std::exception {
protected:
    std::string error_message;
public:
    GPU_Exception(const std::string & _error_message): error_message(_error_message){
	}
    virtual ~GPU_Exception() throw() {}
    virtual const char * what() const throw()
	{
		return error_message.c_str();
	}
};

#define cudaErrorCheck(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line)
{
    if (code != cudaSuccess)
    {

	std::string error_message = "GPU Error: " + std::string(cudaGetErrorString(code)) + " in file "\
	+ std::string(file) + " at line " + std::to_string(line);
        throw GPU_Exception(error_message);
    }
}

static const char *cublasErrorFromCode(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
            
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
            
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
            
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
            
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
            
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
            
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
            
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    
    return "<unknown>";
}
#define cublasErrorCheck(ans) { cublasAssert((ans), __FILE__, __LINE__); }
inline void cublasAssert(cublasStatus_t error, const char * file, int line){
    if(error != CUBLAS_STATUS_SUCCESS){
        
       
   	std::string error_message = "GPU Error: " + std::string(cublasErrorFromCode(error)) + " in file "\
	+ std::string(file) + " at line " + std::to_string(line);
        throw GPU_Exception(error_message);     
    }
    
}

#endif /* gpu-exception_h */

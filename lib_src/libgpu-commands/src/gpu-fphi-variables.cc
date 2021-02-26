#include "gpu-fphi-variables.h"
#include <iostream>
using namespace std;
void GPU_FPHI_Shared_Variables::allocate_gpu_memory(){
	cudaErrorCheck(cudaMalloc((void**)&eigenvalues, sizeof(float)*pitch));
	cudaErrorCheck(cudaMalloc((void**)&eigenvectors_transposed, sizeof(float)*pitch*pitch));
	if(use_covariates) cudaErrorCheck(cudaMalloc((void**)&hat_matrix, sizeof(float)*pitch*pitch));
}

void GPU_FPHI_Shared_Variables::free_gpu_memory(){
	if(eigenvalues) cudaErrorCheck(cudaFree(eigenvalues));
	eigenvalues = 0;
	if(eigenvectors_transposed) cudaErrorCheck(cudaFree(eigenvectors_transposed));
	eigenvectors_transposed = 0;
	if(use_covariates) {
		if(hat_matrix) cudaErrorCheck(cudaFree(hat_matrix));
		hat_matrix = 0;
	}
}

GPU_FPHI_Shared_Variables::GPU_FPHI_Shared_Variables(const int _gpu_id,const int _pitch,const  int _n_subjects, float * const _eigenvalues, float * const _eigenvectors_transposed, float * _hat_matrix,const float _ZTZI_0, \
				const float _ZTZI_1,const float _ZTZI_2, const float _ZTZI_3):gpu_id(_gpu_id),pitch(_pitch),n_subjects(_n_subjects),\
				ZTZI_0(_ZTZI_0),ZTZI_1(_ZTZI_1),ZTZI_2(_ZTZI_2),ZTZI_3(_ZTZI_3){
	
	//gpu_id = _gpu_id;
	//cudaErrorCheck(cudaDeviceSet(gpu_id));
	//pitch = _pitch;
	//n_subjects = _n_subjects;
	//ZTZI_0 = _ZTZI_0;
	//ZTZI_1 = _ZTZI_1;
	//ZTZI_2 = _ZTZI_2;
	//ZTZI_3 = _ZTZI_3;
	if(_hat_matrix == 0) 
		use_covariates = false;
	else
		use_covariates = true;
	
	allocate_gpu_memory();
	cudaErrorCheck(cudaMemcpy(eigenvalues, _eigenvalues, sizeof(float)*pitch, cudaMemcpyHostToDevice));
	cudaErrorCheck(cudaMemcpy(eigenvectors_transposed, _eigenvectors_transposed, sizeof(float)*pitch*pitch, cudaMemcpyHostToDevice));
	
	if(use_covariates){
		cudaErrorCheck(cudaMemcpy(hat_matrix, _hat_matrix, sizeof(float)*pitch*pitch, cudaMemcpyHostToDevice));
	}
}
GPU_FPHI_Shared_Variables::~GPU_FPHI_Shared_Variables(){
	free_gpu_memory();
}
	
void GPU_FPHI_Stream_Variables::allocate_gpu_memory(){
	
	cudaErrorCheck(cudaMalloc((void**)&trait_matrix, sizeof(float)*pitch*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&temp_trait_matrix, sizeof(float)*pitch*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&sigma_e_and_a, sizeof(float2)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&theta, sizeof(float2)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&sigma, sizeof(float)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&boolean_score, sizeof(bool)*max_batch_size));

}

void GPU_FPHI_Stream_Variables::free_gpu_memory(){
	if(trait_matrix) cudaErrorCheck(cudaFree(trait_matrix));
	trait_matrix = 0;
	if(temp_trait_matrix) cudaErrorCheck(cudaFree(temp_trait_matrix));
	temp_trait_matrix = 0;
	if(sigma_e_and_a) cudaErrorCheck(cudaFree(sigma_e_and_a));
	sigma_e_and_a = 0;
	if(theta) cudaErrorCheck(cudaFree(theta));
	theta = 0;
	if(sigma) cudaErrorCheck(cudaFree(sigma));
	sigma = 0;
	if(boolean_score) cudaErrorCheck(cudaFree(boolean_score));
	boolean_score = 0;		
}
GPU_FPHI_Stream_Variables::GPU_FPHI_Stream_Variables(GPU_Data * const _gpu_data, const int _pitch, const int _n_subjects, const int _max_batch_size):gpu_data(_gpu_data), \
	pitch(_pitch), n_subjects(_n_subjects), max_batch_size(_max_batch_size){
	
	//gpu_data = gpu_data;
	///pitch = _pitch;
	//n_subjects = _n_subjects;
	//max_batch_size = _max_batch_size;
	allocate_gpu_memory();
	cpu_trait_matrix = new float[pitch*max_batch_size];
	memset(cpu_trait_matrix, 0, sizeof(float)*pitch*max_batch_size);
}
GPU_FPHI_Stream_Variables::~GPU_FPHI_Stream_Variables(){
	delete [] cpu_trait_matrix;
	free_gpu_memory();
}
void GPU_FPHI_Stream_Variables::copy_trait_data_to_gpu(const int batch_size){
	gpu_data->n_data_sets = batch_size;
	current_batch_size = batch_size;
	
	
	cudaErrorCheck(cudaMemcpyAsync(trait_matrix, cpu_trait_matrix, sizeof(float)*pitch*batch_size, cudaMemcpyHostToDevice, gpu_data->stream));
}
void GPU_FPHI_Results::allocate_gpu_memory(){
	cudaErrorCheck(cudaMalloc((void**)&h2r, sizeof(float)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&chi_squared, sizeof(float)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&SE, sizeof(float)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&score, sizeof(float)*max_batch_size));

}

void GPU_FPHI_Results::free_gpu_memory(){
	if(h2r) cudaErrorCheck(cudaFree(h2r));
	h2r = 0;
	if(chi_squared) cudaErrorCheck(cudaFree(chi_squared));
	chi_squared = 0;
	if(SE) cudaErrorCheck(cudaFree(SE));
	SE = 0;
	if(score) cudaErrorCheck(cudaFree(score));
	score = 0;		
}			
GPU_FPHI_Results::GPU_FPHI_Results(float * const _cpu_h2r, float * const _cpu_chi_squared, float * const _cpu_SE,const int _max_batch_size):max_batch_size(_max_batch_size),cpu_h2r(_cpu_h2r),cpu_chi_squared(_cpu_chi_squared),cpu_SE(_cpu_SE){
	//cpu_h2r = _cpu_h2r;
	//cpu_chi_squared = _cpu_chi_squared;
	//cpu_SE = _cpu_SE;
	//max_batch_size = _max_batch_size;
	allocate_gpu_memory();
}
GPU_FPHI_Results::~GPU_FPHI_Results(){
	free_gpu_memory();
}

void GPU_FPHI_Results::copy_results_to_cpu(const int start_index, const int batch_size, cudaStream_t stream){
	cudaErrorCheck(cudaMemcpyAsync(cpu_h2r + start_index, h2r, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream));
	cudaErrorCheck(cudaMemcpyAsync(cpu_chi_squared + start_index, chi_squared, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream));
	cudaErrorCheck(cudaMemcpyAsync(cpu_SE + start_index, SE, sizeof(float)*batch_size, cudaMemcpyDeviceToHost, stream));
}

GPU_FPHI_Estimator::~GPU_FPHI_Estimator(){
	


}
GPU_FPHI_Estimator::GPU_FPHI_Estimator( double * const _raw_phenotype_buffer, const  int _n_phenotypes, const int _n_devices, const int _n_streams,\
					const int _blockSize, const int _pitch,const int _n_subjects,const int _max_batch_size,const bool _verbose): raw_phenotype_buffer(_raw_phenotype_buffer),\
					 n_phenotypes(_n_phenotypes), n_devices(_n_devices), n_streams(_n_streams), blockSize(_blockSize), pitch(_pitch), n_subjects(_n_subjects),\
					 max_batch_size(_max_batch_size), verbose(_verbose){
 	//n_streams = _n_streams;
 	//n_devices = _n_devices;
 	//n_phenotypes = _n_phenotypes;
 	n_phenotypes_left = n_phenotypes;
 	//blockSize = _blockSize;
 	//pitch = _pitch;
 	//n_subjects = _n_subjects;
 	//max_batch_size = _max_batch_size;

 	//raw_phenotype_buffer = _raw_phenotype_buffer;   
 	//trait_matrices = new float*[n_devices*n_streams];
 	//cout << "cpu trait matrix: pitch " << pitch << "  max batch size " << max_batch_size << endl; 
 	//for(int i = 0 ; i < n_devices*n_streams; i++){
 	//	cout << "allocation made\n";
 	//	trait_matrices[i] = new float[pitch*max_batch_size];
 	//	memset(trait_matrices[i], 0, sizeof(float)*pitch*max_batch_size);
 	//}
 	
 	set_function_pointer();	
 	current_col_index = 0;	
    			
    			
}
void GPU_FPHI_Estimator::Thread_Launch(const int gpu_index, const int stream_index, GPU_Data * const gpu_data,  GPU_FPHI_Stream_Variables * const stream_variables,  GPU_FPHI_Shared_Variables * const shared_variables, GPU_FPHI_Results * const  results, const char ** error_message, bool & break_loop){

	try{
		cudaErrorCheck(cudaSetDevice(gpu_data->gpu_id));
	}catch(GPU_Exception & e){
		*error_message = e.what();
		break_loop = true;
		return;
	}

	//float * local_trait_matrix = trait_matrices[gpu_index*n_streams + stream_index];
	int start_index;
	int batch_size = 0;
	if(break_loop) return;
	do{
		#pragma omp critical
		{
		  start_index = fill_trait_matrix(stream_variables->cpu_trait_matrix, batch_size);
		}
		
		if(batch_size != 0){
			try{
				stream_variables->copy_trait_data_to_gpu(batch_size);
				Run_GPU_FPHI_Function_Pointer(gpu_data, shared_variables, stream_variables, results);
				results->copy_results_to_cpu(start_index, batch_size, gpu_data->stream);
				cudaErrorCheck(cudaStreamSynchronize(gpu_data->stream));
			}catch(GPU_Exception & e){
				*error_message = e.what();
				break_loop= true;
				return;
			}
		}
		if(verbose){
		#pragma omp critical
			{
				cout << "Number of left to traits compute: " << n_phenotypes_left << " Percent Complete: " << floor(double(n_phenotypes  - n_phenotypes_left)/double(n_phenotypes))*100 << "% \r";
				cout.flush();
			}
		}
		if(break_loop) return;
	
	}while(batch_size != 0);
	if(verbose && gpu_index == 0 && stream_index == 0){
		cout << endl;
	}


}
int GPU_FPHI_Estimator::fill_trait_matrix(float * const cpu_trait_matrix, int & batch_size){
	if(n_phenotypes_left <= 0){
		batch_size = 0;
		 return -1;
	}
	int current_start = current_col_index;
	if(n_phenotypes_left < max_batch_size)
		batch_size = n_phenotypes_left;
	else
		batch_size = max_batch_size;
	//float * current_trait_matrix = trait_matrices[gpu_index*n_streams + stream_index];
	for(int col = 0; col < batch_size; col++){
		for(int row =0 ;row < n_subjects; row++){
			cpu_trait_matrix[col*pitch + row] = (float)raw_phenotype_buffer[(col+current_start)*n_subjects + row];
		}
	}
	  
	current_col_index += batch_size;
	n_phenotypes_left -= batch_size;
	
	return current_start;
}
/*
GPU_FPHI_Trait_Matirx::GPU_FPHI_Trait_Matrix(const int _gpu_id, double * _raw_phenotype_buffer, const int _max_batch_size, const int _pitch, const int _n_subjects){
	raw_phenotype_buffer = _raw_phenotype_buffer;
	gpu_id = _gpu_id;
	max_batch_size = _max_batch_size;
	n_subjects = _n_subjects;
	pitch = _pitch;
	cpu_trait_matrix = new float[pitch*max_batch_size];
	memset(cpu_trait_matrix, 0 , pitch*max_batch_size);
	cudaErrorCheck(cudaSetDevice(gpu_id));
	allocate_gpu_memory();
}

int GPU_FPHI_Trait_Matirx::fill_trait_matrix(int & current_col_index, int & n_phenotypes_left, int & batch_size){
	if(n_phenotypes_left <= 0){
		batch_size = 0;
		return -1;
	}
	int current_start = current_col_index;
	if(n_phenotype_left < max_batch_size){
		batch_size = n_phenotype_left;
	else



	~GPU_FPHI_Trait_Matrix();*/


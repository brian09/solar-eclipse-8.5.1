#include "gpu-pedigree-data.h"
#include <omp.h>
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip> 
void GPU_Pedigree_Data::allocate_gpu_data(){
	int device;
	const int n_rows = (subset_size == 0) ? n_subjects : subset_size; 
	cudaErrorCheck(cudaGetDevice (&device));
	if(device != gpu_id) cudaErrorCheck(cudaSetDevice(gpu_id));
	cudaErrorCheck(cudaMalloc((void**)&gpu_snp_data, sizeof(float)*max_batch_size*pitch));
	/*if(subset_pitch != 0){
		cudaErrorCheck(cudaMalloc((void**)&gpu_snp_data_subset, sizeof(int)*max_batch_size*subset_pitch));
		cudaErrorCheck(cudaMemset(gpu_snp_data_subset, 0, sizeof(int)*max_batch_size*subset_pitch));
	}*/
	cudaErrorCheck(cudaMalloc((void**)&gpu_frequencies, sizeof(float)*max_batch_size));
	cudaErrorCheck(cudaMalloc((void**)&gpu_empirical_pedigree, sizeof(float)*n_rows*n_rows));
//	cudaErrorCheck(cudaMalloc((void**)&missing_snps, sizeof(bool)*pitch*max_batch_size));
//	cudaErrorCheck(cudaMalloc((void**)&gpu_missing_snp_count, sizeof(unsigned)*array_size));
	
	cudaErrorCheck(cudaMemset(gpu_empirical_pedigree, 0, sizeof(float)*n_rows*n_rows));
//	cudaErrorCheck(cudaMemset(gpu_missing_snp_count, 0, sizeof(unsigned)*array_size));
	cudaErrorCheck(cudaStreamCreate(&stream));
	cublasErrorCheck(cublasCreate(&handle));
	cublasErrorCheck(cublasSetStream(handle, stream));
}
void GPU_Pedigree_Data::free_gpu_data(){
	int device;
	cudaErrorCheck(cudaGetDevice (&device));
	if(device != gpu_id) cudaErrorCheck(cudaSetDevice(gpu_id));
	cudaErrorCheck(cudaFree(gpu_snp_data));
	cudaErrorCheck(cudaFree(gpu_frequencies));
	/*if(subset_pitch != 0){
		cudaErrorCheck(cudaFree(gpu_snp_data_subset));
	}*/
	cudaErrorCheck(cudaFree(gpu_empirical_pedigree));
//	cudaErrorCheck(cudaFree(gpu_missing_snp_count));	
//	cudaErrorCheck(cudaFree(missing_snps));
	cublasErrorCheck(cublasDestroy(handle));
	cudaErrorCheck(cudaStreamDestroy(stream));
}
void GPU_Pedigree_Data::copy_snp_data_to_gpu(){
	cudaErrorCheck(cudaMemcpyAsync(gpu_snp_data,cpu_snp_data, sizeof(float)*pitch*current_batch_size, cudaMemcpyHostToDevice, stream));
	cudaErrorCheck(cudaMemcpyAsync(gpu_frequencies,cpu_frequencies, sizeof(float)*current_batch_size, cudaMemcpyHostToDevice, stream));
}
GPU_Pedigree_Data::~GPU_Pedigree_Data(){
	free_gpu_data();
	cudaErrorCheck(cudaFreeHost(cpu_snp_data));
	cudaErrorCheck(cudaFreeHost(cpu_frequencies));
	//delete [] cpu_snp_data;
}
void GPU_Pedigree_Data::enable_peer_access(const int peer_gpu_id){
	int device;
	cudaErrorCheck(cudaGetDevice (&device));
	if(device != gpu_id) cudaErrorCheck(cudaSetDevice(gpu_id));
	cudaErrorCheck(cudaDeviceEnablePeerAccess(peer_gpu_id, 0));
}
void GPU_Pedigree_Data::disable_peer_access(const int peer_gpu_id){
	int device;
	cudaErrorCheck(cudaGetDevice (&device));
	if(device != gpu_id) cudaErrorCheck(cudaSetDevice(gpu_id));
	cudaErrorCheck(cudaDeviceDisablePeerAccess(peer_gpu_id));
}
GPU_Pedigree_Data::GPU_Pedigree_Data(const int _gpu_id, const float _alpha, const int _blockSize,\
					const int _max_batch_size, const int _total_snps, const int _n_subjects, const int _pitch,\ 
					const int _array_size, const int * const _subset_index_map = 0,\
					 const int _subset_size = 0):gpu_id(_gpu_id),alpha(_alpha),blockSize(_blockSize),max_batch_size(_max_batch_size)\
					 ,total_snps(_total_snps), n_subjects(_n_subjects),pitch(_pitch), array_size(_array_size),\
					  subset_index_map(_subset_index_map), subset_size(_subset_size){
					  
	allocate_gpu_data();
	cudaErrorCheck(cudaHostAlloc((void**)&cpu_snp_data,sizeof(float)*max_batch_size*pitch,cudaHostAllocWriteCombined));
	cudaErrorCheck(cudaHostAlloc((void**)&cpu_frequencies,sizeof(float)*max_batch_size,cudaHostAllocWriteCombined));
//	cpu_snp_data = new int[max_batch_size*pitch];
//	cpu_frequencies = new float[max_batch_size];
	memset(cpu_snp_data,0,sizeof(float)*max_batch_size*pitch);
	
	
}


int GPU_Pedigree_Context::fill_snp_data_array(GPU_Pedigree_Data * const stream_data ){
	if(n_snps_left == 0) return 0;
	
	int batch_size = max_batch_size;
	if(batch_size > n_snps_left) batch_size = n_snps_left;
	std::string freq_str;
	if(subset_size == 0){
		for(int col = 0; col < batch_size; col++){
			freq_stream >> freq_str;
			stream_data->cpu_frequencies[col] = stof(freq_str);
			pio_next_row(plink_file, buffer);
			for(int row = 0; row < n_subjects; row++) stream_data->cpu_snp_data[col*pitch + row] = buffer[row];
			//pio_next_row(plink_file,stream_data->cpu_snp_data + col*pitch);
		/*	for(int row = 0; row < n_subjects; row++){
				stream_data->cpu_snp_data[col*pitch + row] = buffer[row];
			}*/
		}
	}else{
		for(int col = 0; col < batch_size; col++){
			freq_stream >> freq_str;
			stream_data->cpu_frequencies[col] = stof(freq_str);		
			pio_next_row(plink_file,buffer);
			for(int row = 0; row < subset_size; row++){
				stream_data->cpu_snp_data[col*pitch + row] = buffer[stream_data->subset_index_map[row]];
			}
		}
	}		

	stream_data->current_batch_size = batch_size;
	n_snps_left -= batch_size;
	return batch_size;
}
GPU_Pedigree_Context::GPU_Pedigree_Context(pio_file_t * const _plink_file, const char * frequency_filename, const int _max_batch_size, const int _thread_size,\
				  const int _n_subjects, const int _snp_stride, const int _pitch, const int _subset_size = 0):plink_file(_plink_file),\
				  max_batch_size(_max_batch_size),thread_size(_thread_size),\
				  n_subjects(_n_subjects),snp_stride(_snp_stride),pitch(_pitch), \
				  subset_size(_subset_size){
				  
	freq_stream.open(frequency_filename);
	std::string first_line;
	std::getline(freq_stream,first_line);			  
	if(subset_size == 0){
		array_size = n_subjects*(n_subjects+1)/2;
	}else{
		array_size = subset_size*(subset_size+1)/2;
	}
    const int n_rows = (subset_size == 0) ? n_subjects : subset_size; 	  
	temp_cpu_empirical_pedigree = new float[n_rows*n_rows];
//	temp_cpu_missing_snp_count = new unsigned[array_size];
	buffer = new snp_t[n_subjects];
	set_blockSize_gpu_kernel_pointer();
}
GPU_Pedigree_Context::~GPU_Pedigree_Context(){
	freq_stream.close();
	delete [] buffer;
	delete [] temp_cpu_empirical_pedigree;
//	delete [] temp_cpu_missing_snp_count;
}
std::string GPU_Pedigree_Context::combine_gpu_results(GPU_Pedigree_Data ** gpu_pedigree_data, float * const cpu_empirical_pedigree, const int n_gpus){
	std::string error_message;
	for(int gpu_index = 0; gpu_index < n_gpus; gpu_index++){
		try{
			copy_results_to_cpu(gpu_pedigree_data[gpu_index]);
		}catch(GPU_Exception & e){
			error_message = "Failed to copy GPU pedigree data to CPU on GPU ID " + std::to_string(gpu_pedigree_data[gpu_index]->gpu_id);
			return error_message;
		}
		int array_index = 0;
		const int n_rows = (subset_size == 0) ? n_subjects : subset_size;
		for(int col = 0 ; col < n_rows; col++){
		    for(int row = col; row < n_rows; row++){
		        cpu_empirical_pedigree[array_index] += temp_cpu_empirical_pedigree[col*n_rows + row];
		      //  cpu_missing_snp_count[array_index] += temp_cpu_missing_snp_count[array_index];
		        array_index++;
		    }
	    }   
	}
	
	
/*	for(int gpu_index = 1; gpu_index < n_gpus; gpu_index++){
		try{
			gpu_pedigree_data[0]->enable_peer_access(gpu_pedigree_data[gpu_index]->gpu_id);
		}catch(GPU_Exception & e){
			error_message = "Failed to enable peer access for GPU ID " + std::to_string(gpu_pedigree_data[gpu_index]->gpu_id);
			return error_message;
		}
	}
	Combine_GPU_Results_Kernel(gpu_pedigree_data, n_gpus);
	try{
		cudaErrorCheck(cudaStreamSynchronize(gpu_pedigree_data[0]->stream));
	}catch(GPU_Exception & e){
		error_message = "Failure at stream synchronization after calling GPU kernels to combine GPU data";
	}
	for(int gpu_index = 1; gpu_index < n_gpus; gpu_index++){
		try{
			gpu_pedigree_data[0]->disable_peer_access(gpu_pedigree_data[gpu_index]->gpu_id);
		}catch(GPU_Exception & e){
			error_message = "Failed to disable peer access for GPU ID " + std::to_string(gpu_pedigree_data[gpu_index]->gpu_id);
			return error_message;
		}
	}*/	
	return error_message;

}

static std::chrono::high_resolution_clock::time_point start;
static int string_length = 0; 
std::string GPU_Pedigree_Context::gpu_thread_pedigree_creation(GPU_Pedigree_Data * gpu_pedigree_data){
	std::string error_message;
	try{
		cudaErrorCheck(cudaSetDevice(gpu_pedigree_data->gpu_id));
	}catch(GPU_Exception & e){
		error_message = "Failed to set device " + std::to_string(gpu_pedigree_data->gpu_id);
		stop_loop = true;
		return error_message;
	}
	int current_batch_size = 0;
	#pragma omp critical
	{
		current_batch_size = fill_snp_data_array(gpu_pedigree_data);
	}	
	while(current_batch_size != 0 && stop_loop == false){

		gpu_pedigree_data->copy_snp_data_to_gpu();
		Call_GPU_Kernels(gpu_pedigree_data);
		try{
			cudaErrorCheck(cudaStreamSynchronize(gpu_pedigree_data->stream));
					
			//n_snps_computed += current_batch_size;
			
		}catch(GPU_Exception & e){
			error_message = "Error occurred at stream synchronization, likely caused by a problem with a GPU kernel";
			stop_loop = true;
			return error_message;
		}
		#pragma omp critical
		{
			n_snps_computed += current_batch_size;
			std::chrono::high_resolution_clock::time_point next = std::chrono::high_resolution_clock::now();
			std::string progress_str = "Progress: " + std::to_string(100.f*n_snps_computed/n_snps) + "% complete amount of time passed " +std::to_string(std::chrono::duration_cast<std::chrono::seconds>( next - start ).count()) + "s\r";
			if(progress_str.length() > string_length){
				string_length = progress_str.length();
			}
			std::cout << std::setw(string_length) << progress_str;//floor(100.0*n_snps_computed/n_snps) << "% complete\r";
			
			std::cout.flush();			
			current_batch_size = fill_snp_data_array(gpu_pedigree_data);
		}		
	}
	
	return error_message;
}
/*
std::string GPU_Pedigree_Context::gpu_thread_pedigree_creation(GPU_Pedigree_Data * gpu_pedigree_data){
	std::string error_message;
	try{
		cudaErrorCheck(cudaSetDevice(gpu_pedigree_data->gpu_id));
	}catch(GPU_Exception & e){
		error_message = "Failed to set device " + std::to_string(gpu_pedigree_data->gpu_id);
		stop_loop = true;
		return error_message;
	}
	while(n_snps_left > 0 && stop_loop == false){
		int current_batch_size;
		#pragma omp critical
		{
			current_batch_size = fill_snp_data_array(gpu_pedigree_data);
		}
		if(current_batch_size == 0) break;
		gpu_pedigree_data->copy_snp_data_to_gpu();
		Call_GPU_Kernels(gpu_pedigree_data);
		try{
			cudaErrorCheck(cudaStreamSynchronize(gpu_pedigree_data->stream));
					
			n_snps_computed += current_batch_size;
		}catch(GPU_Exception & e){
			error_message = "Error occurred at stream synchronization, likely caused by a problem with a GPU kernel";
			stop_loop = true;
			return error_message;
		}
	}
	return error_message;
}*/		
std::string GPU_Pedigree_Context::run_pedigree_creation(GPU_Pedigree_Data ** gpu_pedigree_data, const int n_gpus, const int _n_snps){
	std::string error_message;
	const char * errmsg = 0;
	stop_loop = false;
	n_snps_computed = 0;
	int last_n_snps_left = 0;
	n_snps = _n_snps;
	n_snps_left = n_snps;
	if(n_gpus == 1){
		/*volatile int snp_data_read_status = 0;
		int next_batch_size = 0;
		int current_batch_size = 0;
		omp_set_num_threads(2);
		#pragma omp parallel
		{
			const int thread_index = omp_get_thread_num();
			switch(thread_index){
				case 0:
				{
				int total_read = 0;
				do{
					while(snp_data_read_status && !stop_loop);
					if(stop_loop) break;
					next_batch_size = fill_snp_data_array(gpu_pedigree_data[0]);
					total_read += next_batch_size;
					snp_data_read_status = 1;
					std::cout << "Progress: " << std::setw(14) << floor(100.0*n_snps_computed/n_snps) << "% complete\r";
					std::cout.flush();
				}while(total_read < n_snps && !stop_loop);
				break;
				}
				case 1:
				{
				while(n_snps_computed < n_snps && !stop_loop){
					while(!snp_data_read_status && !stop_loop);
					if(stop_loop) break;
					current_batch_size = next_batch_size;
					try{
						cudaErrorCheck(cudaMemcpy(gpu_pedigree_data[0]->gpu_snp_data, gpu_pedigree_data[0]->cpu_snp_data, sizeof(int)*gpu_pedigree_data[0]->pitch*current_batch_size, cudaMemcpyHostToDevice));
					}catch(GPU_Exception & e){
						error_message = "Error copying SNP data to GPU";
						stop_loop = true;
						break;
					}
					snp_data_read_status = 0;
					Call_GPU_Kernels(gpu_pedigree_data[0]);
					try{
						cudaErrorCheck(cudaStreamSynchronize(gpu_pedigree_data[0]->stream));
						n_snps_computed += current_batch_size;
			//n_snps_computed += current_batch_size;
			
					}catch(GPU_Exception & e){
						error_message = "Error occurred at stream synchronization, likely caused by a problem with a GPU kernel";
						stop_loop = true;
					}
				}
				break;
				}
			}
		}*/							
		start = std::chrono::high_resolution_clock::now();		
		std::string local_error_message = gpu_thread_pedigree_creation(gpu_pedigree_data[0]);
		if(error_message.length() == 0 && local_error_message.length() != 0){
			error_message = local_error_message;
		}
	}else{
				
		omp_set_num_threads(n_gpus);
		start = std::chrono::high_resolution_clock::now();		
		#pragma omp parallel
		{
			const int gpu_index =omp_get_thread_num();
			std::string local_error_message = gpu_thread_pedigree_creation(gpu_pedigree_data[gpu_index]);
			if(error_message.length() == 0 && local_error_message.length() != 0){
				error_message = local_error_message;
			}
		}
	}
	if(errmsg && error_message.length()== 0){
		error_message = std::string(errmsg);
	}	
							
	std::cout << std::setw(string_length + 1) << "\n";
	//std::cout.flush();
	//std::cout << std::endl;
	return error_message;	
}			
			

					
/*			
std::string GPU_Pedigree_Context::run_pedigree_creation(GPU_Pedigree_Data ** gpu_pedigree_data, const int n_gpus, const int _n_snps){
	std::string error_message;
	const char * errmsg = 0;
	stop_loop = false;
	n_snps_computed = 0;
	if(n_gpus == 1){
		GPU_Pedigree_Data * local_gpu_pedigree_data = gpu_pedigree_data[0];
		int last_n_snps_left = 0;
		n_snps = _n_snps;
		n_snps_left = n_snps;
		
		std::chrono::steady_clock::time_point t1;
		std::chrono::steady_clock::time_point t2;
	
		std::chrono::duration<double> time_span;
		double time_diff = 0;
		double time_diff_sum = 0.0;
		double time_diff_average = 0;
		unsigned time_diffs_measured = 0;
		unsigned total_n_snps_left = n_snps;
		unsigned seconds = 0;
		unsigned minutes = 0;
		unsigned hours = 0;
		unsigned days = 0;
		unsigned max_string_width = 0;
		std::string time_string;
		cudaErrorCheck(cudaSetDevice(local_gpu_pedigree_data->gpu_id));
		omp_set_num_threads(3);
		int current_batch_size;
		volatile bool read_snp_data = true;
		volatile bool copy_snp_data_to_gpu = true;
		volatile bool run_gpu_calculations = false;
	#pragma omp parallel 
		{
			const int thread_index =omp_get_thread_num();

			switch(thread_index){
				case 0:
				while(n_snps_left > 0 && stop_loop == false){
					//std::cout << "case 0 next iter\n";
					while(!read_snp_data){
					}
					current_batch_size = fill_snp_data_array(local_gpu_pedigree_data);

					read_snp_data = false;
					while(!copy_snp_data_to_gpu){
					}
					try{
						cudaErrorCheck(cudaMemcpy(local_gpu_pedigree_data->gpu_snp_data, local_gpu_pedigree_data->cpu_snp_data, sizeof(int)*current_batch_size*local_gpu_pedigree_data->pitch, cudaMemcpyHostToDevice));
					
					}catch(GPU_Exception & e){
					
						errmsg = "Error occurred when copy SNP data to GPU";
						stop_loop = true;
					}
					if(stop_loop) break;
					copy_snp_data_to_gpu = false;
					run_gpu_calculations = true;
				
				}
				break;
			
				case 1:
				while(n_snps_computed != n_snps){
					//std::cout << "case 1 next iter\n";
					while(!run_gpu_calculations){
					}
					const int local_batch_size = current_batch_size;
				
					if(local_batch_size == 0) break;
				
					Call_GPU_Kernels(local_gpu_pedigree_data);
				
					read_snp_data = true;
					try{
						cudaErrorCheck(cudaStreamSynchronize(local_gpu_pedigree_data->stream));
					
						n_snps_computed += local_batch_size;
					}catch(GPU_Exception & e){
						errmsg = "Error occurred at stream synchronization, likely caused by a problem with a GPU kernel";
						stop_loop = true;
					}
					if(stop_loop) break;
					run_gpu_calculations = false;
					copy_snp_data_to_gpu = true;
				}
				break;
			
				case 2:
				bool time_diff_created = false;
			
				while(n_snps_computed != n_snps && stop_loop == false){
					//std::cout << "case 2 next iter \n";
					if(time_diff_created == false){
						int last_n_snps_computed;
						while(n_snps_computed != 0){
						}
						last_n_snps_computed = n_snps_computed;
						t1 = std::chrono::steady_clock::now();
						int diff = 0;
						while( diff == 0 && n_snps_computed != n_snps){
							diff = n_snps_computed - last_n_snps_computed;
						}
						t2 = std::chrono::steady_clock::now();
						time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);	
						time_diff = diff/time_span.count();
						time_diff_created = true;
					}
					days = hours = minutes = 0;
					const int current_snps_left = n_snps - n_snps_computed;
					seconds = ceil(current_snps_left/time_diff); 
					if(seconds >= 60){
						minutes = floor(seconds/60.f);
						if(minutes >= 60){
							hours = floor(minutes/60.f);
							if(hours >= 24){
								days = floor(hours/24.f);
								time_string = std::to_string(days) + " days " + std::to_string(hours % 24) + " hrs " + std::to_string(minutes % 60) + " mins " + std::to_string(seconds % 60) + " secs";
							}else{
								time_string = std::to_string(hours) + " hrs " + std::to_string(minutes % 60) + " mins " + std::to_string(seconds % 60) + " secs";
							}
						}else{
							time_string = std::to_string(minutes) + " mins " + std::to_string(seconds % 60) + " secs";
						}
					}else{
						time_string = std::to_string(seconds) + " secs";
					}					
					std::string output_string = std::to_string(current_snps_left) + " SNPs Left, Estimated Time Until Completion: " + time_string;
					if(max_string_width == 0) max_string_width = output_string.length();
					if(max_string_width < output_string.length()) max_string_width = output_string.length();
					std::cout << std::setw(max_string_width);
					std::cout << output_string << "\r";
					std::cout.flush();
				}
				break;				
				
		
			}
		}
	}else{
		
		
		int last_n_snps_left = 0;
		n_snps = _n_snps;
		n_snps_left = n_snps;
		
		std::chrono::steady_clock::time_point t1;
		std::chrono::steady_clock::time_point t2;
	
		std::chrono::duration<double> time_span;
		double time_diff = 0;
		double time_diff_sum = 0.0;
		double time_diff_average = 0;
		unsigned time_diffs_measured = 0;
		unsigned total_n_snps_left = n_snps;
		unsigned seconds = 0;
		unsigned minutes = 0;
		unsigned hours = 0;
		unsigned days = 0;
		unsigned max_string_width = 0;
		std::string time_string;	
		omp_set_num_threads(n_gpus + 1);
		#pragma omp parallel
		{
			const int thread_index =omp_get_thread_num();
			switch(thread_index){
				case 0:
				{
				bool time_diff_created = false;
			
				while(n_snps_computed != n_snps && stop_loop == false){
				
					if(time_diff_created == false){
						int last_n_snps_computed;
						while(n_snps_computed != 0){
						}
						last_n_snps_computed = n_snps_computed;
						t1 = std::chrono::steady_clock::now();
						int diff = 0;
						while( diff < n_gpus*max_batch_size && n_snps_computed != n_snps){
							diff = n_snps_computed - last_n_snps_computed;
						}
						t2 = std::chrono::steady_clock::now();
						time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);	
						time_diff = diff/time_span.count();
						time_diff_created = true;
					}
					days = hours = minutes = 0;
					const int current_snps_left = n_snps - n_snps_computed;
					seconds = ceil(current_snps_left/time_diff); 
					if(seconds >= 60){
						minutes = floor(seconds/60.f);
						if(minutes >= 60){
							hours = floor(minutes/60.f);
							if(hours >= 24){
								days = floor(hours/24.f);
								time_string = std::to_string(days) + " days " + std::to_string(hours % 24) + " hrs " + std::to_string(minutes % 60) + " mins " + std::to_string(seconds % 60) + " secs";
							}else{
								time_string = std::to_string(hours) + " hrs " + std::to_string(minutes % 60) + " mins " + std::to_string(seconds % 60) + " secs";
							}
						}else{
							time_string = std::to_string(minutes) + " mins " + std::to_string(seconds % 60) + " secs";
						}
					}else{
						time_string = std::to_string(seconds) + " secs";
					}					
					std::string output_string = std::to_string(current_snps_left) + " SNPs Left, Estimated Time Until Completion: " + time_string;
					if(max_string_width == 0) max_string_width = output_string.length();
					if(max_string_width < output_string.length()) max_string_width = output_string.length();
					std::cout << std::setw(max_string_width);
					std::cout << output_string << "\r";
					std::cout.flush();
				}
				}
				break;
			
				default:
				;
				const int gpu_index = thread_index - 1;
				std::string local_error_message = gpu_thread_pedigree_creation(gpu_pedigree_data[gpu_index]);
				if(error_message.length() == 0 && local_error_message.length() != 0){
					error_message = local_error_message;
				}
				break;
				
			}
		}
	}
	if(errmsg && error_message.length()== 0){
		error_message = std::string(errmsg);
	}	
							
	std::cout << "                                                                                     \r";
	std::cout.flush();
	std::cout << std::endl;
	return error_message;
}*/

/*
const char * GPU_Pedigree_Data::run_pedigree_creation(GPU_Stream_Pedigree_Data * stream_data, const int n_streams, const int _n_snps){
	const char * errmsg = 0;
	bool stop_loop = false;
	int last_n_snps_left = 0;
	n_snps = _n_snps;
	n_snps_left = n_snps;
	std::chrono::steady_clock::time_point t1;
	std::chrono::steady_clock::time_point t2;
	
	std::chrono::duration<double> time_span;
	double time_diff = 0;
	double time_diff_sum = 0.0;
	double time_diff_average = 0;
	unsigned time_diffs_measured = 0;
	unsigned total_n_snps_left = n_snps;
	//double seconds = 0;
	//double minutes = 0;
	///double hours = 0;
	//double days = 0;
		//double time_diff = 0;
	unsigned seconds = 0;
	unsigned minutes = 0;
	unsigned hours = 0;
	unsigned days = 0;
	unsigned max_string_width = 0;
	std::string time_string;
	try{
		cudaErrorCheck(cudaMemset(gpu_empirical_pedigree, 0, sizeof(float)*array_size));
	}catch(GPU_Exception & e){
		return "Failed to initialize GPU empirical pedigree values to zero";
	}
	try{
		cudaErrorCheck(cudaMemset(gpu_missing_snp_count, 0, sizeof(unsigned)*array_size));
	}catch(GPU_Exception & e){
		return "Failed to initialize GPU missing SNP count values to zero";
	}
	omp_set_num_threads(n_streams);
#pragma omp parallel 
	{
		const int thread_index =omp_get_thread_num();
		while(n_snps_left != 0 && !stop_loop){
			int batch_size;
			
			#pragma omp critical
			{
				batch_size = fill_snp_data_array(stream_data[thread_index]);
			}
			if(batch_size == 0 || stop_loop) break;
			
			
			try{
				stream_data[thread_index]->copy_snp_data_to_gpu();
			}catch(GPU_Exception & e){
				std::string error_message = "Error occurred when copying SNP data to GPU on stream ID " + std::to_string(thread_index);
				errmsg = error_message.c_str();
				stop_loop = true;
			}
			if(stop_loop) break;
			Call_GPU_Kernels(stream_data[thread_index], gpu_empirical_pedigree, gpu_missing_snp_count);
			try{
				cudaErrorCheck(cudaStreamSynchronize(stream_data[thread_index]->stream));
			}catch(GPU_Exception & e){
				std::string error_message = "Error occurred at stream synchronization on stream ID " + std::to_string(thread_index) + ", likely caused by problem with a GPU kernel";
				errmsg = error_message.c_str();
				stop_loop = true;
			}
			#pragma omp critical
				{
					total_n_snps_left -= batch_size;
				}
			if(stop_loop) break;
			//#pragma omp critical
				if(thread_index == 0){
				const int current_snps_left = total_n_snps_left;
				if(last_n_snps_left != 0){

					
					if(time_diffs_measured != 3) {
						t1 = t2;
						t2 = std::chrono::steady_clock::now();
						time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);					
						time_diff = (last_n_snps_left - current_snps_left)/time_span.count();
						time_diff_sum += time_diff;
						time_diff_average = time_diff_sum/++time_diffs_measured;
					}
					last_n_snps_left = current_snps_left;
					minutes = hours = days = 0;
					seconds = ceil(current_snps_left/time_diff_average);
					if(seconds >= 60){
						minutes = floor(seconds/60.f);
						if(minutes >= 60){
							hours = floor(minutes/60.f);
							if(hours >= 24){
								days = floor(hours/24.f);
								time_string = std::to_string(days) + " days " + std::to_string(hours % 24) + " hrs " + std::to_string(minutes % 60) + " mins " + std::to_string(seconds % 60) + " secs";
							}else{
								time_string = std::to_string(hours) + " hrs " + std::to_string(minutes % 60) + " mins " + std::to_string(seconds % 60) + " secs";
							}
						}else{
							time_string = std::to_string(minutes) + " mins " + std::to_string(seconds % 60) + " secs";
						}
					}else{
						time_string = std::to_string(seconds) + " secs";
					}
				}else{
					t2 = std::chrono::steady_clock::now();
					last_n_snps_left = current_snps_left;
				}
				if(time_string.length() == 0){
					std::cout << current_snps_left << " SNPs Left   \r";
					std::cout.flush();
				}else{
					std::string output_string = std::to_string(current_snps_left) + " SNPs Left, Estimated Time Until Completion: " + time_string;
					if(max_string_width == 0) max_string_width = output_string.length();
					if(max_string_width < output_string.length()) max_string_width = output_string.length();
					std::cout << std::setw(max_string_width);
					std::cout << output_string << "\r";
					std::cout.flush();
				}
				}
		}
	}
	std::cout << "                                                                                     \r";
	std::cout.flush();
	std::cout << std::endl;
	return errmsg;
}*/
void GPU_Pedigree_Context::copy_results_to_cpu(GPU_Pedigree_Data * pedigree_data){
	cudaErrorCheck(cudaSetDevice(pedigree_data->gpu_id));
	const int n_rows = (subset_size == 0) ? n_subjects : subset_size; 
	cudaErrorCheck(cudaMemcpy(temp_cpu_empirical_pedigree, pedigree_data->gpu_empirical_pedigree, sizeof(float)*n_rows*n_rows, cudaMemcpyDeviceToHost));	
//	cudaErrorCheck(cudaMemcpy(temp_cpu_missing_snp_count, pedigree_data->gpu_missing_snp_count, sizeof(unsigned)*array_size, cudaMemcpyDeviceToHost));
}

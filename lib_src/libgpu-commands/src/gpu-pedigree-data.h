#ifndef GPU_PEDIGREE_DATA_H
#define GPU_PEDIGREE_DATA_H
#include "gpu-exception.h"
#include "cuda_runtime.h"
#include "plinkio.h"
#include <fstream>
#include <string>
#include "cublas_v2.h"
class GPU_Pedigree_Data{
private:
	void allocate_gpu_data();
	void free_gpu_data();
public:
	GPU_Pedigree_Data(const int _gpu_id,  const float _alpha, const int _blockSize, const int _max_batch_size,\
				 const int _total_snps, const int _n_subjects, const int _pitch, const int _array_size, \
				 const int * const _subset_index_map, const int _subset_size);
	~GPU_Pedigree_Data();
	const int gpu_id;
	const int max_batch_size;
	const int n_subjects;
	const int pitch;
	const int blockSize;
	int total_snps;
	//const int subset_pitch;
	const int array_size;
	const int subset_size;
	const int * const subset_index_map;
	int current_batch_size;
	float * gpu_frequencies;
	float * cpu_frequencies;
	float * gpu_empirical_pedigree;
//	bool * missing_snps;
//	unsigned * gpu_missing_snp_count;
	const float alpha;
	float * cpu_snp_data;
	float * gpu_snp_data;
//	int * gpu_snp_data_subset;
	cudaStream_t stream;
	cublasHandle_t handle;
	void copy_snp_data_to_gpu();
	void enable_peer_access(const int peer_gpu_id);
	void disable_peer_access(const int peer_gpu_id);

};
class GPU_Pedigree_Context{
private:
   	void (*Call_GPU_Kernels)(GPU_Pedigree_Data * const gpu_pedigree_data);
   	void (*Combine_GPU_Results_Kernel)(GPU_Pedigree_Data ** gpu_pedigree_data, const int _n_gpus);
	int current_snp_data_index;
	int fill_snp_data_array(GPU_Pedigree_Data * const gpu_pedigree_data );
	template<int blockSize>
	void set_gpu_kernel_pointer();
	void set_blockSize_gpu_kernel_pointer();
	std::string gpu_thread_pedigree_creation(GPU_Pedigree_Data * gpu_pedigree_data);
	
	
public:
	snp_t * buffer;
	pio_file_t * const plink_file;
	int n_snps;
	volatile int n_snps_left;
	volatile int n_snps_computed;
	bool stop_loop;
	int array_size;
	int max_batch_size;
	int thread_size;
	const int n_subjects;
	const int subset_size;
	const int pitch;
	int snp_stride;
	float * temp_cpu_empirical_pedigree;
//	unsigned * temp_cpu_missing_snp_count;
	std::ifstream freq_stream;
	
	
	GPU_Pedigree_Context(pio_file_t * const _plink_file,const char * frequency_filename, const int _max_batch_size,\
					 const int _thread_size,  const int _n_subjects, const int _snp_stride, \
					 const int _pitch, const int _subset_size);
	~GPU_Pedigree_Context();
	std::string run_pedigree_creation(GPU_Pedigree_Data ** gpu_pedigree_data, const int _n_gpus, const int _n_snps);
	std::string combine_gpu_results(GPU_Pedigree_Data ** gpu_pedigree_data, float * const cpu_empirical_pedigree, const int n_gpus);
	void copy_results_to_cpu(GPU_Pedigree_Data * gpu_pedigree_data);


};


#endif

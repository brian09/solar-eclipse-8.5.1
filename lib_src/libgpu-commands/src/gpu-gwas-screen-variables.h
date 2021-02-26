#include "gpu-gwas-estimator-context.h"
#include "gpu-exception.h"
#include "cuda_runtime.h"
class GPU_GWAS_Screen_Data{
public:
	int pitch;
	int n_subjects;
	int max_batch_size;
	int current_batch_size;
	float * Sigma_Y;
	float * Sigma;
	float * test_statistic_denom;
	float * chi_squared;
	float * beta_data;
	
	
	float * cpu_chi_squared;
	float * cpu_beta_data;
	
	
	GPU_GWAS_Screen_Data(float * _Sigma_Y, float * _Sigma,   const int _n_subjects , const int _max_batch_size , const int _pitch);
	~GPU_GWAS_Screen_Data();
	void copy_results_to_cpu(const int batch_size, cudaStream_t stream);
	void allocate_gpu_memory();
	void free_gpu_memory();
};

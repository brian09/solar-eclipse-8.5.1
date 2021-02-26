#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "gpu-exception.h"
#include "gpu-data.h"

class GPU_FPHI_Shared_Variables{
private:
	void allocate_gpu_memory();
	void free_gpu_memory();
public:
	const int gpu_id;
	const int pitch;
	const int n_subjects;
	float *  eigenvalues;
	float *  eigenvectors_transposed;
	float * hat_matrix;
	bool use_covariates;
	const float ZTZI_0;
	const float ZTZI_1;
	const float ZTZI_2;
	const float ZTZI_3;
	GPU_FPHI_Shared_Variables(const int _gpu_id,const int _pitch,const int _n_subjects, float * const _eigenvalues, float * const _eigenvectors_transposed, float * hat_matrix,const float _ZTZI_0, \
				const float _ZTZI_1,const float _ZTZI_2, const float _ZTZI_3);
	~GPU_FPHI_Shared_Variables();
	static inline int Static_GPU_Memory_Cost(const int _pitch, const bool _use_covariates){
		int covariate_memory = 0;
		if(_use_covariates){
			covariate_memory = sizeof(float)*_pitch*_pitch;
		}
		return sizeof(float)*(_pitch + _pitch*_pitch) + covariate_memory;
	}
};

class GPU_FPHI_Stream_Variables{
private:
	void allocate_gpu_memory();
	void free_gpu_memory();
public:
	GPU_Data * const gpu_data;
	const int pitch;
	const int n_subjects;
	const int max_batch_size;
	int current_batch_size;
	float * trait_matrix;
	float * temp_trait_matrix;
	float * cpu_trait_matrix;
	float2 * sigma_e_and_a;
	float * sigma;
	float2 * theta;
	bool * boolean_score;
	GPU_FPHI_Stream_Variables(GPU_Data * const  _gpu_data,const int _pitch,const int _n_subjects,const int _max_batch_size);
	void copy_trait_data_to_gpu(const int batch_size);
	~GPU_FPHI_Stream_Variables();
	static inline int Adjustable_GPU_Memory_Cost(const int _pitch){
		return sizeof(float)*(2*_pitch + 1) + sizeof(bool) + sizeof(float2)*2;
	}
};

class GPU_FPHI_Results{
private:
 	void allocate_gpu_memory();
 	void free_gpu_memory();
public:
	const int max_batch_size;
	float * h2r;
	float * chi_squared;
	float * SE;
	float * score;
	
	float * const cpu_h2r;
	float * const cpu_chi_squared;
	float * const cpu_SE;
	
	GPU_FPHI_Results(float * _cpu_h2r, float * _cpu_chi_squared, float * _cpu_SE,const int _max_batch_size);
	~GPU_FPHI_Results();
	void copy_results_to_cpu(const int start_index, const int batch_size,cudaStream_t stream);
	static inline int Adjustable_GPU_Memory_Cost(){
		return sizeof(float)*4;
	}
};
/*
class GPU_FPHI_Trait_Matrix{
private:
	void allocate_gpu_memory();
	void free_gpu_memory();
public:
	int gpu_id;
	double * raw_phenotype_buffer;
	float * cpu_trait_matrix;
	float * trait_matrix;
	float * temp_trait_matrix;
	int max_batch_size;
	int pitch;
	int n_subjects;
	GPU_FPHI_Trait_Matrix(const int _gpu_id, double * _raw_phenotype_buffer, const int _max_batch_size, const int _n_subjects);
	int fill_trait_matrix(int & current_col_index, int & n_phenotypes_left, int & batch_size);
	~GPU_FPHI_Trait_Matrix();
	*/
class GPU_FPHI_Estimator{
private:
    void set_function_pointer();
    void (*Run_GPU_FPHI_Function_Pointer)(GPU_Data * const gpu_data, GPU_FPHI_Shared_Variables * const shared_variables, GPU_FPHI_Stream_Variables * const stream_variables, GPU_FPHI_Results * const results);
public:
    const int n_phenotypes;
    int n_phenotypes_left;
    int current_col_index;
    const int n_devices;
    const int n_streams;
    const int blockSize;
    const int pitch;
    const int n_subjects;
    const int max_batch_size;
    const bool verbose;
    GPU_FPHI_Estimator( double * const _raw_phenotype_buffer, const int _n_phenotypes,const  int _n_devices,const int _n_streams,const int _blockSize,const int _pitch,const int _n_subjects, const int _max_batch_size, const bool _verbose);
    ~GPU_FPHI_Estimator();
    double * const raw_phenotype_buffer;
   
    int fill_trait_matrix( float * const , int & batch_size);
    void Thread_Launch(const int, const int, GPU_Data * const gpu_data,  GPU_FPHI_Stream_Variables * const stream_variables,  GPU_FPHI_Shared_Variables * const shared_variables, GPU_FPHI_Results * const  results, const char ** error_message, bool & break_loop);
};

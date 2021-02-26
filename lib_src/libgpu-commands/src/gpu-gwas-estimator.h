//
//  gpu-gwas-estimator.h
//  
//
//  Created by Brian Donohue on 8/11/18.
//

#ifndef gpu_gwas_estimator_h
#define gpu_gwas_estimator_h
#include "gpu-gwas-estimator-context.h"
#include <vector>
#include <fstream>
#include <string>
#include "gpu-gwas-settings.h"

typedef struct gwas_data{
	double beta;
	double pvalue;
	double SE;
	double chi;
	double SD;
    double h2r;
    double loglik;
	gwas_data & operator = (const gwas_data & var) {
		beta = var.beta;
		pvalue = var.pvalue;
		SE = var.SE;
		chi = var.chi;
		SD= var.SD;
        h2r = var.h2r;
        loglik = var.loglik;
	}
	
}gwas_data;

class GPU_GWAS_Estimator{
private:
    
    void (*Run_GPU_GWAS_Function_Pointer)(GPU_Data * const, GPU_GWAS_Estimator_Context * const , GPU_GWAS_Estimator_Data * const ,\
                                          const GPU_GWAS_TYPE * const);
    void (*Run_GPU_GWAS_Demean_And_Multiply_SNP_Function_Pointer)(GPU_Data * const,GPU_GWAS_Estimator_SNP_Data * const );
    GPU_GWAS_Estimator_Data ** results;
    GPU_GWAS_Estimator_Context ** context;
    GPU_GWAS_Estimator_SNP_Data ** snp_data;
    GPU_Data ** gpu_data;
    std::vector<int> gpu_id_list;
    GPU_GWAS_TYPE * parameters;
  //  GPU_GWAS_TYPE * loglik;

    GPU_GWAS_TYPE **  GPU_Y_component_D_mean_column_component_A;
    //GPU_GWAS_TYPE **  GPU_mean_column;
    GPU_GWAS_TYPE **  GPU_eigenvalues; 
    pio_file_t * plink_file;
  //  GPU_GWAS_TYPE **  GPU_raw_components_A_D;
    std::ofstream  output_stream;
    gwas_data null_result;
    std::string trait_name;
    int max_batch_size;
    int precision;
    int n_subjects;
    int n_snps;
    int pitch;
    int blockSize;
    int n_streams;
    bool open_output_stream;
    bool verbose;
    void set_function_pointers(const int);

    GPU_GWAS_TYPE ** GPU_eigenvectors_transposed;

    void Thread_Launch(const int gpu_index, const int stream_index);
    void allocate_gpu_memory();
    void free_gpu_memory();  
    void write_data_to_file(GPU_GWAS_TYPE * , int , int ); 
    void initialize_output_stream(const char *);
    void close_output_stream();
public:  
     GPU_GWAS_Estimator(std::vector<int> _gpu_id_list, pio_file_t * const _plink_file, const int * const index_map,const GPU_GWAS_TYPE * const CPU_mean_column_component_A,\
                                    const GPU_GWAS_TYPE * const CPU_eigenvectors_transposed, const GPU_GWAS_TYPE * const  CPU_eigenvalues,\
				     int _precision, \
                                    int _max_batch_size, int _blockSize, int total_snps, int _n_subjects,\
				    int _pitch, int _n_streams, const bool _verbose); 
    void run(std::string, const GPU_GWAS_TYPE * const, const char *,gwas_data);
    ~GPU_GWAS_Estimator();
    
    static inline size_t Static_GPU_Memory_Cost(const unsigned _pitch){
        return sizeof(GPU_GWAS_TYPE)*(_pitch*_pitch + _pitch*2);
    }
    	   
};
/*

class Multi_Trait_GPU_GWAS_Estimator:GPU_GWAS_Estimator
{
private:
    unsigned n_snps_computed;
    GPU_GWAS_TYPE * CPU_SNP_Data;
    GPU_GWAS_TYPE ** GPU_SNP_Data;
    GPU_GWAS_TYPE ** CPU_Pinned_SNP_Data;
    void allocate_gpu_memory();
    void free_gpu_memory();
    void Process_SNP_Data_Thread_Launch(const unsigned stream_index, GPU_GWAS_Estimator_SNP_Data * const snp_data);
    void Thread_Launch(const int );
public:
    
    void Process_SNP_Data(pio_file_t * plink_file, const GPU_GWAS_TYPE * const CPU_eigenvectors_transposed, const int * const index_map);
    Multi_Trait_GPU_GWAS_Estimator(std::vector<int> _gpu_id_list,const GPU_GWAS_TYPE * const CPU_eigenvalues,\
                        GPU_GWAS_TYPE * _h2r, GPU_GWAS_TYPE * _loglik, GPU_GWAS_TYPE * _beta_snp, \
                        GPU_GWAS_TYPE * _beta_se,  GPU_GWAS_TYPE * _variance,\
                        int _precision, unsigned _max_batch_size, unsigned \
                        _dataset_thread_count, unsigned _blockSize,\
                       unsigned _total_snps, unsigned _n_subjects, unsigned _pitch, unsigned _n_streams);
    
    
    
    ~Multi_Trait_GPU_GWAS_Estimator();
    void run(const GPU_GWAS_TYPE * const CPU_Y);


    static inline size_t Static_GPU_Memory_Cost(const unsigned _n_subjects, const unsigned _pitch){
        return sizeof(GPU_GWAS_TYPE)*(_n_subjects*3 + _pitch*2);
    }
    static inline size_t Adjustable_GPU_Memory_Cost(const unsigned _pitch){
        return sizeof(GPU_GWAS_TYPE)*_pitch;
    }
};

class Single_Trait_GPU_GWAS_Estimator:GPU_GWAS_Estimator
{
private:

    GPU_GWAS_TYPE ** GPU_eigenvectors_transposed;

    void Thread_Launch(const int gpu_index, const int stream_index);
    
    void allocate_gpu_memory();
    void free_gpu_memory();
    
public:
    
    Single_Trait_GPU_GWAS_Estimator(std::vector<int> _gpu_id_list, pio_file_t * const plink_file, const int * const index_map,
                                    const GPU_GWAS_TYPE * const CPU_Y,const GPU_GWAS_TYPE * const CPU_eigenvectors_transposed, const GPU_GWAS_TYPE * const  CPU_eigenvalues,\
				    GPU_GWAS_TYPE * _parameters, GPU_GWAS_TYPE * _loglik, int _precision, \
                                    int _max_batch_size, int _blockSize, int total_snps, int _n_subjects,\
				    int _pitch, int _n_streams);
    void run();
    ~Single_Trait_GPU_GWAS_Estimator();
    
    static inline size_t Static_GPU_Memory_Cost(const unsigned _pitch){
        return sizeof(GPU_GWAS_TYPE)*(_pitch*_pitch + _pitch*2);
    }
    
    
};

*/

#endif /* gpu_gwas_estimator_h */


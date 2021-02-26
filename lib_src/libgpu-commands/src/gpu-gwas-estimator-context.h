//
//  gpu-gwas-estimator-context.h
//  
//
//  Created by Brian Donohue on 8/18/18.
//

#ifndef gpu_gwas_estimator_context_h
#define gpu_gwas_estimator_context_h

#include <stdio.h>
#include <vector>
#include "plinkio.h"
#include "gpu-exception.h"
#include "gpu-gwas-settings.h"
#include "gpu-data.h"


class GPU_GWAS_Estimator_Data{
private:
    void allocate_gpu_memory();
    void free_gpu_memory();
public:
    GPU_Data *  gpu_data;
    
    int max_batch_size;
    GPU_GWAS_TYPE * parameters;
    //GPU_GWAS_TYPE * loglik;

    
    GPU_GWAS_TYPE * cpu_parameters;
  //  GPU_GWAS_TYPE * cpu_loglik;

    
//    bool is_memory_pinned;
    
    
  //  GPU_GWAS_TYPE * cpu_pinned_parameters;
  //  GPU_GWAS_TYPE * cpu_pinned_loglik;

    
    
    ~GPU_GWAS_Estimator_Data();
    void copy_results_to_cpu(const int);
   // void unpin_cpu_results();
   // void pin_cpu_results(const int , const int);
    GPU_GWAS_Estimator_Data(GPU_Data *  const _gpu_data, \
			    const int _max_batch_size);
    
    static inline size_t Adjustable_GPU_Memory_Cost(){
        return sizeof(GPU_GWAS_TYPE)*5;
    }
};

class GPU_GWAS_Estimator_SNP_Data{
private:
    const int *  index_map;
    
    
    GPU_GWAS_TYPE * CPU_SNP_data;
    
    snp_t * snp_buffer;

    
    void allocate_gpu_memory();
    void free_gpu_memory();
public:
    pio_file_t *   plink_file;
    int  n_snps;
    int start;
    int next_start;
    int next_n_snps;
    int n_subjects;
    int pitch;
    int max_batch_size;
    int blockSize;
    int total_snps;
    GPU_Data *  gpu_data;
    GPU_GWAS_TYPE * SNP_data;
    GPU_GWAS_TYPE * temp_SNP_data;
    const GPU_GWAS_TYPE *  eigenvectors_transposed;
    GPU_GWAS_Estimator_SNP_Data(const GPU_GWAS_TYPE  * const _eigenvectors, pio_file_t * const file,GPU_Data * const _gpu_data, \
                                const int _blockSize, const int _n_subjects,const int _total_snps, const int batch_size,const int _pitch, \
                                const int * const _index_map);
    ~GPU_GWAS_Estimator_SNP_Data();
    inline int get_start() { return start; }
    inline int get_n_snps() { return n_snps; }
    inline int get_next_start() {return next_start; }
    inline int get_next_n_snps() {return next_n_snps; }
    const GPU_GWAS_TYPE * get_SNP_data() const { return SNP_data; }
    int prepare_snp_data();
    void copy_SNP_data_to_gpu();
    static inline size_t Adjustable_GPU_Memory_Cost(const unsigned _pitch){
        return sizeof(GPU_GWAS_TYPE)*2*_pitch;
    }
    
    static inline size_t Shared_GPU_Memory_Cost(const unsigned _pitch){
    	return sizeof(GPU_GWAS_TYPE)*_pitch*_pitch;
    }
    
};

class GPU_GWAS_Estimator_Context{
private:
    void free_gpu_memory();
    void allocate_gpu_memory();
public:
    GPU_Data * gpu_data;
    GPU_GWAS_TYPE * omega;
    const GPU_GWAS_TYPE * Y_component_D_mean_column_component_A;
    //const double * mean_column;
    const GPU_GWAS_TYPE * eigenvalues;
   
    // double *  raw_components_A_D;
   
    GPU_GWAS_TYPE * raw_beta_components_B_C_E;
  
  //  double * components_B_C_E;
  //  GPU_GWAS_TYPE * SNP_data;


    
    //    size_t & n_data_sets;
    int n_subjects;
    int max_batch_size;
    int blockSize;
    int pitch;
    int precision;
    ~GPU_GWAS_Estimator_Context();
    GPU_GWAS_Estimator_Context(GPU_Data * const _gpu_data, GPU_GWAS_TYPE * const _mean_column_component_A, GPU_GWAS_TYPE  * const _eigenvalues,\
                               const int _n_subjects, const int _max_batch_size, const int _blockSize,\
                               const int _pitch, const int _precision);
    
    static inline size_t Adjustable_GPU_Memory_Cost(const int _pitch){
        return sizeof(GPU_GWAS_TYPE)*(_pitch*3 + _pitch);
    }
    static inline size_t Shared_GPU_Memory_Cost(const int _pitch){
    	return sizeof(GPU_GWAS_TYPE)*(5*_pitch);
    }
};
#endif /* gpu_gwas_estimator_context_h */

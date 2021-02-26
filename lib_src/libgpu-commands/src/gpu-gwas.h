//
//  gpu-gwas.h
//  
//
//  Created by Brian Donohue on 8/21/18.
//

#ifndef gpu_gwas_h
#define gpu_gwas_h

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>
#include "cuda_runtime.h"
#include "solar-trait-reader.h"
#include "plinkio.h"
#include "gpu-gwas-settings.h"
/*typedef struct gwas_data{
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
*/

class GPU_GWAS_Exception: virtual public std::exception{
private:
    std::string error_message;
public:
    GPU_GWAS_Exception(const std::string & error): error_message(error){
    }
    virtual ~GPU_GWAS_Exception() throw() {}
    virtual const char * what() const throw()
	{
		return error_message.c_str();
	}
};

class GPU_GWAS{
private:

    
    std::vector<std::string> plink_ids;
    pio_file_t * plink_file;
    Solar_Trait_Reader * reader;
    std::vector<int> gpu_id_list;
    int n_sets;
    int set_index;
    int n_traits;
    int n_subjects;
    int n_snps;
    int max_batch_size;
    int blockSize;
    int scale;
    int dataset_thread_count;
    int precision;
    int pitch;
    int n_streams;
    int defined_blockSize;
    int defined_batch_size;

   // GPU_GWAS_TYPE * loglik;
    bool verbose;
    bool run_multi_trait;
    const char * output_filename;
public:
    int process_next_trait_set();
    std::chrono::milliseconds calibrate_blockSize(const int , int & , double, int );
    GPU_GWAS(std::vector<int> , std::vector<std::string> , const char * , const char *,\
             const char * ,  const int , const int, const int, const int, const bool);
    
    ~GPU_GWAS();
    inline int get_n_sets() {return n_sets; };
    inline int get_current_set_index() {return set_index; };
    
};
#endif /* GPU_GWAS_hpp */

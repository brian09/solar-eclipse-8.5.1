//
//  pedifromsnps.cpp
//  
//
//  Created by Brian Donohue on 8/25/18.
//


#include <cstdlib>
#include <iostream>
#include "plinkio.h"
#include <fstream>
#include "solar.h"
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <mutex>
#include <iomanip>  
#include <new>
#include <omp.h>
#include <forward_list>
#include <iterator>
using namespace std;
//#ifdef USE_FORTRAN
//extern "C" void pedcalc_ (unsigned char * snp_data, float * empirical_pedigree, int * snp_count, int * n_subjects, int * batch_size, int * array_size, float  * alpha);
extern "C" void corrpedcalcwrapperone_ (int *  snp_count,float * empirical_pedigree, float * frequencies,\
unsigned char * snp_data,int * n_subjects,float  * alpha,int * array_size, int * batch_size);
extern "C" void corrpedcalcwrappertwo_ (float *  variance_sum,float * empirical_pedigree, float * frequencies,\
unsigned char * snp_data,int * n_subjects,int * array_size,int * batch_size);
extern "C" void corrpedcalcidlistwrapperone_ (int * snp_count, float * empirical_pedigree,  float * frequencies, unsigned char * snp_data, \
     float * alpha, int * n_subjects, int * array_size,int * max_batch_size);
extern "C" void corrpedcalcidlistwrappertwo_ (float * variance_sum, float * empirical_pedigree, float * frequencies, unsigned char * snp_data, \
       int * n_subjects, int * array_size,int * max_batch_size);     
//extern "C" void kingpedcalcwrapper_(unsigned char * snp_data,int * squared_difference_sum,float * variance_sum,\
     int * n_subjects,int * array_size,int * max_batch_size);
extern "C" void kingpedcalcwrapper_(unsigned char * snp_data,int * squared_difference_sum,int * hetero_col_count,\
     int *hetero_row_count,int * n_subjects,int * array_size,int * max_batch_size);
extern "C" void  kingpedcalcidlistwrapper_(unsigned char * snp_data, int * squared_difference_sum,int * hetero_col_count,\
     		int * hetero_row_count, int * n_subjects, int * array_size, int * batch_size);  
extern "C" void kingpedhomocalcwrapper_(unsigned char * snp_data, int * squared_difference_sum, float * variance_sum, int * n_subjects,\
			 int * array_size, int * batch_size);     		  		      
//#endif 
static inline int create_index(const int i, const int j, const int N)
{
   if (i <= j)
      return i * N - (i - 1) * i / 2 + j - i;
   else
      return j * N - (j - 1) * j / 2 + i - j;
}
static void write_epedigree_to_file(const char * name_buffer, float * epedigree, pio_file_t * plink_file, const int n_subjects,  int * index_map){
    ofstream output_stream(name_buffer);
    
    pio_sample_t * sample_i;
    pio_sample_t * sample_j;
    output_stream << "IDA,IDB,KIN\n"; 
    if(index_map == 0 ){
    	for(int j = 0; j < n_subjects ;j++){
        	sample_j = fam_get_sample(&plink_file->fam_file, j);
        	output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  epedigree[create_index(j,j, n_subjects)] << endl;
        	for(int i = j + 1; i < n_subjects ;i++){
            		sample_i = fam_get_sample(&plink_file->fam_file, i);
			const int index = create_index(i,j, n_subjects);
            		output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  epedigree[index] << endl;
            		output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  epedigree[index] << endl;
        	}
    	}
    }else{
	for(int j = 0 ; j < n_subjects; j++){
		sample_j = fam_get_sample(&plink_file->fam_file, index_map[j]);
		output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  epedigree[create_index(j,j, n_subjects)] << endl;
        	for(int i = j + 1; i < n_subjects ;i++){
            		sample_i = fam_get_sample(&plink_file->fam_file, index_map[i]);
			const int index = create_index(i,j, n_subjects);
            		output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  epedigree[index] << endl;
            		output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  epedigree[index] << endl;
        	}
	}
    }		
    
    output_stream.close();
}
static void write_epedigree_to_file_2(const char * name_buffer, int * epedigree, pio_file_t * plink_file, const int n_subjects,  int * index_map){
    ofstream output_stream(name_buffer);
    
    pio_sample_t * sample_i;
    pio_sample_t * sample_j;
    output_stream << "IDA,IDB,KIN\n"; 
    if(index_map == 0 ){
    	for(int j = 0; j < n_subjects ;j++){
        	sample_j = fam_get_sample(&plink_file->fam_file, j);
        	output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  epedigree[create_index(j,j, n_subjects)] << endl;
        	for(int i = j + 1; i < n_subjects ;i++){
            		sample_i = fam_get_sample(&plink_file->fam_file, i);
			const int index = create_index(i,j, n_subjects);
            		output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  epedigree[index] << endl;
            		output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  epedigree[index] << endl;
        	}
    	}
    }else{
	for(int j = 0 ; j < n_subjects; j++){
		sample_j = fam_get_sample(&plink_file->fam_file, index_map[j]);
		output_stream << sample_j->iid << "," << sample_j->iid <<  "," <<  epedigree[create_index(j,j, n_subjects)] << endl;
        	for(int i = j + 1; i < n_subjects ;i++){
            		sample_i = fam_get_sample(&plink_file->fam_file, index_map[i]);
			const int index = create_index(i,j, n_subjects);
            		output_stream << sample_j->iid << "," << sample_i->iid <<  "," <<  epedigree[index] << endl;
            		output_stream << sample_i->iid << "," << sample_j->iid <<  "," <<  epedigree[index] << endl;
        	}
	}
    }		
    
    output_stream.close();
}
static bool plink_buffer_size_equals_n_subjects;
static mutex mtx; 
static volatile int N_SNPS_LEFT;
static volatile int N_SNPS_COMPUTED = 0 ;

class Empirical_Pedigree_Thread {
private:
	int status = 2;
	int max_batch_size;
	void allocate_memory();
		
	void free_memory();
		
public:
	void set_arrays_to_zero();
	const int get_status(){return status;}
	const char * check_error_status(){
		switch(status){
			case 0:
			return 0;
			case 1:
			return "Failure in allocating memory";
			case 2:
			return "Memory has not been allocated";
			case 3:
			return "Memory has been freed through class function call";		
			case 4:
			return "Runtime error encountered during calculation";
			default:
			return "An unknown error has occurred";
		}
		return 0;
	}	
	
	void flag_runtime_error() { status = 4;}
	const int get_max_batch_size() {return max_batch_size; }

		
	pio_file_t * const plink_file;
	const int * const index_map;
	const int plink_buffer_size;
	const int n_subjects;
	const int array_size;
	const float alpha;
	const bool use_method_one;
	float * frequencies = 0;
	float * empirical_pedigree = 0;
	int * snp_count = 0;
	float * variance_sum = 0;
	
	snp_t * buffer = 0;
	
	void change_max_batch_size(const int new_max_batch_size){
		if(status == 0){
			max_batch_size=new_max_batch_size;
			if(buffer){
				snp_t * new_buffer = 0;
				if(n_subjects == plink_buffer_size){
					new_buffer = new(nothrow) snp_t[max_batch_size*plink_buffer_size];
				}else{
					new_buffer = new(nothrow) snp_t[max_batch_size*n_subjects];
				}
				if(new_buffer){
					delete [] buffer;
					buffer = new_buffer;	
				}else{
					status = 1;
					delete [] buffer;
					buffer = 0;
				}
			}
		}
	}	
	Empirical_Pedigree_Thread(pio_file_t * _plink_file, const int * const _index_map,\
		const int _n_subjects, const int _plink_buffer_size, const float _alpha,\
		const int _max_batch_size, const bool method):plink_file(_plink_file),\
		index_map(_index_map),n_subjects(_n_subjects),\
		plink_buffer_size(_plink_buffer_size), array_size(_n_subjects*(_n_subjects + 1)/2),\
		alpha(_alpha),max_batch_size(_max_batch_size),use_method_one(method){
		allocate_memory();
		if(status == 0){
			set_arrays_to_zero();
		} 
	}
	
	~Empirical_Pedigree_Thread(){
		free_memory();

	}




};
void Empirical_Pedigree_Thread::set_arrays_to_zero(){
	
	if(empirical_pedigree) memset(empirical_pedigree, 0, sizeof(float)*array_size);
	if(snp_count && use_method_one) memset(snp_count, 0, sizeof(int)*array_size);
	if(variance_sum && !use_method_one) memset(variance_sum, 0, sizeof(float)*array_size);
	
}
void Empirical_Pedigree_Thread::allocate_memory(){
	if(!empirical_pedigree){
		
		empirical_pedigree = new(nothrow) float[array_size];
		if(!empirical_pedigree) {
			status = 1;
			return;
		}
	}
	if(!buffer){
		if(n_subjects != plink_buffer_size){
			buffer = new(nothrow) snp_t[max_batch_size*n_subjects];
		}else{
			buffer = new(nothrow) snp_t[max_batch_size*plink_buffer_size];
		}
		if(!buffer){
			status = 1;
			return;
		}
	}
	if(!frequencies){
		frequencies = new(nothrow) float[max_batch_size];
		if(!frequencies){
			status = 1;
			return;
		}
	}	
	if(!snp_count && use_method_one){
		snp_count = new(nothrow) int[array_size];
		if(!snp_count){
			status = 1;
			return;
		}

	}
	if(!variance_sum && !use_method_one){
		variance_sum = new(nothrow) float[array_size];
		if(!variance_sum){
			status = 1;
			return;
		}

	}	
	if(status == 2) status = 0;
}
	
void Empirical_Pedigree_Thread::free_memory(){

	if(frequencies){
		delete [] frequencies;
		frequencies = 0;
	}
	
	if(empirical_pedigree){
		delete [] empirical_pedigree;
		empirical_pedigree = 0;
	}
	if(buffer){
		delete [] buffer;
		buffer = 0;
	}
	if(snp_count && use_method_one){
		delete [] snp_count;
		snp_count = 0;
	}
	if(variance_sum && !use_method_one){
		delete [] variance_sum;
		variance_sum = 0;
	}	
	if(status == 0 ) status = 3;
}
//#endif
class King_Empirical_Pedigree_Thread {
private:
	int status = 2;
	int max_batch_size;
	void allocate_memory();
		
	void free_memory();
		
public:
	void set_arrays_to_zero();
	const int get_status(){return status;}
	const char * check_error_status(){
		switch(status){
			case 0:
			return 0;
			case 1:
			return "Failure in allocating memory";
			case 2:
			return "Memory has not been allocated";
			case 3:
			return "Memory has been freed through class function call";		
			case 4:
			return "Runtime error encountered during calculation";
			default:
			return "An unknown error has occurred";
		}
		return 0;
	}	
	
	void flag_runtime_error() { status = 4;}
	const int get_max_batch_size() {return max_batch_size; }

		
	pio_file_t * const plink_file;
	const int * const index_map;
	const int plink_buffer_size;
	const int n_subjects;
	const int array_size;

	int * hetero_col_count = 0;
	int * hetero_row_count = 0;
//	int * hetero_count_i = 0;
//	int * hetero_count_j = 0;
	int * squared_difference_sum = 0;
	float * empirical_pedigree = 0;
	float * variance_sum = 0;
	snp_t * buffer = 0;
	
	void change_max_batch_size(const int new_max_batch_size){
		if(status == 0){
			max_batch_size=new_max_batch_size;
			if(buffer){
				snp_t * new_buffer = 0;
				if(n_subjects == plink_buffer_size){
					new_buffer = new(nothrow) snp_t[max_batch_size*plink_buffer_size];
				}else{
					new_buffer = new(nothrow) snp_t[max_batch_size*n_subjects];
				}
				if(new_buffer){
					delete [] buffer;
					buffer = new_buffer;	
				}else{
					status = 1;
					delete [] buffer;
					buffer = 0;
				}
			}
		}
	}	
	King_Empirical_Pedigree_Thread(pio_file_t * _plink_file, const int * const _index_map,\
		const int _n_subjects, const int _plink_buffer_size, \
		const int _max_batch_size):plink_file(_plink_file),\
		index_map(_index_map),n_subjects(_n_subjects),\
		plink_buffer_size(_plink_buffer_size), array_size(_n_subjects*(_n_subjects + 1)/2),\
		max_batch_size(_max_batch_size){
		allocate_memory();
		if(status == 0){
			set_arrays_to_zero();
		} 
	}
	
	~King_Empirical_Pedigree_Thread(){
		free_memory();

	}




};
void King_Empirical_Pedigree_Thread::set_arrays_to_zero(){
	
	if(empirical_pedigree) memset(empirical_pedigree, 0, sizeof(float)*array_size);
	if(variance_sum) memset(variance_sum, 0, sizeof(float)*array_size);
	if(squared_difference_sum) memset(squared_difference_sum, 0, sizeof(int)*array_size);
	if(hetero_col_count) memset(hetero_col_count, 0, sizeof(int)*array_size);
	if(hetero_row_count) memset(hetero_row_count, 0, sizeof(int)*array_size);
	
}
void King_Empirical_Pedigree_Thread::allocate_memory(){
	if(!empirical_pedigree){
		
		empirical_pedigree = new(nothrow) float[array_size];
		if(!empirical_pedigree) {
			status = 1;
			return;
		}
	}
	if(!squared_difference_sum){
		
		squared_difference_sum = new(nothrow) int[array_size];
		if(!squared_difference_sum) {
			status = 1;
			return;
		}
	}	
	if(!hetero_row_count){
		
		hetero_row_count = new(nothrow) int[array_size];
		if(!hetero_row_count) {
			status = 1;
			return;
		}
	}
	if(!hetero_col_count){
		
		hetero_col_count = new(nothrow) int[array_size];
		if(!hetero_col_count) {
			status = 1;
			return;
		}
	}
			
	if(!buffer){
		if(n_subjects != plink_buffer_size){
			buffer = new(nothrow) snp_t[max_batch_size*n_subjects];
		}else{
			buffer = new(nothrow) snp_t[max_batch_size*plink_buffer_size];
		}
		if(!buffer){
			status = 1;
			return;
		}
	}
	if(!variance_sum){
		variance_sum = new(nothrow) float[array_size];
		if(!variance_sum){
			status = 1;
			return;
		}

	}
	if(status == 2) status = 0;
}
	
void King_Empirical_Pedigree_Thread::free_memory(){

	if(empirical_pedigree){
		delete [] empirical_pedigree;
		empirical_pedigree = 0;
	}
	if(hetero_col_count){
		delete [] hetero_col_count;
		hetero_col_count =0;
	}
	if(hetero_row_count){
		delete [] hetero_row_count;
		hetero_row_count =0;
	}
	if(squared_difference_sum){
		delete [] squared_difference_sum;
		squared_difference_sum = 0;
	}
		
	if(buffer){
		delete [] buffer;
		buffer = 0;
	}
	if(variance_sum){
		delete [] variance_sum;
		variance_sum = 0;
	}
	if(status == 0 ) status = 3;
}

static int string_length = 0; 
static std::chrono::high_resolution_clock::time_point start;
static void fill_buffer_thread_king(snp_t * buffer,  int & current_batch_size, const int batch_size = 0, \
		 pio_file_t * input_file = 0, const int * _index_map = 0,  const int _map_size = 0, const int total_snps = 0){
	const std::lock_guard<std::mutex> lock(mtx);
	static snp_t * raw_buffer = 0; 

	if(current_batch_size == -1 && raw_buffer){
		delete [] raw_buffer;

		return;
	}
	static volatile int n_snps_computed = 0;
	static pio_file_t * plink_file = 0;
	static int max_batch_size = 0;
	static int plink_buffer_size = 0;
	static int map_size = 0;
	static int n_snps_left;
	static int n_snps;
	static const int * index_map = 0;
	
	if(input_file){
		plink_file = input_file;
		plink_buffer_size = plink_file->bed_file.header.num_samples;
		if(total_snps != 0){
			n_snps = total_snps;
		}else{
			n_snps = plink_file->bed_file.header.num_loci;
		}
		if(_index_map != 0){
			index_map = _index_map;
			map_size  = _map_size;
			if(raw_buffer == 0 ) raw_buffer = new snp_t[plink_buffer_size];
		}
		
		max_batch_size = batch_size;

		
		
		n_snps_computed = 0;
		n_snps_left = N_SNPS_LEFT;
		
		return;
	}
	
		
	if(current_batch_size != 0) n_snps_computed += current_batch_size;
	std::chrono::high_resolution_clock::time_point next = std::chrono::high_resolution_clock::now();
	std::string progress_string = "Progress: " + std::to_string((float)100.f*n_snps_computed/n_snps) + "% complete amount of time passed: " + std::to_string(std::chrono::duration_cast<std::chrono::seconds>( next - start ).count()) + "s\r";
	if(progress_string.length() > string_length){
		string_length = progress_string.length();
	}
	cout << setw(string_length) << progress_string;
	cout.flush();
	current_batch_size = 0;
	if(!N_SNPS_LEFT) return;
	
	current_batch_size = max_batch_size;
	if(N_SNPS_LEFT < current_batch_size){
		current_batch_size = N_SNPS_LEFT;
	}
	if(index_map != 0){	
	
		for(int col = 0; col < current_batch_size; col++){
			pio_next_row(plink_file, raw_buffer);

			for(int row = 0 ; row < map_size; row++){
				buffer[row + col*map_size] = raw_buffer[index_map[row]];
			}
			
		}
	}else{

		for(int col = 0; col < current_batch_size; col++){

			pio_next_row(plink_file, buffer + col*plink_buffer_size);
		}
	}		
	N_SNPS_LEFT -= current_batch_size;

}
static void fill_buffer_thread(snp_t * buffer, float * frequencies, int & current_batch_size, const int batch_size = 0, const char * frequencies_filename= 0,\
		 pio_file_t * input_file = 0, const int * _index_map = 0,  const int _map_size = 0, const int total_snps = 0){
	const std::lock_guard<std::mutex> lock(mtx);
	static snp_t * raw_buffer = 0; 
	static ifstream freq_stream;
	if(current_batch_size == -1 && raw_buffer){
		delete [] raw_buffer;
		if(freq_stream.is_open()){
			freq_stream.close();
		}
		return;
	}
	static volatile int n_snps_computed = 0;
	static pio_file_t * plink_file = 0;
	static int max_batch_size = 0;
	static int plink_buffer_size = 0;
	static int map_size = 0;
	static int n_snps_left;
	static int n_snps;
	static const int * index_map = 0;
	
	if(input_file){
		plink_file = input_file;
		plink_buffer_size = plink_file->bed_file.header.num_samples;
		if(total_snps != 0){
			n_snps = total_snps;
		}else{
			n_snps = plink_file->bed_file.header.num_loci;
		}
		if(_index_map != 0){
			index_map = _index_map;
			map_size  = _map_size;
			if(raw_buffer == 0 ) raw_buffer = new snp_t[plink_buffer_size];
		}
		
		max_batch_size = batch_size;
		if(freq_stream.is_open() == false){
			freq_stream.open (frequencies_filename, std::ifstream::in);
			std::string first_line;
			getline(freq_stream, first_line);
		}
		
		
		n_snps_computed = 0;
		n_snps_left = N_SNPS_LEFT;
		
		return;
	}
	
		
	if(current_batch_size != 0) n_snps_computed += current_batch_size;
	std::chrono::high_resolution_clock::time_point next = std::chrono::high_resolution_clock::now();
	std::string progress_string = "Progress: " + std::to_string((float)100.f*n_snps_computed/n_snps) + "% complete amount of time passed: " + std::to_string(std::chrono::duration_cast<std::chrono::seconds>( next - start ).count()) + "s\r";
	if(progress_string.length() > string_length){
		string_length = progress_string.length();
	}
	cout << setw(string_length) << progress_string;
	cout.flush();
	current_batch_size = 0;
	if(!N_SNPS_LEFT) return;
	
	current_batch_size = max_batch_size;
	if(N_SNPS_LEFT < current_batch_size){
		current_batch_size = N_SNPS_LEFT;
	}
	if(index_map != 0){	
		std::string freq_str;
		for(int col = 0; col < current_batch_size; col++){
			pio_next_row(plink_file, raw_buffer);
			freq_stream >> freq_str;
			frequencies[col] = stof(freq_str);
			for(int row = 0 ; row < map_size; row++){
				buffer[row + col*map_size] = raw_buffer[index_map[row]];
			}
			
		}
	}else{
		std::string freq_str;
		for(int col = 0; col < current_batch_size; col++){
			freq_stream >> freq_str;
			frequencies[col] = stof(freq_str);
			pio_next_row(plink_file, buffer + col*plink_buffer_size);
		}
	}		
	N_SNPS_LEFT -= current_batch_size;

}

extern "C" void fillbufferthread_ (unsigned char * buffer,float * frequencies, int * current_batch_size){


	fill_buffer_thread((snp_t*)buffer, frequencies, *current_batch_size);
	
}
extern "C" void fillbufferthreadking_ (unsigned char * buffer, int * current_batch_size){


	fill_buffer_thread_king((snp_t*)buffer,  *current_batch_size);
	
}

static string calculate_correlation_empirical_pedigree_with_id_list(pio_file_t * plink_file, const char * frequency_filename, const char * output_filename,  int *  index_map, float alpha, const int per_chromosome,\
								 int n_subjects, int batch_size , const int n_threads, const bool normalize, const bool use_method_one){
    string errmsg;
    int plink_buffer_size = plink_file->bed_file.header.num_samples;
    int array_size =    n_subjects*(n_subjects + 1)/2; 
    const int n_snps= plink_file->bed_file.header.num_loci;	
    float * empirical_pedigree = new(nothrow) float[array_size];
    float * variance_sum = 0;
    int * snp_count = 0;
    if(!empirical_pedigree) return string("Failed to allocate memory for storing empirical pedigree");
    memset(empirical_pedigree, 0, sizeof(float)*array_size);
    if(use_method_one){
    	snp_count = new(nothrow) int[array_size];
    	if(!snp_count) return string("Failed to allocate memory for storing empirical pedigree");
    }else{
    	variance_sum = new(nothrow) float[array_size];
    	if(!variance_sum) return string("Failed to allocate memory for storing empirical pedigree");    
    }   
   
    

    Empirical_Pedigree_Thread ** thread_data = new(nothrow) Empirical_Pedigree_Thread*[n_threads];
    if(!thread_data){
    	if(use_method_one)
    		delete [] snp_count;
    	else
    		delete [] variance_sum;
    	delete [] empirical_pedigree;
    	return string("Failed to allocate memory for storing CPU thread data");
    }
    	
    for(int index = 0; index < n_threads; index++){

    	thread_data[index] = new(nothrow) Empirical_Pedigree_Thread(plink_file, index_map, n_subjects,\
    					 plink_buffer_size, alpha, batch_size, use_method_one);
    	
    	if(!thread_data[index] ){
    		delete [] empirical_pedigree;
    		if(use_method_one)
    			delete [] snp_count;
    		else
    			delete [] variance_sum;
    		for(int i =0; i < index; i++) delete thread_data[i];
    		delete [] thread_data;
    		errmsg = "Failed to allocate memory for CPU thread ID: " + to_string(index);
    		return errmsg;
    	}
	const char * error_status = 0;
	if( error_status = thread_data[index]->check_error_status()){
    		delete [] empirical_pedigree;
    		if(use_method_one)
    			delete [] snp_count;
    		else
    			delete [] variance_sum;
    		for(int i =0; i < index; i++) delete thread_data[i];
    		delete [] thread_data;
    		errmsg = "Error: "  + string(error_status) + " CPU Thread ID: " + to_string(index);
    		return errmsg;
	}
    	
     }
 

     if(!per_chromosome){
     
     	N_SNPS_COMPUTED = 0;
    	N_SNPS_LEFT = n_snps;
    	vector<thread> cpu_threads;
	int dummy_var = 0;
	
 	fill_buffer_thread(0, 0, dummy_var,   batch_size , frequency_filename, plink_file, index_map,  n_subjects);
 	if(use_method_one){
 		start = std::chrono::high_resolution_clock::now();
 		for(int thread_index = 0; thread_index < n_threads; thread_index++){
 			cpu_threads.push_back(thread(corrpedcalcidlistwrapperone_,thread_data[thread_index]->snp_count,\
 			thread_data[thread_index]->empirical_pedigree,thread_data[thread_index]->frequencies,thread_data[thread_index]->buffer,\
     		 	 &alpha, &n_subjects, &array_size, &batch_size));
     		}
     		
     	}else{
     		start = std::chrono::high_resolution_clock::now();
 		for(int thread_index = 0; thread_index < n_threads; thread_index++){
 			cpu_threads.push_back(thread(corrpedcalcidlistwrappertwo_,thread_data[thread_index]->variance_sum,\
 			thread_data[thread_index]->empirical_pedigree,thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	 &n_subjects, &array_size, &batch_size));
     		}
     	}     			 	
	/*#pragma omp parallel
	{
		const int thread_index =omp_get_thread_num();
		if(use_method_one){
			corrpedcalcidlistwrapperone_ (thread_data[thread_index]->snp_count,thread_data[thread_index]->empirical_pedigree, thread_data[thread_index]->buffer,\
     		 	 &alpha, index_map, &plink_buffer_size, &n_subjects, &array_size, &batch_size);   
     		 }else{
 			corrpedcalcidlistwrappertwo_(thread_data[thread_index]->variance_sum,thread_data[thread_index]->empirical_pedigree, thread_data[thread_index]->buffer,\
     		 	index_map, &plink_buffer_size, &n_subjects, &array_size, &batch_size); 
     		 }   		 	
     	}*/		
	if(use_method_one){	
		memset(snp_count, 0, sizeof(int)*array_size);
	}else{	
		memset(variance_sum, 0, sizeof(float)*array_size);
	}
     	for(int thread_index = 0; thread_index < n_threads; thread_index++){
     		cpu_threads[thread_index].join();
     		if(use_method_one){
     			for(int index = 0; index < array_size; index++){
     				empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     				snp_count[index] += thread_data[thread_index]->snp_count[index];
     			}
     		}else{
     			for(int index = 0; index < array_size; index++){
     				empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     				variance_sum[index] += thread_data[thread_index]->variance_sum[index];
     			}
     		}     			
     		delete thread_data[thread_index];
     	}

     	delete [] thread_data;
     	for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] /= ((use_method_one) ? snp_count[index] : variance_sum[index]);//snp_count[index];
    	 }
    	if(normalize){
     		float norms[n_subjects];
     		for(int index = 0; index < n_subjects; index++){
     			norms[index] = sqrt(empirical_pedigree[create_index(index,index,n_subjects)]);
     		//empirical_pedigree[create_index(index,index, n_subjects)] = 1.f;
     		}
     		for(int col = 0; col < n_subjects; col++){
     			for(int row = col ; row < n_subjects; row++){
     				empirical_pedigree[create_index(row,col, n_subjects)] /= (norms[row]*norms[col]);
     			}
   		} 
   	}
    	write_epedigree_to_file(output_filename, empirical_pedigree, plink_file,  n_subjects,  0);  
     	string final_message =  "Empirical pedigree creation is complete";
     	cout << final_message << setw(string_length + 1) << "\n";
     	dummy_var = -1;
     	fill_buffer_thread(0, 0, dummy_var);
     	
   
    }else{
 	pio_locus_t * locus;
        locus = bim_get_locus(&plink_file->bim_file, 0);
        unsigned char chromosome = locus->chromosome;
        int snp_start = 0;
        int snp_end = 0;
        int last_batch_size = batch_size;
        while(snp_start < plink_file->bed_file.header.num_loci){
        
            int snp = snp_start;
            while(locus->chromosome == chromosome && snp < plink_file->bed_file.header.num_loci) locus = bim_get_locus(&plink_file->bim_file, snp++);
            snp_end = snp;
            int snp_batch_size = snp_end - snp_start;
            N_SNPS_COMPUTED = 0;
            N_SNPS_LEFT = snp_batch_size;
            int current_batch_size = last_batch_size;
            if(snp_batch_size < n_threads*last_batch_size){
            	current_batch_size = floor(snp_batch_size/n_threads);
            	for(int i = 0; i < n_threads; i++){
            		thread_data[i]->change_max_batch_size(current_batch_size);
            		if(thread_data[i]->get_status()){
            			delete [] empirical_pedigree;
    				if(use_method_one)
    					delete [] snp_count;
    				else
    					delete [] variance_sum;            			
            			const char * error_status_message = thread_data[i]->check_error_status();
            			for(int j = 0; j < n_threads; j++){
            				delete thread_data[j];
            			}
            			delete [] thread_data;
            			errmsg = "Error: "  + string(error_status_message) + " CPU Thread ID: " + to_string(i);
            			return errmsg;
            		}
            	}
            }
            
            for(int index =0; index < n_threads;index++){
            	thread_data[index]->set_arrays_to_zero();
            }
            memset(empirical_pedigree, 0, sizeof(float)*array_size);
	    if(use_method_one){	
		memset(snp_count, 0, sizeof(int)*array_size);
	    }else{	
		memset(variance_sum, 0, sizeof(float)*array_size);
	    }
            cout << "Starting calculation for chromosome: " << (unsigned)chromosome << " containing " <<  snp_batch_size << " loci\n";
            cout << "Total loci left to calculate: " <<   plink_file->bed_file.header.num_loci - snp_start << endl;
            
    	    vector<thread> cpu_threads;
	    int dummy_var = 0;
	    
 	    fill_buffer_thread(0, 0, dummy_var,   batch_size , frequency_filename, plink_file, index_map,  n_subjects, snp_batch_size);
 	    if(use_method_one){
 	    	start = std::chrono::high_resolution_clock::now();
 		for(int thread_index = 0; thread_index < n_threads; thread_index++){
 			cpu_threads.push_back(thread(corrpedcalcidlistwrapperone_,thread_data[thread_index]->snp_count,\
 			thread_data[thread_index]->empirical_pedigree,thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	 &alpha,  &n_subjects, &array_size, &batch_size));
     		}
     	    }else{
     	    	start = std::chrono::high_resolution_clock::now();
 		for(int thread_index = 0; thread_index < n_threads; thread_index++){
 			cpu_threads.push_back(thread(corrpedcalcidlistwrappertwo_,thread_data[thread_index]->variance_sum,\
 			thread_data[thread_index]->empirical_pedigree,thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	 &n_subjects, &array_size, &batch_size));
     		}
     	    }  	              
     	          	   
     	    for(int thread_index = 0; thread_index < n_threads; thread_index++){
     		cpu_threads[thread_index].join();
     		for(int index = 0; index < array_size; index++){
     			empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     			if(use_method_one){
     				snp_count[index] += thread_data[thread_index]->snp_count[index];
     			}else{	
     				variance_sum[index] += thread_data[thread_index]->variance_sum[index];
     			}
     		}
     		
     	    }
     	    
     	    for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] /= ((use_method_one) ? snp_count[index] : variance_sum[index]);//variance_sum_or_snp_count[index];

    	    } 
    	    if(normalize){
     	    	float norms[n_subjects];
     	    	for(int index = 0; index < n_subjects; index++){
     			norms[index] = sqrt(empirical_pedigree[create_index(index,index,n_subjects)]);
     		//empirical_pedigree[create_index(index,index, n_subjects)] = 1.f;
     	    	}
     	    	for(int col = 0; col < n_subjects; col++){
     			for(int row = col ; row < n_subjects; row++){
     				empirical_pedigree[create_index(row,col, n_subjects)] /= (norms[row]*norms[col]);
     			}
   	    	} 
   	    }
   	    string str_output_filename = string(output_filename) + ".chr" + to_string(chromosome) + ".csv";
    	    write_epedigree_to_file(str_output_filename.c_str(), empirical_pedigree, plink_file,  n_subjects,  0);

     	    cout << "Empirical pedigree creation is complete for chromosome " << (unsigned)chromosome  << endl;     	    
     	      
            if(snp_end != plink_file->bed_file.header.num_loci){
               chromosome = locus->chromosome;
            }
            snp_start = snp_end;
            
        }
        for(int index =0 ;index < n_threads; index++){
        	delete thread_data[index];
        }
        delete [] thread_data;
        int dummy_var = -1;
	fill_buffer_thread(0, 0, dummy_var);
        
        
    }
    delete [] empirical_pedigree;
    if(use_method_one)
        delete [] snp_count;
    else
        delete [] variance_sum;
        
    return errmsg;
}

static string calculate_robust_king_empirical_pedigree_id_list(pio_file_t * plink_file, const char * output_filename,  int * index_map, const int per_chromosome,  int n_subjects,  int batch_size, const int n_threads){
	string error_message;
	int plink_buffer_size = plink_file->bed_file.header.num_samples;
	int array_size = n_subjects*(n_subjects + 1)/2;
	int n_snps = plink_file->bed_file.header.num_loci;
	float * empirical_pedigree = new(nothrow) float[array_size];
	int * hetero_row_count = new(nothrow) int[array_size];
	int * hetero_col_count = new(nothrow) int[array_size];
    	if(!empirical_pedigree) return string("Failed to allocate memory for storing empirical pedigree");
    	memset(empirical_pedigree, 0, sizeof(float)*array_size);
    	
//    	float * variance_sum = new(nothrow) float[array_size];
//     	if(!variance_sum) return string("Failed to allocate memory for storing variance sums");
//    	memset(variance_sum, 0, sizeof(float)*array_size);
    	King_Empirical_Pedigree_Thread ** thread_data = new(nothrow) King_Empirical_Pedigree_Thread*[n_threads]; 
     	if(!thread_data){
 //   		delete [] variance_sum;
    		delete [] empirical_pedigree;
    		return string("Failed to allocate memory for storing CPU thread data");
   	 }   	  	
   	
 	for(int index = 0; index < n_threads; index++){
 	
 		thread_data[index] = new(nothrow) King_Empirical_Pedigree_Thread(plink_file, index_map, n_subjects, plink_buffer_size, batch_size);
 		if(!thread_data[index]){
 			delete [] empirical_pedigree;
 			delete [] hetero_row_count;
 			delete [] hetero_col_count;
 			for(int i = 0; i < index; i++) delete thread_data[index];
 			delete [] thread_data;
 			error_message = "Failed to allocate memory for CPU thread ID: " +to_string(index);
 			return error_message;
 		}
 		const char * error_status = 0;
		if( error_status = thread_data[index]->check_error_status()){
    			delete [] empirical_pedigree;
 			delete [] hetero_row_count;
 			delete [] hetero_col_count;
    			for(int i =0; i < index; i++) delete thread_data[i];
    			delete [] thread_data;
    			error_message = "Error: "  + string(error_status) + " CPU Thread ID: " + to_string(index);
    			return error_message;
		}
	}
    		
     if(!per_chromosome){
     
     	N_SNPS_COMPUTED = 0;
    	N_SNPS_LEFT = n_snps;

    	vector<thread> cpu_threads;
	int dummy_var = 0;
 	
 	fill_buffer_thread_king(0, dummy_var,   batch_size , plink_file, index_map,  n_subjects);
 	start = std::chrono::high_resolution_clock::now();
 	for(int thread_index = 0; thread_index < n_threads; thread_index++){
 		cpu_threads.push_back(thread(kingpedcalcidlistwrapper_,(unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count,\
     		thread_data[thread_index]->hetero_row_count, &n_subjects, &array_size, &batch_size));
     	}
     	

    	memset(hetero_row_count, 0, sizeof(int)*array_size);
        memset(hetero_col_count, 0, sizeof(int)*array_size);
 /*	fill_buffer(0, dummy_var,  n_subjects, batch_size, plink_file);  
 	#pragma omp parallel
 	{
 		const int thread_index =omp_get_thread_num(); 
 	     	kingpedcalcidlistwrapper_((unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count,\
     		thread_data[thread_index]->hetero_row_count, index_map,&plink_buffer_size, &n_subjects, &array_size,  &batch_size);      		
     	}      	
         */
     	for(int thread_index = 0; thread_index < n_threads; thread_index++){
     		cpu_threads[thread_index].join();
     		for(int index = 0; index < array_size; index++){
     			empirical_pedigree[index] += thread_data[thread_index]->squared_difference_sum[index];
     			hetero_row_count[index] += thread_data[thread_index]->hetero_row_count[index];
     			hetero_col_count[index] += thread_data[thread_index]->hetero_col_count[index];
     			//variance_sum[index] += thread_data[thread_index]->variance_sum[index];
     		}
     		delete thread_data[thread_index];
     	}
     	delete [] thread_data;
     	for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] = 1.0 - empirical_pedigree[index]/(2.f*min(hetero_row_count[index], hetero_col_count[index]));
    	 }

    	write_epedigree_to_file(output_filename, empirical_pedigree, plink_file,  n_subjects,  0);  
     	
     	cout << "King Empirical pedigree creation is complete"  << "\n";
     	delete [] empirical_pedigree;
 	delete [] hetero_row_count;
 	delete [] hetero_col_count;
 	dummy_var = -1;
 	fill_buffer_thread_king(0, dummy_var);
     	return error_message;
     }else{
 	pio_locus_t * locus;
        locus = bim_get_locus(&plink_file->bim_file, 0);
        unsigned char chromosome = locus->chromosome;
        int snp_start = 0;
        int snp_end = 0;
        int last_batch_size = batch_size;
        while(snp_start < plink_file->bed_file.header.num_loci){
        
            int snp = snp_start;
            while(locus->chromosome == chromosome && snp < plink_file->bed_file.header.num_loci) locus = bim_get_locus(&plink_file->bim_file, snp++);
            snp_end = snp;
            int snp_batch_size = snp_end - snp_start;
            N_SNPS_COMPUTED = 0;
            N_SNPS_LEFT = snp_batch_size;
            int current_batch_size = last_batch_size;
            if(snp_batch_size < n_threads*last_batch_size){
            	current_batch_size = floor(snp_batch_size/n_threads);
            	for(int i = 0; i < n_threads; i++){
            		thread_data[i]->change_max_batch_size(current_batch_size);
            		if(thread_data[i]->get_status()){
            			delete [] empirical_pedigree;
 				delete [] hetero_row_count;
 				delete [] hetero_col_count;
            			const char * error_status_message = thread_data[i]->check_error_status();
            			for(int j = 0; j < n_threads; j++){
            				delete thread_data[j];
            			}
            			delete [] thread_data;
            			error_message = "Error: "  + string(error_status_message) + " CPU Thread ID: " + to_string(i);
            			return error_message;
            		}
            	}
            }
            
            for(int index =0; index < n_threads;index++){
            	thread_data[index]->set_arrays_to_zero();
            }
            
           
            cout << "Starting calculation for chromosome: " << (unsigned)chromosome << " containing " <<  snp_batch_size << " loci\n";
            cout << "Total loci left to calculate: " <<   plink_file->bed_file.header.num_loci - snp_start << endl;
            
    	    vector<thread> cpu_threads;
	    int dummy_var = 0;
 	    
 	    fill_buffer_thread_king(0, dummy_var,   batch_size , plink_file, index_map,  n_subjects, snp_batch_size);
 	    start = std::chrono::high_resolution_clock::now();
 	    for(int thread_index = 0; thread_index < n_threads; thread_index++){
 		cpu_threads.push_back(thread(kingpedcalcidlistwrapper_,(unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count,\
     		thread_data[thread_index]->hetero_row_count,  &n_subjects, &array_size, &batch_size));
     	    }            
            
            memset(empirical_pedigree, 0, sizeof(float)*array_size);
    	    memset(hetero_row_count, 0, sizeof(int)*array_size);
            memset(hetero_col_count, 0, sizeof(int)*array_size);


     	    for(int thread_index = 0; thread_index < n_threads; thread_index++){
     		cpu_threads[thread_index].join();
     		for(int index = 0; index < array_size; index++){
     			empirical_pedigree[index] += thread_data[thread_index]->squared_difference_sum[index];
     			hetero_row_count[index] += thread_data[thread_index]->hetero_row_count[index];
     			hetero_col_count[index] += thread_data[thread_index]->hetero_col_count[index];
     		}
     		
     	    }
     	    
     	    for(int index = 0; index < array_size; index++){
     	    	empirical_pedigree[index] = 1.0 - empirical_pedigree[index]/(2.f*min(hetero_row_count[index], hetero_col_count[index]));
     		

    	    } 

   	    string str_output_filename = string(output_filename) + ".chr" + to_string(chromosome) + ".csv";
    	    write_epedigree_to_file(str_output_filename.c_str(), empirical_pedigree, plink_file,  n_subjects,  0);  
     	 
     	    cout << "King Empirical pedigree creation is complete for chromosome "  << endl;     	    
     	      
            if(snp_end != plink_file->bed_file.header.num_loci){
               chromosome = locus->chromosome;
            }
            snp_start = snp_end;
            
        }
        for(int index =0 ;index < n_threads; index++){
        	delete thread_data[index];
        }
        delete [] thread_data;
        delete [] empirical_pedigree;   
	delete [] hetero_row_count;
	delete [] hetero_col_count;
	int dummy_var = -1;
	fill_buffer_thread_king(0, dummy_var);	
	return error_message;
    	
     }
}




    		
static string calculate_robust_king_empirical_pedigree(pio_file_t * plink_file, const char * output_filename, const int per_chromosome, int batch_size, const int n_threads){
	string error_message;
	int n_subjects = plink_file->bed_file.header.num_samples;
	int array_size = n_subjects*(n_subjects + 1)/2;
	const int n_snps = plink_file->bed_file.header.num_loci;
	float * empirical_pedigree = new(nothrow) float[array_size];
	int * hetero_row_count = new(nothrow) int[array_size];
	int * hetero_col_count = new(nothrow) int[array_size];
    	if(!empirical_pedigree) return string("Failed to allocate memory for storing empirical pedigree");
    	memset(empirical_pedigree, 0, sizeof(float)*array_size);
    	
 //   	float * variance_sum = new(nothrow) float[array_size];
 //    	if(!variance_sum) return string("Failed to allocate memory for storing variance sums");
 //   	memset(variance_sum, 0, sizeof(float)*array_size);
    	King_Empirical_Pedigree_Thread ** thread_data = new(nothrow) King_Empirical_Pedigree_Thread*[n_threads]; 
     	if(!thread_data){
 		delete [] hetero_row_count;
 		delete [] hetero_col_count;
    		delete [] empirical_pedigree;
    		return string("Failed to allocate memory for storing CPU thread data");
   	 }   	  	
   
 	for(int index = 0; index < n_threads; index++){
 		thread_data[index] = new(nothrow) King_Empirical_Pedigree_Thread(plink_file, 0, n_subjects, n_subjects, batch_size);
 		if(!thread_data[index]){
 			delete [] empirical_pedigree;
 			delete [] hetero_row_count;
 			delete [] hetero_col_count;
 			for(int i = 0; i < index; i++) delete thread_data[index];
 			delete [] thread_data;
 			error_message = "Failed to allocate memory for CPU thread ID: " +to_string(index);
 			return error_message;
 		}
 		const char * error_status = 0;
		if( error_status = thread_data[index]->check_error_status()){
    			delete [] empirical_pedigree;
 			delete [] hetero_row_count;
 			delete [] hetero_col_count;
    			for(int i =0; i < index; i++) delete thread_data[i];
    			delete [] thread_data;
    			error_message = "Error: "  + string(error_status) + " CPU Thread ID: " + to_string(index);
    			return error_message;
		}
	}

	
    	
     if(!per_chromosome){
     
     	N_SNPS_COMPUTED = 0;
    	N_SNPS_LEFT = n_snps;
    	vector<thread> cpu_threads;
	int dummy_var = 0;
 	
 	fill_buffer_thread_king(0, dummy_var,   batch_size , plink_file, 0,  n_subjects);
 	start = std::chrono::high_resolution_clock::now();
 	for(int thread_index = 0; thread_index < n_threads; thread_index++){
 		cpu_threads.push_back(thread(kingpedcalcwrapper_,(unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count,\
     		thread_data[thread_index]->hetero_row_count, &n_subjects, &array_size, &batch_size));
     	}

 	/*fill_buffer(0, dummy_var,  n_subjects, batch_size, plink_file);  
 	#pragma omp parallel
 	{
 		const int thread_index =omp_get_thread_num(); 
 	     	kingpedcalcwrapper_((unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count, thread_data[thread_index]->hetero_row_count, &n_subjects,\
     		&array_size, &batch_size);
     		  
     	} */	
    	memset(hetero_row_count, 0, sizeof(int)*array_size);
        memset(hetero_col_count, 0, sizeof(int)*array_size);
     	for(int thread_index = 0; thread_index < n_threads; thread_index++){
		cpu_threads[thread_index].join();
     		for(int index = 0; index < array_size; index++){
     			empirical_pedigree[index] += thread_data[thread_index]->squared_difference_sum[index];
     			//variance_sum[index] += thread_data[thread_index]->variance_sum[index];
     			hetero_row_count[index] += thread_data[thread_index]->hetero_row_count[index];
     			hetero_col_count[index] += thread_data[thread_index]->hetero_col_count[index];
     		}
     		delete thread_data[thread_index];
     	}
     	delete [] thread_data;
     	for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] = 1.0 - empirical_pedigree[index]/(2.f*min(hetero_row_count[index], hetero_col_count[index]));
    	 }

    	write_epedigree_to_file(output_filename, empirical_pedigree, plink_file,  n_subjects,  0);  
     	
     	cout << "King Empirical pedigree creation is complete\n    ";// << setw(max_string_length + 1) << "\n";
     	delete [] empirical_pedigree;
 	delete [] hetero_row_count;
 	delete [] hetero_col_count;
        dummy_var = -1;
    	fill_buffer_thread_king(0, dummy_var);
     }else{
 	pio_locus_t * locus;
        locus = bim_get_locus(&plink_file->bim_file, 0);
        unsigned char chromosome = locus->chromosome;
        int snp_start = 0;
        int snp_end = 0;
        int last_batch_size = batch_size;
        while(snp_start < plink_file->bed_file.header.num_loci){
        
            int snp = snp_start;
            while(locus->chromosome == chromosome && snp < plink_file->bed_file.header.num_loci) locus = bim_get_locus(&plink_file->bim_file, snp++);
            snp_end = snp;
            int snp_batch_size = snp_end - snp_start;
            N_SNPS_COMPUTED = 0;
            N_SNPS_LEFT = snp_batch_size;
            int current_batch_size = last_batch_size;
            if(snp_batch_size < n_threads*last_batch_size){
            	current_batch_size = floor(snp_batch_size/n_threads);
            	for(int i = 0; i < n_threads; i++){
            		thread_data[i]->change_max_batch_size(current_batch_size);
            		if(thread_data[i]->get_status()){
            			delete [] empirical_pedigree;
 				delete [] hetero_row_count;
 				delete [] hetero_col_count;
            			const char * error_status_message = thread_data[i]->check_error_status();
            			for(int j = 0; j < n_threads; j++){
            				delete thread_data[j];
            			}
            			delete [] thread_data;
            			error_message = "Error: "  + string(error_status_message) + " CPU Thread ID: " + to_string(i);
            			return error_message;
            		}
            	}
            }
            
            for(int index =0; index < n_threads;index++){
            	thread_data[index]->set_arrays_to_zero();
            }
           
           
          

            cout << "Starting calculation for chromosome: " << (unsigned)chromosome << " containing " <<  snp_batch_size << " loci\n";
            cout << "Total loci left to calculate: " <<   plink_file->bed_file.header.num_loci - snp_start << endl;
    	vector<thread> cpu_threads;
	int dummy_var = 0;
 	
 	fill_buffer_thread_king(0, dummy_var,   batch_size , plink_file, 0,  n_subjects, snp_batch_size);
 	start = std::chrono::high_resolution_clock::now();
 	for(int thread_index = 0; thread_index < n_threads; thread_index++){
 		cpu_threads.push_back(thread(kingpedcalcwrapper_,(unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count,\
     		thread_data[thread_index]->hetero_row_count, &n_subjects, &array_size, &batch_size));
     	}  
 	    
	/*
 	#pragma omp parallel
 	{
 		const int thread_index =omp_get_thread_num(); 
 	     	kingpedcalcwrapper_((unsigned char *)thread_data[thread_index]->buffer,\
     		thread_data[thread_index]->squared_difference_sum,thread_data[thread_index]->hetero_col_count, thread_data[thread_index]->hetero_row_count, &n_subjects,\
     		&array_size, &batch_size);  
     	}  */	               

     	    memset(empirical_pedigree, 0, sizeof(float)*array_size);  
  	    memset(hetero_row_count, 0, sizeof(int)*array_size);
       	    memset(hetero_col_count, 0, sizeof(int)*array_size);    		 	   
     	    for(int thread_index = 0; thread_index < n_threads; thread_index++){
     	
     		for(int index = 0; index < array_size; index++){
     			empirical_pedigree[index] = thread_data[thread_index]->squared_difference_sum[index];
     			hetero_row_count[index] = thread_data[thread_index]->hetero_row_count[index];
     			hetero_col_count[index] = thread_data[thread_index]->hetero_col_count[index];
     		}
     		
     	    }
     	    
     	    for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] = 1.0 - empirical_pedigree[index]/(2.f*min(hetero_row_count[index], hetero_col_count[index]));

    	    } 

   	    string str_output_filename = string(output_filename) + ".chr" + to_string(chromosome) + ".csv";
    	    write_epedigree_to_file(str_output_filename.c_str(), empirical_pedigree, plink_file,  n_subjects,  0);  
     	   // cout << std::setw(max_string_length + 3) << endl;
     	    cout << "King Empirical pedigree creation is complete for chromosome " << (unsigned)chromosome <<  "    " << endl;     	    
     	      
            if(snp_end != plink_file->bed_file.header.num_loci){
               chromosome = locus->chromosome;
            }
            snp_start = snp_end;
            
        }
        for(int index =0 ;index < n_threads; index++){
        	delete thread_data[index];
        }
        delete [] thread_data;
        delete [] empirical_pedigree;    
	delete [] hetero_row_count;
	delete [] hetero_col_count;	
	
    	int dummy_var = -1;
    	fill_buffer_thread_king(0, dummy_var);
     } 
     return error_message;
}
static volatile size_t TOTAL_SNPS_COMPUTED = 0;
static void load_snp_buffer(snp_t * buffer, pio_file_t * plink_file, size_t & current_batch_size, size_t max_batch_size){
	const std::lock_guard<std::mutex> lock(mtx);
	static pio_status_t status = PIO_OK;
	if(current_batch_size != 0){
		TOTAL_SNPS_COMPUTED += current_batch_size;
	}
	cout << "Progress: " << 100.f*TOTAL_SNPS_COMPUTED/plink_file->bed_file.header.num_loci << "% complete\r";
	cout.flush();	
	current_batch_size = 0;
	if(status == PIO_END) return;
	
	while(status == PIO_OK && current_batch_size < max_batch_size){
		status = pio_next_row(plink_file, buffer + current_batch_size*plink_file->bed_file.header.num_samples);	
		current_batch_size++;
	}
	
}	
static string calculate_correlation_empirical_pedigree(pio_file_t * plink_file,const char * frequencies_filename, const char * output_filename,  float alpha, const int per_chromosome,  int batch_size , \
							 int n_threads, const bool normalize, const bool use_method_one){
    string errmsg;
    int n_subjects = plink_file->bed_file.header.num_samples;
    int array_size =    n_subjects*(n_subjects + 1)/2; 
    const int n_snps= plink_file->bed_file.header.num_loci;
   
    float * variance_sum = 0;
    int * snp_count = 0;
    float * empirical_pedigree = new(nothrow) float[array_size];
    if(!empirical_pedigree) return string("Failed to allocate memory for storing empirical pedigree");
    memset(empirical_pedigree, 0, sizeof(float)*array_size);
    if(use_method_one){
    	snp_count = new(nothrow) int[array_size];
    	if(!snp_count) return string("Failed to allocate memory for storing empirical pedigree");
    }else{
    	variance_sum = new(nothrow) float[array_size];
    	if(!variance_sum) return string("Failed to allocate memory for storing empirical pedigree");    
    }	


    Empirical_Pedigree_Thread ** thread_data = new(nothrow) Empirical_Pedigree_Thread*[n_threads];
    if(!thread_data){
    	if(use_method_one) 
    		delete [] snp_count;
    	else
    		delete [] variance_sum;
    	delete [] empirical_pedigree;
    	return string("Failed to allocate memory for storing CPU thread data");
    }
    	
    for(int index = 0; index < n_threads; index++){

    	thread_data[index] = new(nothrow) Empirical_Pedigree_Thread(plink_file, 0, n_subjects,\
    					 n_subjects, alpha, batch_size, use_method_one);
   	
    	if(!thread_data[index] ){
    		delete [] empirical_pedigree;
    		if(use_method_one) 
    			delete [] snp_count;
    		else
    			delete [] variance_sum;
    		for(int i =0; i < index; i++) delete thread_data[i];
    		delete [] thread_data;
    		errmsg = "Failed to allocate memory for CPU thread ID: " + to_string(index);
    		return errmsg;
    	}
	const char * error_status = 0;
	if( error_status = thread_data[index]->check_error_status()){
    		delete [] empirical_pedigree;
    		if(use_method_one) 
    			delete [] snp_count;
    		else
    			delete [] variance_sum;
    		for(int i =0; i < index; i++) delete thread_data[i];
    		delete [] thread_data;
    		errmsg = "Error: "  + string(error_status) + " CPU Thread ID: " + to_string(index);
    		return errmsg;
	}
    	
     }
     	
     if(!per_chromosome){
     
     	N_SNPS_COMPUTED = 0;
    	N_SNPS_LEFT = n_snps;
    	int dummy_var = 0;
    	
	fill_buffer_thread(0, 0, dummy_var,   batch_size , frequencies_filename, plink_file, 0,  n_subjects);
    	vector<thread> cpu_threads;
    	if(use_method_one){
    		start = std::chrono::high_resolution_clock::now();
     		for(int thread_index = 0; thread_index < n_threads; thread_index++){
     			cpu_threads.push_back(thread(corrpedcalcwrapperone_, thread_data[thread_index]->snp_count, \
     			thread_data[thread_index]->empirical_pedigree,thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	&n_subjects, &alpha, &array_size, &batch_size));
     		}
     	}else{
     		start = std::chrono::high_resolution_clock::now();
     		for(int thread_index = 0; thread_index < n_threads; thread_index++){
     			cpu_threads.push_back(thread(corrpedcalcwrappertwo_, thread_data[thread_index]->variance_sum, \
     			thread_data[thread_index]->empirical_pedigree, thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	&n_subjects, &array_size, &batch_size));
     		}     		
	}
    				
	if(use_method_one){	
		memset(snp_count, 0, sizeof(int)*array_size);
	}else{	
		memset(variance_sum, 0, sizeof(float)*array_size);
	}
	
     	for(int thread_index = 0; thread_index < n_threads; thread_index++){
     		cpu_threads[thread_index].join();
     		if(use_method_one){
     			for(int index = 0; index < array_size; index++){
     				empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     				snp_count[index] += thread_data[thread_index]->snp_count[index];
     			}
     		}else{
     			for(int index = 0; index < array_size; index++){
     				empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     				variance_sum[index] += thread_data[thread_index]->variance_sum[index];
     			}
     		}     			
     		delete thread_data[thread_index];
     	}
     	
     	delete [] thread_data;
     
     	for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] /= ((use_method_one) ? snp_count[index] : variance_sum[index]);//snp_count[index];
    	 }
    	 if(normalize){
     		float norms[n_subjects];
     		for(int index = 0; index < n_subjects; index++){
     			norms[index] = sqrt(empirical_pedigree[create_index(index,index,n_subjects)]);
     		//empirical_pedigree[create_index(index,index, n_subjects)] = 1.f;
     		}
     		for(int col = 0; col < n_subjects; col++){
     			for(int row = col ; row < n_subjects; row++){
     				empirical_pedigree[create_index(row,col, n_subjects)] /= (norms[row]*norms[col]);
     			}
   		}
   	}
    //	write_epedigree_to_file_2("snp_counts_2.csv", variance_sum_or_snp_count, plink_file,  n_subjects,  0);  
     	write_epedigree_to_file(output_filename, empirical_pedigree, plink_file,  n_subjects,  0); 
     	cout << "\n"; 
     	cout << "Empirical pedigree creation is complete                   \n";// << setw(max_string_length + 1) << "\n";
     	delete [] empirical_pedigree;
    	if(use_method_one) 
    		delete [] snp_count;
    	else
    		delete [] variance_sum;
   	dummy_var = -1;
	fill_buffer_thread(0, 0, dummy_var);
    }else{
 	pio_locus_t * locus;
        locus = bim_get_locus(&plink_file->bim_file, 0);
        unsigned char chromosome = locus->chromosome;
        int snp_start = 0;
        int snp_end = 0;
        int last_batch_size = batch_size;
        while(snp_start < plink_file->bed_file.header.num_loci){
        
            int snp = snp_start;
            while(locus->chromosome == chromosome && snp < plink_file->bed_file.header.num_loci) locus = bim_get_locus(&plink_file->bim_file, snp++);
            snp_end = snp;
            int snp_batch_size = snp_end - snp_start;
            N_SNPS_COMPUTED = 0;
            N_SNPS_LEFT = snp_batch_size;
    	//    Empirical_Pedigree_Thread * dummy_ptr = 0;
    	//    int dummy_var = 0;  
     	//    fill_buffer_2<Empirical_Pedigree_Thread>(dummy_ptr, dummy_var, true);            
            int current_batch_size = last_batch_size;
            if(snp_batch_size < n_threads*last_batch_size){
            	current_batch_size = floor(snp_batch_size/n_threads);
            	for(int i = 0; i < n_threads; i++){
            		thread_data[i]->change_max_batch_size(current_batch_size);
            		if(thread_data[i]->get_status()){
            			delete [] empirical_pedigree;
    				if(use_method_one) 
    					delete [] snp_count;
    				else
    					delete [] variance_sum;
            			const char * error_status_message = thread_data[i]->check_error_status();
            			for(int j = 0; j < n_threads; j++){
            				delete thread_data[j];
            			}
            			delete [] thread_data;
            			errmsg = "Error: "  + string(error_status_message) + " CPU Thread ID: " + to_string(i);
            			return errmsg;
            		}
            	}
            }
            
            for(int index =0; index < n_threads;index++){
            	thread_data[index]->set_arrays_to_zero();
            }
            memset(empirical_pedigree, 0, sizeof(float)*array_size);
	    if(use_method_one){	
		memset(snp_count, 0, sizeof(int)*array_size);
	    }else{	
		memset(variance_sum, 0, sizeof(float)*array_size);
	    }
            cout << "Starting calculation for chromosome: " << (unsigned)chromosome << " containing " <<  snp_batch_size << " loci\n";
            cout << "Total loci left to calculate: " <<   plink_file->bed_file.header.num_loci - snp_start << endl;
    	    int dummy_var = 0;
    	    
            fill_buffer_thread(0, 0, dummy_var,   batch_size , frequencies_filename, plink_file, 0,  n_subjects, snp_batch_size);
    	    vector<thread> cpu_threads;
    	    if(use_method_one){
    	    	start = std::chrono::high_resolution_clock::now();
     		for(int thread_index = 0; thread_index < n_threads; thread_index++){
     			cpu_threads.push_back(thread(corrpedcalcwrapperone_, thread_data[thread_index]->snp_count, \
     			thread_data[thread_index]->empirical_pedigree, thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	&n_subjects, &alpha, &array_size,  &batch_size));
     		}
     	    }else{
     	    	start = std::chrono::high_resolution_clock::now();
     		for(int thread_index = 0; thread_index < n_threads; thread_index++){
     			cpu_threads.push_back(thread(corrpedcalcwrappertwo_, thread_data[thread_index]->variance_sum, \
     			thread_data[thread_index]->empirical_pedigree,thread_data[thread_index]->frequencies, thread_data[thread_index]->buffer,\
     		 	&n_subjects, &array_size, &batch_size));
     		}     		
	    }

  	    
     	    for(int thread_index = 0; thread_index < n_threads; thread_index++){
     	    	cpu_threads[thread_index].join();
     		if(use_method_one){
     			for(int index = 0; index < array_size; index++){
     				empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     				snp_count[index] += thread_data[thread_index]->snp_count[index];
     			}
     		}else{
     			for(int index = 0; index < array_size; index++){
     				empirical_pedigree[index] += thread_data[thread_index]->empirical_pedigree[index];
     				variance_sum[index] += thread_data[thread_index]->variance_sum[index];
     			}
     		}
     	    }     	    

     	    
     	    for(int index = 0; index < array_size; index++){
     		empirical_pedigree[index] /= ((use_method_one) ? snp_count[index] : variance_sum[index]);	
    	    } 
    	    if(normalize){
     	    	float norms[n_subjects];
     	    	for(int index = 0; index < n_subjects; index++){
     			norms[index] = sqrt(empirical_pedigree[create_index(index,index,n_subjects)]);
     		//empirical_pedigree[create_index(index,index, n_subjects)] = 1.f;
     	    	}
     	    	for(int col = 0; col < n_subjects; col++){
     			for(int row = col ; row < n_subjects; row++){
     				empirical_pedigree[create_index(row,col, n_subjects)] /= (norms[row]*norms[col]);
     			}
   	    	} 
   	    }
   	    string str_output_filename = string(output_filename) + ".chr" + to_string(chromosome) + ".csv";
    	    write_epedigree_to_file(str_output_filename.c_str(), empirical_pedigree, plink_file,  n_subjects,  0);  
            cout << "\n";
     	    cout << "Empirical pedigree creation is complete for chromosome " << (unsigned)chromosome << "             " << endl;     	    
     	      
            if(snp_end != plink_file->bed_file.header.num_loci){
               chromosome = locus->chromosome;
            }
            snp_start = snp_end;
            
        }
        for(int index =0 ;index < n_threads; index++){
        	delete thread_data[index];
        }
        delete [] thread_data;
        delete [] empirical_pedigree;
        if(use_method_one)
        	delete [] snp_count;
        else
        	delete [] variance_sum;
        
       int dummy_var = -1;
       fill_buffer_thread(0, 0, dummy_var); 
    }

    return errmsg;
}



static void print_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help pedifromsnps");
    
}
static vector<string> read_id_list(const char * id_list_filename){
	vector<string> output;
	ifstream input_stream(id_list_filename);
	if(!input_stream.is_open()) return output;
	string id;
	while(input_stream >> id){
		output.push_back(id);
	}
	return output;
}

extern "C" int pedfromsnpsCmd(ClientData clientData, Tcl_Interp *interp,
                              int argc,const char *argv[]){
    
    int use_king = 0;
    int per_chromosome = 0;
    int batch_size = 500;
    const char * plink_filename = 0;
    const char * output_filename = 0;
    const char * frequency_filename = 0;
    double alpha = -1.0;
    const char * id_list_filename = 0;
    int n_threads = std::thread::hardware_concurrency();
    bool normalize = false;
    bool use_method_one = true;
    
    for(int arg = 1; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--i", case_ins) || \
            !StringCmp(argv[arg], "--input", case_ins) \
            || !StringCmp(argv[arg], "-i", case_ins)\
            || !StringCmp(argv[arg], "-input", case_ins)) && arg + 1 < argc){
            plink_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--o", case_ins) || \
                  !StringCmp(argv[arg], "--output", case_ins) \
                  || !StringCmp(argv[arg], "-o", case_ins)\
                  || !StringCmp(argv[arg], "-output", case_ins)) && arg + 1 < argc){
            output_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "-corr", case_ins) || \
                  !StringCmp(argv[arg], "--corr", case_ins)) && arg + 1 < argc){
            alpha = atof(argv[++arg]);
           
        }else if((!StringCmp(argv[arg], "-freq", case_ins) || \
                  !StringCmp(argv[arg], "--freq", case_ins)) && arg + 1 < argc){
            frequency_filename = argv[++arg];
           
        }else if(!StringCmp(argv[arg], "-king", case_ins) || \
                  !StringCmp(argv[arg], "--king", case_ins) ){
    
            use_king = 1;
        }else if(!StringCmp(argv[arg], "-method_two", case_ins) || \
                  !StringCmp(argv[arg], "--method_two", case_ins) ){
    
            use_method_one = false;
        }else if(!StringCmp(argv[arg], "-normalize", case_ins) || \
                  !StringCmp(argv[arg], "--normalize", case_ins) ){
    
            normalize = true;
        }else if(!StringCmp(argv[arg], "-per-chromo", case_ins) || \
                 !StringCmp(argv[arg], "--per-chromo", case_ins)){
            per_chromosome = 1;
        }else if(!StringCmp(argv[arg], "-help", case_ins) || \
                 !StringCmp(argv[arg], "--help", case_ins) || \
                 !StringCmp(argv[arg], "help", case_ins) ){
            print_help(interp);
            return TCL_OK;
        }else if((!StringCmp(argv[arg], "-id_list", case_ins) || \
                  !StringCmp(argv[arg], "--id_list", case_ins)) && arg + 1 < argc){
            id_list_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "-batch_size", case_ins) || \
                  !StringCmp(argv[arg], "--batch_size", case_ins)) && arg + 1 < argc){
            batch_size = atoi(argv[++arg]);
            if(batch_size == 0){
            	RESULT_LIT("--batch_size must be greater than zero");
            	return TCL_ERROR;
            }
        }else if((!StringCmp(argv[arg], "-n_threads", case_ins) || \
                  !StringCmp(argv[arg], "--n_threads", case_ins)) && arg + 1 < argc){
            cout << "Optimal number of threads for this computer is " << n_threads << " (Default)\n";
            n_threads = atoi(argv[++arg]);
            if(n_threads <= 0 || n_threads >= 10 + std::thread::hardware_concurrency()){
            	cout << "--n_threads value must be greater than 0 and less than " << std::thread::hardware_concurrency() + 10 << std::endl;
            	return TCL_ERROR;
            }
            cout << "Using " << n_threads << " CPU threads instead of default value\n"; 
        }else{
            string error_message = "Invalid argument was entered: " + string(argv[arg]);
            cout << error_message << endl;
           // RESULT_BUF(error_message.c_str());
            return TCL_ERROR;
        }
    }
    
    if(!plink_filename){
         cout << "No plink file set specified with --i"<< endl;
        return TCL_ERROR;
    }
    
    if(!output_filename){
         cout << "No output filename specified with --o"<< endl;
        return TCL_ERROR;
    } 
    bool use_one_loci_per_row_method = true;
    pio_file_t * plink_file = new pio_file_t;
    pio_status_t status;
    status = pio_open(plink_file, plink_filename);
    std::cout << "Number of loci included in GRM: " << plink_file->bed_file.header.num_loci << std::endl;
    std::cout << "Total number of samples in plink file set: " << plink_file->bed_file.header.num_samples << std::endl;
    if(status != PIO_OK){
        if(status == P_FAM_IO_ERROR){
            cout << "Error in loading .fam file"<< endl;
            return TCL_ERROR;
        }else if (status == P_BIM_IO_ERROR){
             cout << "Error in loading .bim file"<< endl;
            return TCL_ERROR;
        }else if (status == P_BED_IO_ERROR){
             cout << "Error in loading .bed file"<< endl;
            return TCL_ERROR;
        }else{
             cout << "Error loading plink file"<< endl;;
            return TCL_ERROR;
        }
    }
    int version = plink_file->bed_file.header.version;
    if (plink_file->bed_file.header.snp_order == BED_UNKNOWN_ORDER){
        pio_close(plink_file);
        cout << "Error in the .bed snp order. Retry creation of file using a different version of plink" << endl;
        
        return TCL_ERROR;
        
    }
    
    if (plink_file->bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
    	use_one_loci_per_row_method = false;
    	cout << "Using one loci per row method\n";
    } else{
    	cout << "Using one sample per row method\n";  
    }  
       
    /*else if (plink_file->bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
        pio_close(plink_file);
        printf("In order to read efficiently the transpose of specified plink file must be performed\n");
        string transpose_filename = string(plink_filename) + string(".trans");
        string message = "Filename of transposed plink file is " + transpose_filename + "\n";
        printf(message.c_str());
        status = pio_transpose(plink_filename, transpose_filename.c_str());
        if(status != PIO_OK){
            RESULT_LIT("Error in creating transpose");
            return TCL_ERROR;
        }
        
        status = pio_open(plink_file, transpose_filename.c_str());
        
        if(status != PIO_OK){
            printf("Error in opening transposed plink file\n");
            return TCL_ERROR;
        }
    }*/
    vector<string> id_vector; 
    int * index_map = 0;
    int index_map_size = 0;
    if(id_list_filename){
	  id_vector = read_id_list(id_list_filename);
	  if(id_vector.size() == 0 ){
		cout << "No IDs read from id list option" << endl;
		return TCL_ERROR;
	  }
	  vector<int> id_index_map;
	  pio_sample_t * sample; 
	  for(int index = 0; index < pio_num_samples(plink_file); index++){
		sample = pio_get_sample(plink_file, index);
		string str_id = string(sample->iid);
		vector<string>::iterator find_iter = find(id_vector.begin(), id_vector.end(), str_id);
		if(find_iter != id_vector.end()){
			id_index_map.push_back(index);
			id_vector.erase(find_iter);
		}
		
	}
	if(id_index_map.size() == 0){
		cout  << "None of the IDs listed in ID list file that was specified were found in the plink data set" << endl;
		return TCL_ERROR;
	}
	index_map = new int[id_index_map.size()];
	for(int index = 0 ; index < id_index_map.size(); index++){
		index_map[index] = id_index_map[index];
	}
	index_map_size = id_index_map.size();
	
   }    	         
    if(index_map_size != 0){
    	std::cout << "Building GRM using sample subset read from ID list file " << id_list_filename << std::endl;
    	std::cout << "Number of samples used in sample subset: " << index_map_size << std::endl;
    }
    const char * errmsg = 0;
    if(use_king){
    	cout << "Creating GRM using Robust King Method\n";
       // create_empirical_pedigree_allele_sharing( plink_file,  output_filename, per_chromosome, index_map, index_map_size);
      /* if(use_king_homogenous){
       		string error_message = calculate_king_homogenous_empirical_pedigree(plink_file, output_filename, per_chromosome, batch_size,  n_threads);
                if(error_message.length() != 0){
            	    cout << errmsg << endl;
            	   return TCL_ERROR;
                }
        }  */     		
	if(index_map_size == 0){
		string error_message = calculate_robust_king_empirical_pedigree(plink_file, output_filename, per_chromosome,  batch_size,  n_threads);
                if(error_message.length() != 0){
            	    cout << errmsg << endl;
            	   return TCL_ERROR;
                }		
	}else{
		
		
		string error_message = calculate_robust_king_empirical_pedigree_id_list(plink_file,  output_filename,  index_map,  per_chromosome,  index_map_size,  batch_size,  n_threads);
                if(error_message.length() != 0){
            	    cout << errmsg << endl;
            	   return TCL_ERROR;
                }
	}
    }else{
    	if(!frequency_filename){
    		std::cout << "Correlation GRM requires a frequency file specified with --freq option and computed using calculate_plink_freq command\n";
    		pio_close(plink_file);
    		return TCL_ERROR;
    	}
    	std::ifstream test_stream(frequency_filename);
    	if(test_stream.is_open() == false){
    		std::cout << "Failed to open " << frequency_filename << " frequency file\n";
    		pio_close(plink_file);
    		return TCL_ERROR;
    	}
    	std::string first_line;
    	getline(test_stream,first_line);
    	const int n_freq = stoi(first_line);
    	if(n_freq != plink_file->bed_file.header.num_loci){
    		std::cout << "Number of frequencies from frequency file and number of loci from plink file do not match\n";
    		std::cout << "Number of frequencies read from " << frequency_filename << ": " << n_freq << std::endl;
    		std::cout << "Number of loci read from " << plink_filename  << ": " << plink_file->bed_file.header.num_loci << std::endl; 
    		test_stream.close();
    		pio_close(plink_file);
    		return TCL_ERROR;
    	}
    	test_stream.close();
        if(use_method_one){
        	cout << "Creating GRM using Correlation Method One (Default Method)\n";
        	
        }else{
        	cout << "Creating GRM using Correlation Method Two\n"; 
        }
 	if(normalize) std::cout << "Final values will be normalized so diagonal elements are all one and off diagonal elements are bounded by one and negative one\n";
	if(index_map_size == 0){
	   

	    string error_message  = calculate_correlation_empirical_pedigree(plink_file, frequency_filename, output_filename, alpha, per_chromosome,\
	    				 batch_size, n_threads, normalize, use_method_one);
	    				
            if(error_message.length() != 0 ){
            	cout << error_message << endl;
            	return TCL_ERROR;
            }	    
	}else{

	    string error_message  = calculate_correlation_empirical_pedigree_with_id_list(plink_file, frequency_filename, output_filename,  index_map, alpha, per_chromosome,\
	    				 index_map_size, batch_size, n_threads, normalize, use_method_one);
            if(error_message.length() != 0 ){
            	cout << error_message << endl;
            	return TCL_ERROR;
            }	
	}
      
        
    }
    if(index_map) delete [] index_map;
    pio_close(plink_file);
    delete plink_file;
  //  RESULT_LIT("Empirical pedigree creation is complete");
    if(errmsg){
    	cout << errmsg << endl;
    	return TCL_ERROR;
    }else{
   	 return TCL_OK;
   }
}



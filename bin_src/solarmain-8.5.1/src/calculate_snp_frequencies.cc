#include "plinkio.h"
#include "solar.h"
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <omp.h>
using namespace std;

static void print_calculate_loci_frequencies_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help plink_freq");
    
}

static void transpose_plink_file_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help transpose_plink_file");
    
}


static void calculate_frequency_chunk(const int thread_index,pio_file_t * plink_file, float * frequencies, const size_t start_index, const size_t end_index){
	snp_t * buffer = new snp_t[plink_file->bed_file.header.num_samples];
	size_t index = 0;
	while(index < start_index){
		pio_skip_row(plink_file);
		++index;
	}

	for(;index < end_index ; index++){
		pio_next_row(plink_file, buffer);
		unsigned sum = 0;
		unsigned n_missing = 0;
		for(size_t subject = 0; subject < plink_file->bed_file.header.num_samples; subject++){
			const snp_t value = buffer[subject];
			if(value != 3){
				sum += value;
			}else{
				++n_missing;
			}
		}
		if(n_missing != plink_file->bed_file.header.num_samples){
			frequencies[index] = (float)sum/(2.f*(plink_file->bed_file.header.num_samples - n_missing));
		}else{
			frequencies[index] = -1.f;
		}

	}
	delete [] buffer;
	cout << "thread: " << thread_index << " is done\n";
}
static void calculate_loci_frequencies(pio_file_t * plink_file, const char * plink_filename, const char * output_filename){
	const size_t n_samples = plink_file->bed_file.header.num_samples;
	const size_t n_loci = plink_file->bed_file.header.num_loci;
	std::cout << "n_loci: " << n_loci << " n_samples: " << n_samples << endl;
	float * frequencies = new float[n_loci];

	snp_t * buffer = new snp_t[n_samples];
	int largest_width = 0;
	for(size_t index = 0; index < n_loci; index++){
		pio_next_row(plink_file, buffer);
		unsigned sum = 0;
		unsigned n_missing = 0;
		for(size_t subject = 0; subject < n_samples; subject++){
			const snp_t value = buffer[subject];
			if(value != 3){
				sum += value;
			}else{
				++n_missing;
			}
		}
		if(n_missing != n_samples){
			frequencies[index] = (float)sum/(2.f*(n_samples - n_missing));
		}else{
			frequencies[index] = -1.f;
		}

	}
	delete [] buffer;				
	ofstream output_stream(output_filename);
	output_stream << n_loci << "\n";
	output_stream << frequencies[0];
	for(int loci = 1; loci < n_loci; loci++){
		output_stream << " " << frequencies[loci];
	}	
	output_stream.close();
	delete [] frequencies;
}	


extern "C" int calculate_loci_frequencies_command(ClientData clientData, Tcl_Interp *interp,
                              int argc,const char *argv[]){
        const char * plink_filename = 0;
        const char * output_filename = 0;                    
	for(int arg = 1; arg < argc; arg++){                            
        	if((!StringCmp(argv[arg], "--plink", case_ins) \
           	 || !StringCmp(argv[arg], "-plink", case_ins)) && arg + 1 < argc){
            		plink_filename = argv[++arg];
        	}else if((!StringCmp(argv[arg], "--o", case_ins) || \
                  	!StringCmp(argv[arg], "--output", case_ins) \
                 	|| !StringCmp(argv[arg], "-o", case_ins)\
                  	|| !StringCmp(argv[arg], "-output", case_ins)) && arg + 1 < argc){
           		output_filename = argv[++arg];
        	}else if(!StringCmp(argv[arg], "-help", case_ins) || \
                 	!StringCmp(argv[arg], "--help", case_ins) || \
                 	!StringCmp(argv[arg], "help", case_ins) ){
            		print_calculate_loci_frequencies_help(interp);
           		 return TCL_OK;
        	}else{
            		string error_message = "Invalid argument was entered: " + string(argv[arg]);
            		cout << error_message << endl;
           // RESULT_BUF(error_message.c_str());
            		return TCL_ERROR;
        	}
        }
        if(!plink_filename){
        	RESULT_LIT("plink filename was not specified");
        	return TCL_ERROR;
        }
         if(!output_filename){
        	RESULT_LIT("output filename was not specified");
        	return TCL_ERROR;
        }       
        pio_file_t * plink_file = new pio_file_t;
    	pio_status_t status;
    	status = pio_open(plink_file, plink_filename);

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
    	if(plink_file->bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
    		RESULT_LIT("plink file SNP order is one sample per row, use transpose command so SNP order is one locus per row");
    		return TCL_ERROR;
    	}
    	calculate_loci_frequencies(plink_file, plink_filename, output_filename);
    	pio_close(plink_file);
    	delete plink_file;
    	return TCL_OK;
}
                   
    	

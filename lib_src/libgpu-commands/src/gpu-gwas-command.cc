//
//  gpu_gwas_command.cc
//
//
//  Created by Brian Donohue on 7/16/18.
//
//
//#include "GPU_GWAS_Estimator.h"
#include <stdio.h>
#include "solar.h"
#include <string>
#include <vector>
#include <dlfcn.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <chrono>
#include <exception>
using namespace std;
#include "gpu-gwas-settings.h"
#include "gpu-gwas.h"
#include "gpu-exception.h"
#include "gpu-selection.h"



const char * gpu_gwas_screen(vector<int> gpu_ids, vector<string> trait_list, const char*  phenotype_filename, const char* plink_filename, const char * evd_data_filename, int defined_max_batch_size, int defined_thread_size, int n_streams, const bool verbose);
vector<string> read_list_of_traits(const char * list_filename){
    ifstream input(list_filename);
    vector<string> trait_list;
    if(input.is_open() != true){
        return trait_list;
    }
    string trait;
    while(!input.eof()) {
        input >> trait;
        trait_list.push_back(trait);
    }
    input.close();
    return trait_list;
}

extern "C" int  GPU_GWAS_Cmd(ClientData clientData, Tcl_Interp *interp,
                             int argc,const char *argv[]){
    const char * plink_filename = 0;
    size_t precision = 6;
    bool select_all_gpus = false;
    const char * gpu_list_string = 0;
    const char * list_filename = 0;
    const char * evd_data_filename = 0;
    const char * output_filename = 0;
    int blockSize = -1;
    int defined_batch_size = -1;
    int stream_count = 1;
    bool use_screen = false;
    bool run_calibration = false;
    double snp_percentage;
    bool verbose = false;
    for(int arg = 1 ; arg < argc; arg++){
        if(!StringCmp(argv[arg], "help", case_ins) || !StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins)){
            //            print_gpu_gwas_help();
            return TCL_OK;
        }else if ((!StringCmp(argv[arg], "--plink", case_ins) || !StringCmp(argv[arg], "-plink", case_ins)) && arg + 1 < argc){
            
            plink_filename = argv[++arg];
            
        }else if ((!StringCmp(argv[arg], "--o", case_ins) || !StringCmp(argv[arg], "-o", case_ins) || !StringCmp(argv[arg], "-out", case_ins) || \
        	!StringCmp(argv[arg], "--out", case_ins) || !StringCmp(argv[arg], "-output", case_ins) || !StringCmp(argv[arg], "--output", case_ins)) && arg + 1 < argc){
            
            output_filename = argv[++arg];
            
        }else if ((!StringCmp(argv[arg], "--precision", case_ins) || !StringCmp(argv[arg], "-precision", case_ins)) && arg + 1 < argc){
            
            precision = stoul(string(argv[++arg]));
            
        }else if ((!StringCmp(argv[arg], "--stream_count", case_ins) || !StringCmp(argv[arg], "-stream_count", case_ins)) && arg + 1 < argc){
            
            stream_count = stoul(string(argv[++arg]));
            if(stream_count <= 0){
            	RESULT_LIT("Defined stream count must be greater than zero");
            	return TCL_ERROR;
            } 
            if(stream_count > MAX_STREAM_COUNT){
            	string error_message = "Defined stream count is greater than stream count limit " + to_string(MAX_STREAM_COUNT);
            	RESULT_BUF(error_message.c_str());
            	return TCL_ERROR;
            }                
            
        }else if ((!StringCmp(argv[arg], "--gpus", case_ins) || !StringCmp(argv[arg], "-gpus", case_ins)) && arg + 1 < argc){
            
            gpu_list_string = argv[++arg];
            
        }else if (!StringCmp(argv[arg], "--all", case_ins) || !StringCmp(argv[arg], "-all", case_ins)){
            
            select_all_gpus = true;
            
        }else if ((!StringCmp(argv[arg], "--list", case_ins) || !StringCmp(argv[arg], "-list", case_ins)) && arg + 1 < argc){
            
            list_filename = argv[++arg];
            
        }else if ((!StringCmp(argv[arg], "--evd_data", case_ins) || !StringCmp(argv[arg], "-evd_data", case_ins)) && arg + 1 < argc){
            
            evd_data_filename = argv[++arg];
            
        }else if ((!StringCmp(argv[arg], "--size", case_ins) || !StringCmp(argv[arg], "-size", case_ins)) && arg + 1 < argc){

            blockSize  = stoul(string(argv[++arg]));
            if(blockSize != 32 && blockSize != 64 && blockSize != 128 && blockSize != 256 && blockSize != 512 && blockSize != 1024){
            	RESULT_LIT("Defined block size must be 32,64,128,256,512, or 1024");
            	return TCL_ERROR;
            }

        }else if ((!StringCmp(argv[arg], "--batch_size", case_ins) || !StringCmp(argv[arg], "-batch_size", case_ins)) && arg + 1 < argc){

            defined_batch_size  = stoul(string(argv[++arg]));
            if(defined_batch_size <= 0){
            	RESULT_LIT("Defined batch size must be greater than zero");
            	return TCL_ERROR;
            }
            if(defined_batch_size > MAX_BATCH_SIZE){
            	string error_message = "Defined batch size is greater than batch size limit " + to_string(MAX_BATCH_SIZE);
            	RESULT_BUF(error_message.c_str());
            	return TCL_ERROR;
            }

        }else if ((!StringCmp(argv[arg], "--calibrate", case_ins) || !StringCmp(argv[arg], "-calibrate", case_ins)) && arg + 1 < argc){
	    run_calibration = true;
            snp_percentage  = stod(string(argv[++arg]));

        }else if (!StringCmp(argv[arg], "--screen", case_ins) || !StringCmp(argv[arg], "-screen", case_ins)){
		use_screen = true;
        }else if (!StringCmp(argv[arg], "--verbose", case_ins) || !StringCmp(argv[arg], "-v", case_ins) || !StringCmp(argv[arg], "--v", case_ins) || !StringCmp(argv[arg], "--verbose", case_ins)){
		verbose = true;
	}else{
            RESULT_LIT("Invalid argument entered");
            return TCL_ERROR;
        }
    }
    if((snp_percentage > 1 || snp_percentage <= 0.0) && run_calibration){
	RESULT_LIT("Calibration percentage must be (0,1]");
	return TCL_ERROR;
     }
    if(!loadedPed() && evd_data_filename == 0){
        RESULT_LIT("No pedigree has been loaded");
        return TCL_ERROR;
    }
    
    if(Trait::Number_Of() == 0 && list_filename == 0){
        RESULT_LIT("No trait has been selected");
        return TCL_ERROR;
    }
    vector<string> trait_list;
    
    if(list_filename){
        trait_list = read_list_of_traits(list_filename);
        if(!trait_list.size()){
            RESULT_LIT("Trait list contained no traits or list file could not be opened");
            return TCL_ERROR;
        }
    }else{
        trait_list.push_back(string(Trait::Name(0)));
    }
    
    /*if(trait_list.size() == 1){
        string str_output_filename = "./" + trait_list[0] + "/gpu_gwas.out";
        output_filename = str_output_filename.c_str();
        if(access(trait_list[0].c_str(), F_OK)  == -1){
            mkdir(trait_list[0].c_str(), 0700);
        }
    }else{
        output_filename = "gpu_gwas.out";
    }*/
    if(!output_filename) output_filename = "gpu-gwas.out"; 
    
    const char * phenotype_filename = Phenotypes::filenames();
    
    
    if(!plink_filename){
        RESULT_LIT("No plink filename was specified");
        return TCL_ERROR;
    }
    // cudaProfilerStart();
    vector<int> gpu_id_list;
    try{
     	gpu_id_list = get_gpu_id_list(gpu_list_string, select_all_gpus);
     }catch(GPU_Exception & e){
     	RESULT_BUF(e.what());
     	return TCL_ERROR;
     }catch(...){
     	RESULT_LIT("Unknown error occurred creating GPU list");
     	return TCL_ERROR;
     }
    
    if(gpu_id_list.size() == 0){
        RESULT_LIT("No GPUs were selected");
        return TCL_ERROR;
    }
 
    if(!evd_data_filename){
    	try {
            load_phi2_matrix(interp);
   	}catch(exception & e){
       	     RESULT_BUF(e.what());
       	     return TCL_ERROR;
   	}catch(...){
             RESULT_LIT("Unkown error occurred when loading phi2 matrix");
             return TCL_ERROR;
    	}
    }
    if(use_screen){
    	const char * errmsg = 0;
    	try{
		errmsg = gpu_gwas_screen(gpu_id_list, trait_list, phenotype_filename, plink_filename,  evd_data_filename,  defined_batch_size,  blockSize,  stream_count,verbose);
	}catch(Solar_Trait_Reader_Exception & e){
		RESULT_BUF(e.what());
		return TCL_ERROR;
	}catch(GPU_Exception & e){
		RESULT_BUF(e.what());
		return TCL_ERROR;
	}catch(exception & e){
		RESULT_BUF(e.what());
		return TCL_ERROR;
	}catch(...){
		RESULT_LIT("An unknown error has occurred");
		return TCL_ERROR;
	}
	if(errmsg){
		RESULT_BUF(errmsg);
		return TCL_ERROR;
	}
    	return TCL_OK;
    }
    if(stream_count == 0 && !use_screen){
	RESULT_LIT("No stream count was specified");
	return TCL_ERROR;
    }
    GPU_GWAS * estimator;
    try{
        estimator = new GPU_GWAS(gpu_id_list, trait_list, phenotype_filename,evd_data_filename,\
             plink_filename, stream_count, precision, blockSize, defined_batch_size, verbose);
    }catch(exception & e){
        RESULT_BUF(e.what());
        return TCL_ERROR;
    }catch(...){
        RESULT_LIT("Unkown error occurred during GPU-GWAS setup");
        return TCL_ERROR;
    }
    
    int n_sets =  estimator->get_n_sets();
    if(n_sets == 0){
	RESULT_LIT("No viable data sets could be read");
    }
    if(run_calibration) n_sets = 1;
    int success = 0;
    if(verbose){
    	cout << "Starting GWAS Estimation\n";
    }
    for(int set = 0; set < n_sets; set++){

	if(!run_calibration){
		if(verbose){
			cout << "Processing trait set " << set + 1 << " out of " << n_sets << " trait sets\n";
		}
		int success;
       		 try{
          		  success = estimator->process_next_trait_set();
       	 	}catch(Solar_Trait_Reader_Exception & e){
            		RESULT_BUF(e.what());
           		 return TCL_ERROR;
         	}catch(GPU_Exception & e){
            		RESULT_BUF(e.what());
           		 return TCL_ERROR;
         	}catch(exception & e){
            		RESULT_BUF(e.what());
           		 return TCL_ERROR;
         	}catch(...){
            		string error_message = "Unknown error occurred while processing set " + to_string(set);
            		RESULT_BUF(error_message.c_str());
            		return TCL_ERROR;
         	 }
		if(success == 1) break;
	}else{
		chrono::milliseconds shortest_time;
		int best_blockSize;
		
		try{
		    shortest_time = estimator->calibrate_blockSize(0, best_blockSize, snp_percentage, defined_batch_size);
		
		}catch(Solar_Trait_Reader_Exception & e){
            		RESULT_BUF(e.what());
           		 return TCL_ERROR;
         	}catch(GPU_Exception & e){
            		RESULT_BUF(e.what());
           		 return TCL_ERROR;
         	}catch(exception & e){
            		RESULT_BUF(e.what());
           		 return TCL_ERROR;
         	}catch(...){
            		string error_message = "Unknown error occurred while processing set " + to_string(set);
            		RESULT_BUF(error_message.c_str());
            		return TCL_ERROR;
         	 }
		cout << "Shortest time: " << shortest_time.count() << " ms\n";
		cout << "Best Block Size: " << best_blockSize << endl;
	}
    }
    try{
        delete estimator;
    }catch(exception & e){
        RESULT_BUF(e.what());
        return TCL_ERROR;
    }catch(...){
        RESULT_LIT("Unknown error occurred during GPU-GWAS cleanup");
        return TCL_ERROR;
    }
    if(success == 1) return TCL_ERROR;
 
    return TCL_OK;
}
/*
extern "C" int  GPU_GWAS_Cmd(ClientData clientData, Tcl_Interp *interp,
                             int argc,const char *argv[]){

    RESULT_LIT("GPU GWAS is not availible on this distribution of solar-eclipse");
    return TCL_OK;
}*/


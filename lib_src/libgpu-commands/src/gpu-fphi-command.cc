#include <string>
#include "solar-trait-reader.h"
#include "solar.h"
#include <exception>
#include <vector>
#include <stdio.h>
#include "gpu-selection.h"
#include "gpu-exception.h"
#include "gpu-fphi-settings.h"
using namespace std;
vector<string> read_list_of_traits(const char * list_filename);
const char * run_gpu_fphi(const char * phenotype_filename, const char * output_filename, const char * evd_base_filename, const char * template_filename, vector<string> raw_trait_list, vector<int> gpu_id_list, \
			 int n_streams,  int defined_blockSize, int defined_batch_size, const bool use_covariates,const bool verbose);
			 
void display_gpu_fphi_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help gpu_fphi");
}			 
extern "C" int gpufphiCmd (ClientData clientData, Tcl_Interp *interp,
		int argc, const char *argv[]){
		
	const char * output_filename = 0;
	const char * template_filename = 0;
	const char * evd_base_filename = 0;
	int n_streams = 1;
	int defined_blockSize = -1;
	int defined_batch_size = -1;
	bool use_covariates = false;
	const char * list_filename = 0;
	bool select_all_gpus = false;
	bool verbose = false;
	const char * gpu_list_string = 0;
	for(int arg = 1; arg < argc; arg++){
		if(!StringCmp(argv[arg], "help", case_ins) || !StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins)){
			display_gpu_fphi_help(interp);
			return TCL_OK;
		}else if((!StringCmp(argv[arg], "-out", case_ins) || !StringCmp(argv[arg], "--out", case_ins) || !StringCmp(argv[arg], "--o", case_ins) ||\
			 !StringCmp(argv[arg], "-o", case_ins) || !StringCmp(argv[arg], "--output", case_ins) || !StringCmp(argv[arg], "-output", case_ins)) && arg + 1 < argc){
			 output_filename = argv[++arg];
		}else if((!StringCmp(argv[arg], "-nifti", case_ins) || !StringCmp(argv[arg], "--nifti", case_ins)) && (arg + 1 < argc || arg < argc) ){
			if(arg + 1 >= argc){
				RESULT_LIT("-nifti option requires the specification of a nifti template filename");
				return TCL_ERROR;
			}
			template_filename = argv[++arg];
		}else if((!StringCmp(argv[arg], "-evd_data", case_ins) || !StringCmp(argv[arg], "--evd_data", case_ins)) && arg + 1 < argc  ){

			evd_base_filename = argv[++arg];
		}else if((!StringCmp(argv[arg], "-n_streams", case_ins) || !StringCmp(argv[arg], "--n_streams", case_ins)) && arg + 1 < argc){ 
			n_streams = atoi(argv[++arg]);
			if(n_streams <= 0){
				RESULT_LIT("Stream number must be at least one");
				return TCL_ERROR;
			}
			if(n_streams > MAX_STREAM_COUNT){
				string error_message = "Maximum stream count is " + to_string(MAX_STREAM_COUNT);
				RESULT_BUF(error_message.c_str());
				return TCL_ERROR;
			}				
		}else if((!StringCmp(argv[arg], "-batch_size", case_ins) || !StringCmp(argv[arg], "--batch_size", case_ins)) && arg + 1 < argc){
			defined_batch_size = atoi(argv[++arg]);
			if(defined_batch_size <= 0){
				RESULT_LIT("Batch Size must be greater than zero");
				return TCL_ERROR;
			}
			if(defined_batch_size > MAX_BATCH_SIZE){
				string error_message = "Maximum batch size is " + to_string(MAX_BATCH_SIZE);
				RESULT_BUF(error_message.c_str());
				return TCL_ERROR;
			}
			
		}else if((!StringCmp(argv[arg], "-thread_size", case_ins) || !StringCmp(argv[arg], "--thread_size", case_ins)) && arg + 1 < argc){
			defined_blockSize = atoi(argv[++arg]);
			if(defined_blockSize != 32 && defined_blockSize != 64 && defined_blockSize != 128 && defined_blockSize != 256 && \
				defined_blockSize != 512 && defined_blockSize != 1024){
				RESULT_LIT("Thread Size must be 32,64,128,256,512 or 1024");
				return TCL_ERROR;
			}
		}else if(!StringCmp(argv[arg], "-use_covs", case_ins) || !StringCmp(argv[arg], "--use_covs", case_ins)){
			use_covariates = true;
		}else if((!StringCmp(argv[arg], "-list", case_ins) || !StringCmp(argv[arg], "--list", case_ins)) && arg + 1 < argc){
			list_filename = argv[++arg];
		}else if ((!StringCmp(argv[arg], "--gpus", case_ins) || !StringCmp(argv[arg], "-gpus", case_ins)) && arg + 1 < argc){
            
           		 gpu_list_string = argv[++arg];
            
        	}else if (!StringCmp(argv[arg], "--all_gpus", case_ins) || !StringCmp(argv[arg], "-all_gpus", case_ins)){
            
           	 	select_all_gpus = true;
            
       	 	}else if (!StringCmp(argv[arg], "--display-gpus", case_ins) || !StringCmp(argv[arg], "-display-gpus", case_ins)){
            
           	 	try{
           	 	    print_usable_devices();
           	 	 }catch(GPU_Exception & e){
           	 	 	RESULT_BUF(e.what());
           	 	 	return TCL_ERROR;
           	 	 }
           	 	 return TCL_OK;
            
       	 	}else if (!StringCmp(argv[arg], "--verbose", case_ins) || !StringCmp(argv[arg], "-verbose", case_ins)){
            
           	 	verbose = true;
            
       	 	}else{
       	 		string error_message = "Invalid option was entered: " + string(argv[arg]);
       	 		RESULT_BUF(error_message.c_str());
       	 		return TCL_ERROR;
       	 	}
       	 }
       	 	
       	 if(!list_filename){
       	 	RESULT_LIT("Please specify a file containing a trait list with --list argument");
       	 	return TCL_ERROR;
       	 }
       	 if(!output_filename){
       	 	RESULT_LIT("Since no output filename was specified the default output filename or base filename will be gpu-fphi.out or -gpu-fphi.nii.gz");
       	 	if(!template_filename){
       	 		output_filename = "gpu-fphi.out";
       	 	}else{
       	 		output_filename = "gpu-fphi";
       	 	}
       	 }
       	 if(!loadedPed()){
       	 	RESULT_LIT("No pedigree is currently loaded");
       	 	return TCL_ERROR;
       	 }
       	 string str_phenotype_filename = Phenotypes::filenames();
       	 if(str_phenotype_filename.length() == 0){
       	 	RESULT_LIT("No phenotype file has been loaded");
       	 	return TCL_ERROR;
       	 }
       	 const char * phenotype_filename = str_phenotype_filename.c_str();
	if(!evd_base_filename){
    		try {
            		load_phi2_matrix(interp);
   		}catch(Solar_Trait_Reader_Exception & e){
       	    		RESULT_BUF(e.what());
       	     		return TCL_ERROR;
   		}catch(...){
             		RESULT_LIT("Unkown error occurred when loading phi2 matrix");
             		return TCL_ERROR;
    		} 
	}
    		      	 	
       	 vector<int> gpu_id_list;
       	 try{
       	 	gpu_id_list = get_gpu_id_list(gpu_list_string, select_all_gpus);
       	 }catch(GPU_Exception & e){
       	 	RESULT_BUF(e.what());
       	 	return TCL_ERROR;
       	 }
       	 vector<string> trait_list = read_list_of_traits(list_filename);
       	 if(trait_list.size() == 0){
       	 	RESULT_LIT("Unable to assemble list of traits from specifed list filename");
       	 	return TCL_ERROR;
       	 }
       	 const char * error_message = 0;
       	 error_message = run_gpu_fphi(phenotype_filename, output_filename, evd_base_filename, template_filename, trait_list, gpu_id_list, \
			  n_streams,   defined_blockSize,  defined_batch_size,  use_covariates, verbose);
	if(error_message){
		RESULT_BUF(error_message);
		return TCL_ERROR;
	}
		
	return TCL_OK;
}
		 
	

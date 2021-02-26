#include "solar.h"
#include "safelib.h"
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <stdio.h>
#include "plinkio.h"
#include <algorithm>

using namespace std;


static void print_help(Tcl_Interp * interp){
  Solar_Eval(interp,"help plink_converter\n");
}


extern "C" int convert_plink_to_csv(struct pio_file_t * input, const char * output_basename, const size_t max_per_file, const int bin_format, const int solar_format,  const int per_chromo, char ** errmsg);


extern "C" int RunPlinkConverter(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[])
{
    const char * plink_filename = 0;
    const char * output_basename = 0;
    int per_chromo = 0;
    int bin_format = 0; 
    int solar_format = 0;
    size_t max_per_file = 0;
    for(int arg = 1; arg < argc; arg++){
        if(!StringCmp("--help", argv[arg], case_ins) || !StringCmp("-help", argv[arg], case_ins)\
           || !StringCmp("help", argv[arg], case_ins)){
            print_help(interp);
            return TCL_OK;
        }else if((!StringCmp("--i", argv[arg], case_ins) || !StringCmp("-i", argv[arg], case_ins))\
                 && arg + 1 < argc){
            plink_filename = argv[++arg];
        }else if((!StringCmp("--o", argv[arg], case_ins) || !StringCmp("-o", argv[arg], case_ins)\
                 || !StringCmp("--out", argv[arg], case_ins) || !StringCmp("-out", argv[arg], case_ins))\
                    && arg + 1 < argc){
            output_basename = argv[++arg];
        }else if(!StringCmp("--bin", argv[arg], case_ins) || !StringCmp("-bin", argv[arg], case_ins)){
            bin_format = 1;
        }else if(!StringCmp("--per-chromo", argv[arg], case_ins) || !StringCmp("-per-chromo", argv[arg], case_ins)){
            per_chromo = 1;
        }else if((!StringCmp("--max", argv[arg], case_ins) || !StringCmp("-max", argv[arg], case_ins))\
                 && arg + 1 < argc){
            max_per_file = strtol(argv[++arg], NULL, 10);
        }else if(!StringCmp("--solar", argv[arg], case_ins) || !StringCmp("-solar", argv[arg], case_ins)){
		solar_format = 1;
	} else{
            RESULT_LIT("Invalid argument was entered");
            return TCL_ERROR;
        }
    }
    
    
    if(!plink_filename){
        RESULT_LIT("A plink file set must be specified with -i <base filename> argument");
        return TCL_ERROR;
    }
    
    
    if(!output_basename){
        RESULT_LIT("The output base filename must be specified with -o <base filename> argument");
        return TCL_ERROR;
    }
    
    
    pio_file_t input;

 // pio_transpose(argv[1], argv[2]);
    pio_status_t status = pio_open(&input, plink_filename);

  if(status != PIO_OK){
      pio_close(&input);
      if(status == P_FAM_IO_ERROR){
	  RESULT_LIT("An error occurred while loading .fam file");
        
      }else if (status == P_BIM_IO_ERROR){
	  RESULT_LIT("An error occurred while loading .bim file");

      }else if (status == P_BED_IO_ERROR){
	  RESULT_LIT("An error occurred while loading .bed file");
      }else{
          RESULT_LIT("An error occurred while loading the plink file set");
      }
      return TCL_ERROR;
    }
    
    
    

    int version = input.bed_file.header.version;
    if (input.bed_file.header.snp_order == BED_UNKNOWN_ORDER){
      pio_close(&input);
      RESULT_LIT("The storage order of the .bed file could not be discerned");
      return TCL_ERROR;
      
    }else if (input.bed_file.header.snp_order == BED_ONE_SAMPLE_PER_ROW){
      pio_close(&input);
      printf("In order to read efficiently, the transpose of %s must be written.\n", plink_filename);
      printf("The transpose will be written to %s_trans\n", plink_filename);
      string trans_input_basename = string(plink_filename) + "_trans";
      pio_status_t status = pio_transpose(plink_filename, trans_input_basename.c_str());
      if(status != PIO_OK){
	  RESULT_LIT("Error in creating transpose");
          return TCL_ERROR;
      }

      status = pio_open(&input, trans_input_basename.c_str());

      if(status != PIO_OK){
          RESULT_LIT("Error in opening transpose.");
          return TCL_ERROR;
      }
    }
    char * errmsg = 0;
    convert_plink_to_csv(&input, output_basename,  max_per_file,  bin_format, solar_format,  per_chromo, &errmsg);
    return TCL_OK;
}



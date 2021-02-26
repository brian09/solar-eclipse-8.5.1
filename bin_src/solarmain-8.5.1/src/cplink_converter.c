#include <stdio.h>
#include <math.h>
#include "plinkio.h"
#include "safelib.h"
/*
static int count_chromosomes(struct pio_file_t * input, size_t n_snps){
	int n_chromo = 1;
	struct pio_locus_t * locus = bim_get_locus(&input->bim_file, 0);
	unsigned char chromosome = locus->chromosome;
	for(size_t snp = 1; snp < n_snps; snp++){
		locus =  bim_get_locus(&input->bim_file, snp);
		if(locus->chromosome != chromosome) {
			n_chromo++;
			chromosome = locus->chromosome;
		}
	}
	return n_chromo;
}
*/
static void write_plink_data(struct pio_file_t * input, snp_t *  snp_data, const char * output_filename, const size_t start, \
                            const size_t n_snps, const int bin_format, const int solar_format){
    const size_t n_samples = input->bed_file.header.num_samples;
    FILE * output_file = fopen(output_filename, "w");
    fprintf(output_file, "id,fid,sex");
    struct pio_locus_t * locus;
    if(bin_format || solar_format){
        FILE * header_file;
        char header_filename[strlen(output_filename) + 12];
        sprintf(header_filename, "%s.header.csv", output_filename);
        header_file = fopen(header_filename, "w");
        fprintf(header_file, "SNP,0,1,2\n");
        for(size_t snp = start; snp < start + n_snps; snp++){
            locus = bim_get_locus(&input->bim_file, snp);
            fprintf(output_file,",snp_%s",locus->name);
            
            fprintf(header_file,"snp_%s,%s%s,%s%s,%s%s\n", locus->name , locus->allele1 , locus->allele1 \
                    , locus->allele1, locus->allele2 ,  locus->allele2 , locus->allele2);
        }
        fprintf(output_file, "\n");
        fclose(header_file);
    }else{
        for(size_t snp = start; snp < start + n_snps; snp++){
            locus = bim_get_locus(&input->bim_file, snp);
            fprintf(output_file,",snp_%s",locus->name);
        }
        fprintf(output_file, "\n");
    }
    struct pio_sample_t * sample;
    for(size_t subject_index = 0; subject_index < n_samples; subject_index++){
        sample = fam_get_sample(&input->fam_file, subject_index);
       
        if(sample->sex == PIO_MALE){
            fprintf(output_file,"%s,%s,M",sample->iid, sample->fid);
        }else if(sample->sex == PIO_FEMALE){
            fprintf(output_file,"%s,%s,F",sample->iid, sample->fid);
        }else{
            fprintf(output_file,"%s,%s,",sample->iid, sample->fid);
        }
     
        snp_t value;
        if(bin_format){
            for(size_t snp = start; snp < start + n_snps; snp++){
                value = snp_data[subject_index + snp*n_samples];
                if(value != 3)
                    fprintf(output_file,",%u",(unsigned)value);
                else
                    fprintf(output_file,",");
            }
        }else if(solar_format){
		for(size_t snp = start; snp < start + n_snps; snp++){
		value = snp_data[subject_index + snp*n_samples];	
		if(value != 3){
			switch(value){
				case 0:
					fprintf(output_file,",1/1");
					break;
				case 1: 
					fprintf(output_file,",1/2");
                                        break;
				case 2:
					fprintf(output_file,",2/2");
                                        break;	
				}
			}else{

			 fprintf(output_file,",");
			}
			}		
           }else{
            for(size_t snp = start; snp < start + n_snps; snp++){
                value = snp_data[subject_index + snp*n_samples];
                if(value != 3){
                    locus = bim_get_locus(&input->bim_file, snp);
                    switch(value){
                        case 0:
                            fprintf(output_file, ",%s/%s",locus->allele1, locus->allele1);
                            break;
                        case 1:
                            fprintf(output_file, ",%s/%s",locus->allele1, locus->allele2);
                            break;
                        case 2:
                            fprintf(output_file, ",%s/%s",locus->allele2, locus->allele2);
                            break;
                        default:
                            fprintf(output_file, ",");
                            break;
                    }
                }else{
                    fprintf(output_file, ",");
                }
            }
        }
        
        fprintf(output_file, "\n");
        
    }
    
    fclose(output_file);
}

int convert_plink_to_csv(struct pio_file_t * input, const char * output_basename, const size_t max_per_file, const int bin_format,const int solar_format, const int per_chromo, char ** errmsg){
    const size_t n_snps = input->bed_file.header.num_loci;
    const size_t n_samples = input->bed_file.header.num_samples;
    
    snp_t * snp_data = (snp_t*)calloc(n_snps*n_samples, sizeof(snp_t));
    if(snp_data == NULL){
        errmsg = "Could not allocate required memory";
        return 1;
    }
    snp_t * buffer = (snp_t*)calloc(n_samples, sizeof(snp_t));
    if(buffer == NULL){
        free(snp_data);
        errmsg = "Could not allocate required memory";
        return 1;
    }
    snp_t * buffer_ptr;
    snp_t * snp_data_ptr = snp_data;
    for(size_t snp = 0 ; snp < n_snps; snp++){
        pio_next_row(input, buffer);
        buffer_ptr = buffer;
        for(size_t row = 0; row < n_samples; row++){
            (*snp_data_ptr++) = (*buffer_ptr++);
        }
    }
    free(buffer);
    struct pio_locus_t * locus;
    if(per_chromo && max_per_file == 0){
       // locus = bim_get_locus(&input->bim_file, 0);
       // unsigned char current_chromosome = locus->chromosome;
        size_t start = 0;
        size_t batch_size;
        size_t snp_index;
        while(start < n_snps){
            locus = bim_get_locus(&input->bim_file, start);
            unsigned char current_chromosome = locus->chromosome;
            snp_index = start;
            while(locus->chromosome == current_chromosome){
                if(++snp_index < n_snps){
                    locus = bim_get_locus(&input->bim_file, snp_index);
                }else{
                    break;
                }
            }
            batch_size = snp_index - start;
            char output_filename[strlen(output_basename) + 12];
            sprintf(output_filename, "%s.chr%u.csv", output_basename, (unsigned) current_chromosome);
            write_plink_data(input, snp_data, (const char *)output_filename, start, \
                             batch_size,  bin_format, solar_format);
            start += batch_size;
        }
        
    }else if(per_chromo){
        size_t start = 0;
        size_t batch_size;
        size_t snp_index;
        size_t file_index;
        size_t next_file_index = 0;
        while(start < n_snps){
            file_index = next_file_index;
            locus = bim_get_locus(&input->bim_file, start);
            unsigned char current_chromosome = locus->chromosome;
            snp_index = start;
            while(locus->chromosome == current_chromosome){
                snp_index++;
                if(snp_index < n_snps && (snp_index - start) < max_per_file){
                    locus = bim_get_locus(&input->bim_file, snp_index);
                }else{
                    break;
                }
            }
            if(snp_index == n_snps){
                batch_size = n_snps - start;
            }else if(snp_index - start == max_per_file){
                batch_size = max_per_file;
                next_file_index = file_index + 1;
            }else{
                next_file_index = 0;
                batch_size = snp_index - start;
            }
            char output_filename[strlen(output_basename) + 20];
            sprintf(output_filename, "%s.chr%u.%u.csv", output_basename, (unsigned) current_chromosome, file_index);
            write_plink_data(input, snp_data, (const char *)output_filename, start, \
                             batch_size,  bin_format, solar_format);
            start += batch_size;
        }
    }else if(max_per_file != 0){
        size_t start = 0;
        size_t batch_size;
        size_t file_index = 0;
        while(start < n_snps){
            if(start + max_per_file < n_snps){
                batch_size = max_per_file;
            }else{
                batch_size = n_snps - start;
            }
            

            char output_filename[strlen(output_basename) + 16];
            sprintf(output_filename, "%s.%u.csv", output_basename,  file_index++);
            write_plink_data(input, snp_data, (const char *)output_filename, start, \
                             batch_size,  bin_format,solar_format);
            start += batch_size;
        }
    }else{
        char output_filename[strlen(output_basename) + 6];
        sprintf(output_filename, "%s.csv", output_basename);
        write_plink_data(input, snp_data, (const char * )output_filename, 0, n_snps, bin_format, solar_format);
    }
    
    free(snp_data);
    
    return 0;
}
/*
void read_and_write_to_csv(struct pio_file_t  * input, const char* output_basename,
			  size_t col_size, size_t max_per_file, int bin_format, int per_chromo){

  unsigned int row_size = input->bed_file.header.num_samples;
  int num_files;

  if(per_chromo){
	  num_files = count_chromosomes(input, col_size);
	  max_per_file = col_size;
  }else if (max_per_file == 0){
      max_per_file = col_size;
      num_files = 1;
  }else{
      num_files = ceil((double)col_size/(double)max_per_file);
  }



  snp_t buffer[row_size];
  struct pio_locus_t * locus;
  FILE *  oss[num_files];

  FILE * oss_header[num_files];
  printf("file count: %i \n", num_files);
  if(per_chromo){
	 int chromo = 0;  
      for(unsigned int i = 0 ; i < num_files; i++){
		char output_filename[200];
      // = output_basename + ".csv";
		sprintf(output_filename, "%s.chr%i.csv",output_basename, chromo);
		oss[i] = fopen(output_filename, "w");
              if(oss[i] == 0) {
				  printf("Error opening file for write, try decreasing the file count.\n");
				  return;
				}		
          if(bin_format){
			  char output_filename_header[200];// = output_basename + ".header.csv";
			  sprintf(output_filename_header, "%s.chr%i.header.csv",output_basename, chromo);              
              oss_header[i] = fopen(output_filename_header, "w");
              if(oss_header[i] == 0) {
				  printf("Error opening file for write, try decreasing the file count.\n");
				  return;
				}              
          }
		
		chromo++;
      }	  

  }	 else if (num_files == 1){
      char output_filename[200];
      // = output_basename + ".csv";
      sprintf(output_filename, "%s.csv",output_basename);
      oss[0] = fopen(output_filename, "w");//.open(output_filename.c_str(), ios::out);
      if(bin_format){
          char output_filename_header[200];// = output_basename + ".header.csv";
          sprintf(output_filename_header, "%s.header.csv",output_basename);
          oss_header[0] = fopen(output_filename_header, "w");//.open(output_filename_header.c_str(), ios::out);
      }
  }else{

      for(unsigned int i = 0 ; i < num_files; i++){
		char output_filename[200];
      // = output_basename + ".csv";
		sprintf(output_filename, "%s_%i.csv",output_basename, i);
		oss[i] = fopen(output_filename, "w");
              if(oss[i] == 0) {
				  printf("Error opening file for write, try decreasing the file count.\n");
				  return;
				}		
          if(bin_format){
			
			  char output_filename_header[200];// = output_basename + ".header.csv";
			  sprintf(output_filename_header, "%s_%i.header.csv",output_basename, i);              
              oss_header[i] = fopen(output_filename_header, "w");
              if(oss_header[i] == 0) {
				  printf("Error opening file for write, try decreasing the file count.\n");
				  return;
				}
          }
		

      }


  }




  fprintf(oss[0],"id,fid,sex,phenotype,");//oss[0] << "id,fid,sex,phenotype,";
  if(bin_format){
	fprintf(oss_header[0],"SNP,0,1,2\n");
  }
  int file = 0;
  short bits;
  locus = bim_get_locus(&input->bim_file, 0);
  unsigned char current_chromo = locus->chromosome;
  unsigned char new_chromo = current_chromo;
  int changed_chromo = 0;

  for(unsigned int i = 0 ; i < col_size ; i++){
	  
      if((max_per_file*(file + 1)) == i || (new_chromo != current_chromo && per_chromo)){
		
	  file++;
	  
	  current_chromo = new_chromo;
	  fprintf(oss[file],"id,fid,sex,phenotype,");
          if(bin_format) {

			fprintf(oss_header[file], "SNP,0,1,2\n");
          }
      }	  
	  if(i != col_size -1 && per_chromo){
		locus = bim_get_locus(&input->bim_file, i + 1);		  
		   new_chromo = locus->chromosome;
	   }	  
		locus = bim_get_locus(&input->bim_file, i);		
		
	  	  fprintf(oss[file], "snp_%s", locus->name);
          if(bin_format)
			fprintf(oss_header[file],"snp_%s,%s%s,%s%s,%s%s\n", locus->name , locus->allele1 , locus->allele1 \
              , locus->allele1, locus->allele2 ,  locus->allele2 , locus->allele2);		  

			if(i == ((file + 1)*(max_per_file) - 1) || (i == col_size - 1)|| (new_chromo != current_chromo && per_chromo)){
    	       fprintf(oss[file],"\n");
    	   }else{
    	       fprintf(oss[file],",");
    	   }      
  }

  struct pio_sample_t * sample;
  for(size_t row = 0 ; row < row_size ; row++){


      sample = fam_get_sample(&input->fam_file, row);
      file = 0;
      fprintf(oss[file], "%s,%s,%i,%f,", sample->iid, sample->fid, sample->sex, sample->phenotype);
	  locus = bim_get_locus(&input->bim_file, 0);     
      current_chromo = locus->chromosome;
      new_chromo = current_chromo;
      for(size_t col = 0 ; col < col_size ; col++){
		   pio_status_t status = pio_next_row(input, buffer);

			locus = bim_get_locus(&input->bim_file, col); 
			
		   if((max_per_file*(file + 1)) == col || (new_chromo != current_chromo && per_chromo)){
			file++;
			sample = fam_get_sample(&input->fam_file, row);
			current_chromo = new_chromo;
			fprintf(oss[file], "%s,%s,%i,%f,", sample->iid, sample->fid, sample->sex, sample->phenotype);
		   }
		   bits = buffer[row];
	
	   
		//   locus = bim_get_locus(&input->bim_file, col);
    	   if((!StringCmp(locus->allele1, "0", 0)
    	       || !StringCmp(locus->allele2, "0", 0))  && bits == 1){
    	       bits = 3;
    	   }else if((!StringCmp(locus->allele1, "0", 0) )
    		  && bits == 0){
    		  bits = 3;
    	   }else if((!StringCmp(locus->allele2, "0", 0))
 	     && bits == 2){
    	       bits = 3;
 	   }


    	   if(bits == 0){

    	       if(bin_format){
    		   fprintf(oss[file], "%i", bits);
    	       }else{
				  fprintf(oss[file], "%s%s", locus->allele1,locus->allele1);
    	       }
    	   }else if (bits == 1){
    	       if(bin_format){
    		   fprintf(oss[file], "%i", bits);
    	       }else{
			   fprintf(oss[file], "%s%s", locus->allele1,locus->allele2);
  
    	       }
    	   }else if (bits == 2){
    	       if(bin_format){
    		   fprintf(oss[file], "%i", bits);
    	       }else{
				 fprintf(oss[file], "%s%s", locus->allele2,locus->allele2);	   
    	       }

			}
		
	  if(col != col_size -1 && per_chromo){
		locus = bim_get_locus(&input->bim_file, col + 1);		  
		   new_chromo = locus->chromosome;
	   }				
		
    	   if(col == ((file + 1)*(max_per_file) - 1) || (col == col_size - 1)|| (new_chromo != current_chromo && per_chromo)){
    	       fprintf(oss[file],"\n");
    	   }else{
    	       fprintf(oss[file],",");
    	   }
      }
      pio_reset_row(input);
  }

  for(int f = 0 ; f < num_files ; f++){
	  fclose(oss[f]);
      if(bin_format)
          fclose(oss_header[f]);
  }

  pio_close(input);

  

}*/

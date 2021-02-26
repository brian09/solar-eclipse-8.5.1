#include <iostream>
#include "solar.h"
#include <fstream>
#include <string>
#include <vector>
#include "plinkio.h"
#include <algorithm>
using namespace std;
extern bool loadedPed ();
extern Pedigree *currentPed;
//extern Pedigree *currentPed;
extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);
static void print_evdCmd_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help create_evd");
}
static void calculate_eigenvectors_and_eigenvalues (double * phi2, double * eigenvectors, double * eigenvalues, int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}
static void write_evd_to_output(vector<string> ids, double * eigenvectors, double * eigenvalues, const char * base_output_filename){
	string filename = string(base_output_filename) + ".ids";
	ofstream id_stream(filename.c_str());
	id_stream << ids[0];
	for(int i = 1; i < ids.size(); i++){
		id_stream << " " << ids[i];
	}
	id_stream.close();
	filename = string(base_output_filename) + ".eigenvalues";

	ofstream eigenvalue_stream(filename.c_str());
	eigenvalue_stream << eigenvalues[0];
	for(int i = 1; i < ids.size(); i++){
		eigenvalue_stream << " " << eigenvalues[i];
	}
	eigenvalue_stream.close();

	filename = string(base_output_filename) + ".eigenvectors";

	ofstream eigenvector_stream(filename.c_str());
	for(int col = 0; col < ids.size(); col++){
		for(int row = 0 ; row < ids.size(); row++){
			eigenvector_stream << eigenvectors[col*ids.size() + row] << " ";
		}
	}
	eigenvector_stream.close();	
}
static vector<string> get_covariate_term_list(){
	Covariate * c;
	
    	vector<string> covariate_terms;
   	  
    	for (int i = 0;( c = Covariate::index(i)); i++){
    		CovariateTerm * cov_term;
        
        	for(cov_term = c->terms(); cov_term; cov_term = cov_term->next){
            	bool found = false;
            
            	for(vector<string>::iterator cov_iter = covariate_terms.begin(); cov_iter != covariate_terms.end(); cov_iter++){
                	if(!StringCmp(cov_term->name, cov_iter->c_str(), case_ins)){
                    		found = true;
                   		 break;
                	}
            	}
            	if(!found){
                	covariate_terms.push_back(string(cov_term->name));
            	}
        	}
      
    	}
    	
    	return covariate_terms;

}
static const char * generate_evd_data(Tcl_Interp * interp, string trait_name, const char * phenotype_filename, const char * base_output_filename, const char * plink_filename, const bool use_covariates){
    vector<string> phenotype_ids;
    vector<string> covariate_terms;
    if(use_covariates){
    	covariate_terms = get_covariate_term_list();
    	if(covariate_terms.size() == 0){
    		cout << "No covariates were found proceeding without the use of covariate data\n";
    	}
    }
    const char * errmsg = 0;
    SolarFile * reader =  SolarFile::open("Create EVD", phenotype_filename, &errmsg);
    if(errmsg){
        return errmsg;
    }
    
    reader->start_setup(&errmsg);
    if(errmsg){
        return errmsg;
    }
    reader->setup("id", &errmsg);
    if(errmsg){
        return errmsg;
    }
    reader->setup(trait_name.c_str(), &errmsg);
    if(errmsg){
        return errmsg;
    }
    
    if(covariate_terms.size() != 0){
    	for(int index = 0; index < covariate_terms.size(); index++){
    		string term = covariate_terms[index];
     		reader->setup(term.c_str(), &errmsg);
   		 if(errmsg){
        		return errmsg;
    		}   		
    	}
    }
    vector<string> skip_ids;
    char ** file_data;

    while (0 != (file_data = reader->get (&errmsg))){
        string current_id = file_data[0];
        bool skip_row = false;
       // all_ids.push_back(string(file_data[0]));
	if(!StringCmp(file_data[1], 0, case_ins)){
		continue;
	}
	if(covariate_terms.size() != 0){
		for(int index = 2; index < covariate_terms.size() + 2; index++){
			if(!StringCmp(file_data[index], 0, case_ins)){
				skip_row = true;
				break;
			}
		}
	}
	if(!skip_row){
		phenotype_ids.push_back(current_id);
	}
	
    }
    delete reader;
    vector<string> shared_id_list;
    if(plink_filename){
    	pio_file_t * plink_file = new pio_file_t;
	if (pio_open(plink_file, plink_filename) != PIO_OK){
		return "Error opening plink file";
	}
	pio_sample_t * sample;
	vector<string> plink_ids;
	const unsigned num_plink_samples = pio_num_samples(plink_file);
	for(unsigned i = 0; i < num_plink_samples; i++){
		sample = pio_get_sample(plink_file, i);
		plink_ids.push_back(string(sample->iid));
	}
	
	for(vector<string>::iterator iter = phenotype_ids.begin(); iter != phenotype_ids.end(); iter++){

		vector<string>::iterator find_iter = find(plink_ids.begin(), plink_ids.end(), *iter);
		if(find_iter != plink_ids.end()){
			shared_id_list.push_back(*iter);
			plink_ids.erase(find_iter);
		}
	}
	pio_close(plink_file);
	delete plink_file;	  
    }else{
    	shared_id_list = phenotype_ids;
    }	

    vector<int> ibdids;
    vector<string> ids;
    SolarFile *  ped_reader = SolarFile::open("Create EVD", "pedindex.out" , &errmsg);

    if(errmsg) return errmsg;
     ped_reader->start_setup(&errmsg);
    if(errmsg){
        return errmsg;
    }   
    ped_reader->setup("id", &errmsg);
    if(errmsg){
        return errmsg;
    }
    int ibdid = 1;
    while (0 != (file_data = ped_reader->get (&errmsg))){
        string current_id = file_data[0];
        vector<string>::iterator find_iter  = find(shared_id_list.begin(), shared_id_list.end(), current_id);
	if(find_iter != shared_id_list.end()){
		ids.push_back(current_id);
		ibdids.push_back(ibdid);
		shared_id_list.erase(find_iter);
	}
	ibdid++;
    } 

    Matrix * solar_phi2 = 0;
    solar_phi2 = Matrix::find("phi2");
    if (!solar_phi2) {
        Solar_Eval(interp, "matrix load phi2.gz phi2");
        solar_phi2 = Matrix::find("phi2");
        if(!solar_phi2){
            return "Matrix could not be loaded from phi2.gz";
        }
    }

    double * phi2 = new double[ids.size()*ids.size()];	
    for(int col = 0; col < ids.size(); col++){
	phi2[col*ids.size() + col] = solar_phi2->get(ibdids[col], ibdids[col]);
	for(int row = col + 1; row < ids.size(); row++){
	    phi2[col*ids.size() + row] =  phi2[row*ids.size() + col] = solar_phi2->get(ibdids[row], ibdids[col]);
	} 
    }
    double * eigenvectors = new double[ids.size()*ids.size()]; 
    double * eigenvalues = new double[ids.size()];

    calculate_eigenvectors_and_eigenvalues (phi2, eigenvectors, eigenvalues, ids.size());

    delete [] phi2;
    
    write_evd_to_output(ids, eigenvectors,  eigenvalues,  base_output_filename);
    string notes_filename = string(base_output_filename) + ".notes";
    ofstream notes_output_stream(notes_filename.c_str());
    notes_output_stream << "Pedigree filename that produced the phi2 matrix used for EVD computation: " << currentPed->filename() << endl;
    notes_output_stream << "Number of IDs: " << ids.size() << endl;
    notes_output_stream << "Phenotype filename used for ID selection: " << phenotype_filename << endl;
    notes_output_stream << "Trait used for ID selection: " << trait_name << endl;
    if(plink_filename){
    	notes_output_stream << "PLINK file set name used for ID selection: " << plink_filename << endl;
    }    
    if(covariate_terms.size() != 0){
    	notes_output_stream << "Phenotype file fields of the selected covariates used for ID selection:";
    	for(int index =0; index < covariate_terms.size(); index++){
    		notes_output_stream << " " << covariate_terms[index];
    	}
    	notes_output_stream << endl;
    }

    
    notes_output_stream.close();

    delete [] eigenvectors;
    delete [] eigenvalues;  
    return 0; 
}
extern "C" int create_evdCmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){

	const char * plink_filename = 0;
	const char * base_output_filename = 0;
  	bool use_covariates = false;
   	 for(unsigned arg = 1; arg < argc; arg++){
		if(!StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins) || !StringCmp(argv[arg], "help", case_ins)){
			print_evdCmd_help(interp);
			return TCL_OK;
		}else if((!StringCmp(argv[arg], "-plink", case_ins) || !StringCmp(argv[arg], "--plink", case_ins)) && arg + 1 < argc){
			plink_filename = argv[++arg];
		}else if((!StringCmp(argv[arg], "-use_covs", case_ins) || !StringCmp(argv[arg], "--use_covs", case_ins))){
			use_covariates = true;
		}else if ((!StringCmp(argv[arg], "-o", case_ins) || !StringCmp(argv[arg], "--o", case_ins) || !StringCmp(argv[arg], "-out", case_ins)\
			  || !StringCmp(argv[arg], "--out", case_ins)) && arg + 1 < argc){
			base_output_filename = argv[++arg];
		}else{
			RESULT_LIT("Invalid argument entered");
			return TCL_ERROR;
		}
	}
    	if(Trait::Number_Of() == 0){
        	RESULT_LIT("No trait has been selected");
       		return TCL_ERROR;
   	}
  	if(!loadedPed()){
        	RESULT_LIT("No pedigree has been loaded");
       		return TCL_ERROR;
    	}
	if(base_output_filename == 0){
		RESULT_LIT("Please enter a base output filename with -o");
		return TCL_ERROR;
	}
	const char * phenotype_filename = 0;
        phenotype_filename = Phenotypes::filenames();
 	if(string(phenotype_filename).length() == 0){
		RESULT_LIT("No phenotype file is currently loaded");
		return TCL_ERROR;
	}
 
	string trait_name = string(Trait::Name(0));
	const char * errmsg = 0;
	errmsg = generate_evd_data(interp, trait_name,phenotype_filename,  base_output_filename, plink_filename, use_covariates);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	return TCL_OK;
}

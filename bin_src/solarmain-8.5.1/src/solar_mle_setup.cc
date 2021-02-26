#include "solar_mle_setup.h"
#include <unordered_map>
#include <vector>
#include "solar.h"
#include <algorithm>
#include <queue>
#include <string>
#include <Eigen/Dense>
#include <utility>
#include <iterator>
#include <iostream>
#include <functional>
#include <fstream>
using namespace std;


extern char* inor_match (char* string);

extern char* zscor_match (char* string);

static solar_mle_setup * CURRENT_LOADER = 0;
double eval_blank_expression(int data_index){
    return MISSING_PHENOTYPE;
}
double eval_inorm_expression(int data_index){
    return  CURRENT_LOADER->inorm_map[CURRENT_LOADER->ids[CURRENT_LOADER->id_index]].at(CURRENT_LOADER->inorm_indices[data_index]);
}

double eval_zscore_expression(int data_index){
    return  (CURRENT_LOADER->term_map[CURRENT_LOADER->ids[CURRENT_LOADER->id_index]].at(CURRENT_LOADER->zscore_indices[data_index]) - CURRENT_LOADER->zscore_means[data_index])/CURRENT_LOADER->zscore_sd[data_index];
}

double eval_default_expression(int data_index){
    return CURRENT_LOADER->term_map[CURRENT_LOADER->ids[CURRENT_LOADER->id_index]].at(data_index);

}

extern "C" void cdfnor_ (int* which, double* p, double* q, double* x, double* mean, double* sd, int* status, double* bound);

extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);

solar_mle_setup::~solar_mle_setup(){
    interp = 0;
    zscore_means.clear();
    zscore_sd.clear();
    inorm_indices.clear();
    zscore_indices.clear();
    ibdids.clear();
    phenotype_names.clear();
    delete _Context;
    for(int i = 0 ;i < n_expressions; i++){
      delete expressions[i];
    }
    delete expressions;
    phenotype_filename.clear();
    ids.clear();
    term_map.clear();
    inorm_map.clear();
    expression_names.clear();
}

solar_mle_setup::solar_mle_setup(vector<string> input_terms, const char * filename,Tcl_Interp * input_interp,  bool ped_switch){
    if(input_terms.size() == 0){
       
       throw Misc_Error("Misc Error: blank list of expressions");
    }
    interp = input_interp;
    _Context = new Context;
    _Context->add("blank", eval_blank_expression, 0);
    expressions = new Expression*[input_terms.size()];
    n_expressions = input_terms.size();
    expression_names = input_terms;
    int result = TCL_OK;
    for(int index = 0; index < input_terms.size(); index++){
       
        result = parse_expression(input_terms[index], index);
        if(result == TCL_ERROR){
            string error_message = "Parse Expression Error with term: " + input_terms[index];
            throw Parse_Expression_Error(error_message);
        }
    }
   
    phenotype_filename = filename;
    use_ped = ped_switch;
    result = read_phenotype_and_fill_maps();
    if(result == TCL_ERROR){
    	string error_message = "Solar File Error filename: " + string(filename);;
        throw Solar_File_Error(error_message);
    }
    n_subjects = ids.size();
    compute_inorm_and_zscore();
   
}

vector<string> solar_mle_setup::get_ids(){
    return ids;
}

vector<string> solar_mle_setup::get_expression_names(){
    return expression_names;
}

void solar_mle_setup::create_output_matrix(){
    id_index = 0;
    Eigen::MatrixXd output_matrix_temp(n_subjects, n_expressions);
    CURRENT_LOADER = this;
    for(int row_idx = 0; row_idx < n_subjects; row_idx++){
        for(int col_idx = 0; col_idx < n_expressions; col_idx++){
            try{
                output_matrix_temp(row_idx, col_idx) = expressions[col_idx]->eval();

            }catch(...){
                CURRENT_LOADER = 0;
                string error_message = "Expression Eval Error id: " + ids[row_idx] + " term: " + expression_names[col_idx] + " row_index: " + to_string(row_idx) + " col_index: " + to_string(col_idx);  
                throw Expression_Eval_Error(error_message);
           }
        }
        id_index++;
    }
   CURRENT_LOADER = 0;
   output_matrix = output_matrix_temp;
}

Eigen::MatrixXd solar_mle_setup::return_output_matrix(){
    if(!output_matrix_created){
	create_output_matrix();
        output_matrix_created = true;
    }
    return output_matrix;
}

void solar_mle_setup::calculate_eigenvectors_and_eigenvalues (double * phi2, double * eigenvalues, double * eigenvectors ,int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}
void solar_mle_setup::compute_eigen_data(){
    Matrix* pm2;
    
    pm2 = Matrix::find("phi2");
    if (!pm2) {
        Solar_Eval(interp, "matrix load phi2.gz phi2");
        pm2 = Matrix::find("phi2");
        if(!pm2){
            throw Misc_Error("Misc Error: Phi2 matrix could not be loaded");
        }
    }
    phi2_matrix.resize(n_subjects, n_subjects);
    double * phi2 = new double[n_subjects*n_subjects];
    double phi2_value;

    for(int col = 0; col < n_subjects; col++){
        for(int row = col; row < n_subjects; row++){
	    try{
            phi2_value = pm2->get(ibdids[row], ibdids[col]);
           
	    }catch(...){
		phi2_value = 0;
	    }
            phi2[col*n_subjects + row] = phi2_value;
            phi2[row*n_subjects + col] = phi2_value;
            phi2_matrix(row,col) = phi2_matrix(col,row) = phi2_value;
        }
    }
    
    double * ptr_eigenvectors = new double[n_subjects*n_subjects];
    double * ptr_eigenvalues = new double[n_subjects];
    calculate_eigenvectors_and_eigenvalues(phi2, ptr_eigenvalues, ptr_eigenvectors, n_subjects);
    
    delete [] phi2;
    
    Eigenvectors = Eigen::Map<Eigen::MatrixXd>(ptr_eigenvectors, n_subjects, n_subjects);
    Eigenvalues = Eigen::Map<Eigen::VectorXd>(ptr_eigenvalues, n_subjects);
    
    delete [] ptr_eigenvectors;
    delete [] ptr_eigenvalues;
    eigen_data_is_computed = true;
}
Eigen::MatrixXd solar_mle_setup::get_eigenvectors(){
    if(!use_ped){
        Eigen::MatrixXd null_matrix;
        return null_matrix;
    }
    
    if(!eigen_data_is_computed){
        compute_eigen_data();
    }
    
    return Eigenvectors;
}
Eigen::MatrixXd solar_mle_setup::get_phi2(){
	if(!use_ped){
		Eigen::MatrixXd null_matrix;
		return null_matrix;
	}
	if(!eigen_data_is_computed){
		compute_eigen_data();
	}
	return phi2_matrix;
}
Eigen::VectorXd solar_mle_setup::get_eigenvalues(){
    if(!use_ped){
        Eigen::VectorXd null_vector;
        return null_vector;
    }
    
    if(!eigen_data_is_computed){
        compute_eigen_data();
    }
    
    return Eigenvalues;
}
/*
struct  call_default{
   solar_mle_setup * loader;
   call_default(solar_mle_setup * i) : loader(i) {}
   double operator()(int data_index) {return loader->term_map[loader->ids[loader->id_index]].at(data_index);}
};

struct  call_zscore{
   solar_mle_setup * loader;
   call_zscore(solar_mle_setup * i) : loader(i) {}
   double operator()(int data_index) {return (loader->term_map[loader->ids[loader->id_index]].at(loader->zscore_indices[data_index]) - loader->zscore_means[data_index])/loader->zscore_sd[data_index];}
};


struct  call_inorm{
   solar_mle_setup * loader;
   call_inorm(solar_mle_setup * i) : loader(i) {}
   double operator()(int data_index) {return  loader->inorm_map[loader->ids[loader->id_index]].at(loader->inorm_indices[data_index]);}
};
*/
int solar_mle_setup::parse_expression(string term_name, int term_index){
  
    Definition * use_def = Definition::Find(term_name.c_str());
    bool use_phen = Phenotypes::available (term_name.c_str());
    if(!use_phen && use_def == 0){
    	string error_message = "Parse Expression Error with term: " + term_name;
    	throw Parse_Expression_Error(error_message);
    }
    if(use_phen){
     
        expressions[term_index] = new Expression(term_name.c_str());
        try{
           expressions[term_index]->bind(_Context);
         }catch(...){
           int pheno_index = 0;
           for(int i = 0; i < phenotype_names.size(); i++, pheno_index++){
		if(!StringCmp(term_name.c_str(), phenotype_names[i].c_str(), case_ins)){
                    break;
                }
      
           }
           if(pheno_index  >= phenotype_names.size()){
              // call_default fptr(this);
       //        function<double(int)> func = fptr();
              // _Context->add(term_name.c_str(), [this](int data_index) {return this->term_map[this->ids[this->id_index]].at(data_index);}, phenotype_names.size());
              _Context->add(term_name.c_str(), call_default(this), phenotype_names.size());
              phenotype_names.push_back(term_name);
            }else{
              _Context->add(term_name.c_str(), call_default(this), pheno_index);
            }
            expressions[term_index]->bind(_Context);
         }
       return TCL_OK;

    }else if(use_def){
        try{
            expressions[term_index] = new Expression (use_def->definition());
        }catch(...){
            return TCL_ERROR;
        }
        bool ok = false;
        while(!ok){
            try{
                expressions[term_index]->bind(_Context);
                ok = true;
            }
            catch (Undefined_Name& un)
            {
               
                char * pname;
                try{
                    if (Phenotypes::available (un.name))
                    {
                       int pheno_index = 0;
          	       for(int i = 0; i < phenotype_names.size(); i++, pheno_index++){
               	           if(!StringCmp(un.name, phenotype_names[i].c_str(), case_ins)){
                             break;
                            }
                         
                        }
          	       if(pheno_index  >= phenotype_names.size()){
             		 _Context->add(un.name, call_default(this), phenotype_names.size());
             		 phenotype_names.push_back(string(un.name));
           	       }else{
             		 _Context->add(un.name, call_default(this), pheno_index);
           	       }
              
                    }else if (0 != (pname = zscor_match (un.name))){
                        if (!Phenotypes::available (pname))
                        {
                            throw Undefined_Name (pname);
                        }
                        int term_names_index = 0;
                        for(int i = 0 ; i < phenotype_names.size(); i++, term_names_index++){
                            if(!StringCmp(phenotype_names[i].c_str(), pname, case_ins)){
                                break;
                            }
                           
                        }
                        
                        if( term_names_index < phenotype_names.size()){
                            _Context->add (un.name, call_zscore(this), zscore_names.size());
                            zscore_names.push_back(string(un.name));
                            zscore_indices.push_back(term_names_index);
                            
                        }else{
                            _Context->add (un.name, call_zscore(this), zscore_names.size());
                            zscore_names.push_back(string(un.name));
                            zscore_indices.push_back(phenotype_names.size());
                            phenotype_names.push_back(string(pname));
                        }
                    }else if (0 != (pname = inor_match (un.name))){
                        if (!Phenotypes::available (pname))
                        {
                            throw Undefined_Name (pname);
                        }
                        int term_names_index = 0;
                        for(int i = 0 ; i < phenotype_names.size(); i++, term_names_index++){
                            if(!StringCmp(phenotype_names[i].c_str(), pname, case_ins)){
                                break;
                            }
                        }
                        
                        if( term_names_index < phenotype_names.size()){
                            _Context->add (un.name, call_inorm(this), inorm_names.size());
                            inorm_names.push_back(string(un.name));
                            inorm_indices.push_back(term_names_index);
                        }else{
                            _Context->add (un.name, call_inorm(this), inorm_names.size());
                            inorm_names.push_back(string(un.name));
                            inorm_indices.push_back(phenotype_names.size());
                            phenotype_names.push_back(string(pname));
                        }
                    
                    }else{
                        throw Undefined_Name(un.name);
                    }
                }catch(Undefined_Name& un_2){
                    return TCL_ERROR;
                }
            }
        }
    }
    return TCL_OK;
    
}


double solar_mle_setup::compute_inverse_normal(double  pct){
    
    double q = 1.0 - pct;
    double mean = 0.0;
    double standard_deviation = 1.0;
    int status = 0;
    int which = 2;
    double bound = 0;
    double x = 0;
    cdfnor_(&which, &pct, &q, &x, &mean, &standard_deviation, &status, &bound);
    return x;
}


vector<double> solar_mle_setup::inormalize(vector< pair<double, size_t> >  data_in){
    vector<double> data_out(data_in.size());
    sort(data_in.begin(), data_in.end(),
         [] (const pair<double, size_t> & lhs, const pair<double, size_t> & rhs) {return lhs.first < rhs.first;
         });
    double pct, z;
    const size_t count = data_in.size();
    size_t position = 1;
    size_t current_id;
    double shared_value, sum = 0.0;
    vector< pair<double,size_t> >::iterator data_iter = data_in.begin();
    while(data_iter != data_in.end()){
        shared_value = data_iter->first;
        queue<size_t> volume_ids;
        while(shared_value == data_iter->first && data_iter != data_in.end()){
            volume_ids.push(data_iter->second);
            data_iter++;
            pct = double(position++)/(count + 1);
            z = compute_inverse_normal(pct);
            sum += z;
        }
        sum /= volume_ids.size();
        while(volume_ids.size() != 0){
            current_id = volume_ids.front();
            volume_ids.pop();
            data_out[current_id] = sum;
        }
        
        sum = 0.0;
        
    }
    
    return data_out;
    
}
void solar_mle_setup::compute_inorm_and_zscore(){
    if(inorm_names.size() != 0){
        vector< pair<double, size_t> > inorm_data(ids.size());
        vector<double> row_values(inorm_names.size());
        
        for(int i = 0; i < ids.size(); i++){
            inorm_map[ids[i]] = row_values;
        }
        for(int index = 0; index < inorm_names.size(); index++){
            
            for(int i = 0 ; i < ids.size(); i++){
                
                inorm_data[i] = pair<double, size_t> (term_map.at(ids[i]).at(inorm_indices[index]), i);
            }
            vector<double> data_out = inormalize(inorm_data);
            for(int i = 0; i < ids.size();i++){
                inorm_map[ids[i]].at(index) =  data_out[i];
            }
        }
    }
    
    if(zscore_names.size() != 0){
        
        
        zscore_means.resize(zscore_names.size());
        fill(zscore_means.begin(), zscore_means.end(), 0.0);
        
        zscore_sd.resize(zscore_names.size());
        fill(zscore_sd.begin(), zscore_sd.end(), 0.0);
        for(int i = 0; i < ids.size(); i++){
            vector<double> row_values = term_map[ids[i]];
            for(int index = 0; index < zscore_names.size(); index++){
                zscore_means[index] += row_values[zscore_indices[index]];
            }
        }
        
        for(int index = 0 ; index < zscore_names.size(); index++){
            zscore_means[index] /= ids.size();
        }
        
        for(int i = 0; i < ids.size(); i++){
            vector<double> row_values = term_map[ids[i]];
            for(int index = 0; index < zscore_names.size(); index++){
                double value = row_values[zscore_indices[index]] - zscore_means[index];
                zscore_sd[index] += value*value;
            }
        }
        
        for(int index = 0 ; index < zscore_names.size(); index++){
            zscore_sd[index] = sqrt(zscore_sd[index]/ids.size());
        }
        
    }
}



int solar_mle_setup::read_phenotype_and_fill_maps(){
    const char * errmsg = 0;
   
    SolarFile * file = SolarFile::open("solar mle setup", phenotype_filename.c_str(), &errmsg);
    if(errmsg){
        return TCL_ERROR;
    }
    
    if(errmsg){
        return TCL_ERROR;
    }
    
    file->start_setup(&errmsg);
    if(errmsg){
        return TCL_ERROR;
    }
    file->setup("id", &errmsg);
    if(errmsg){
        return TCL_ERROR;
    }
    char ** file_data;
    vector<string> phenotype_ids;
    while (0 != (file_data = file->get (&errmsg))){
        
        phenotype_ids.push_back(string(file_data[0]));
        
        
    }
    file->rewind(&errmsg);
    if(errmsg){
        return TCL_ERROR;
    }
    file->start_setup(&errmsg);
    if(errmsg){
        return TCL_ERROR;
    }
    for(int i = 0 ; i < phenotype_names.size(); i++){
        file->setup(phenotype_names[i].c_str(), &errmsg);
        if(errmsg){
            string error_message = "Misc Error: " + string(errmsg);
            throw Misc_Error(error_message);
        }
    }
    vector<double> row_values(phenotype_names.size());
    vector<string>::iterator id_iter = phenotype_ids.begin();
    while (0 != (file_data = file->get (&errmsg))){
        bool skip_row = false;
        for(int index = 0; index < phenotype_names.size(); index++){
           if(!StringCmp(file_data[index], "F", case_ins)){
             row_values[index] = 1.0;
           }else  if(!StringCmp(file_data[index], "M", case_ins)){
                row_values[index] = 0.0;
           }else  if(StringCmp(file_data[index], 0, case_ins)){
                row_values[index] = strtod(file_data[index],NULL);
            }else{
                skip_row = true;
                break;
            }
        }
        
        if(!skip_row){
            term_map[*id_iter++] = row_values;
        }else{
            id_iter = phenotype_ids.erase(id_iter);
        }
        
    }
    if(use_ped){
        SolarFile * ped_file = SolarFile::open("solar mle setup", "pedindex.out", &errmsg);
        if(errmsg){
            return TCL_ERROR;
        }
        ped_file->start_setup(&errmsg);
        if(errmsg){
            return TCL_ERROR;
        }
        ped_file->setup("id", &errmsg);
        if(errmsg){
            return TCL_ERROR;
        }
        int ibdid = 1;
        while (0 != (file_data = ped_file->get (&errmsg))){
            string ped_id = string(file_data[0]);
            id_iter = find(phenotype_ids.begin(), phenotype_ids.end(), ped_id);
            if(id_iter != phenotype_ids.end()){
                ids.push_back(ped_id);
                ibdids.push_back(ibdid);
                phenotype_ids.erase(id_iter);
            }
            ibdid++;
        }
    }else{
        ids = phenotype_ids;
    }
    return TCL_OK;
}

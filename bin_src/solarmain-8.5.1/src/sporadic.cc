//
//  sporadic.cpp
//  
//
//  Created by Brian Donohue on 6/22/17.
//
//
#include "solar.h"
#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <omp.h>
#include <chrono>
#include <cctype>
using namespace std;
static const char * cortical_effects_filename = "cortical-effects.csv";
static const char * subcortical_effects_filename = "subcortical-effects.csv";
static const char * dtifa_effects_filename = "dtifa-effects.csv";
static vector<string> load_headers(const char * header_filename){
    vector<string> headers;
    ifstream headers_in(header_filename);
    if(!headers_in.is_open()) return headers;
    
    string header;
    char c;
    
    while(!headers_in.eof()){
	c = headers_in.get();
      //  if(iscntrl(c)) continue;
	if(isgraph(c) && c != ','){
	   header += c;
	}else{
	  if(header.length() != 0){
	    headers.push_back(header);
	    header.clear();
	  }
	}
    }
    if(header.length() != 0) headers.push_back(header);
    headers_in.close();
    return headers;
}

static int load_traits_per_inclusion(Tcl_Interp * interp, SolarFile * file, vector<string> traits, vector<int> lines_included, double * Y, int start, int end, int n_rows){
    const char * errmsg = 0;
    
    file->rewind(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    for(int col = start; col < end; col++){
        file->setup(traits[col].c_str(), &errmsg);
        if(errmsg){
            RESULT_LIT(errmsg);
            return TCL_ERROR;
        }
    }
    int batch_size = end - start;
    int row = 0;
    int line = 0;
    char ** file_data;
    vector<int>::iterator line_iter = lines_included.begin();
    while (0 != (file_data = file->get (&errmsg))){
        if(*line_iter == line){
            for(int col = 0; col < batch_size; col++){
               if(strlen(file_data[col]) != 0){
                   Y[col*n_rows + row] = strtod(file_data[col], NULL);
               }else{
                   Y[col*n_rows + row] = nan("");
               }
            }
            row++;
            line_iter++;
        }
        line++;
    }
        
    
    return TCL_OK;

}


static int load_traits(Tcl_Interp * interp, SolarFile * file, vector<string> traits, double * Y, int start, int end, int n_rows){
    const char * errmsg = 0;
    
  
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    const char * str_name;
    for(int col = start; col < end; col++){
        str_name = traits[col].c_str();
//        cout << str_name << endl;
        file->setup(str_name, &errmsg);
        if(errmsg){
            RESULT_LIT(errmsg);
            return TCL_ERROR;
        }
    }
    int batch_size = end - start;
    int row = 0;
    char ** file_data;
  
    while (0 != (file_data = file->get (&errmsg))){
        
        for(int col = 0; col < batch_size; col++){
            if(strlen(file_data[col]) != 0){
                Y[col*n_rows + row] = strtod(file_data[col], NULL);
            }else{
                Y[col*n_rows + row] = nan("");
            }
        }
        row++;
    }
    return TCL_OK;
}
static Eigen::MatrixXd create_covariate_matrix(vector<string> ids, vector<string> covariate_terms,
                 unordered_map<string, vector<double> > covariate_term_map, int n_covariates){
    
    
    Eigen::MatrixXd covariate_matrix = Eigen::ArrayXXd::Ones(ids.size(), n_covariates + 1);
  //  covariate_matrix(0, n_covariates) = 0.0;
    unordered_map<string, vector<double> > covariate_map;
    Eigen::MatrixXd covariate_term_matrix(ids.size(), covariate_terms.size());
    int row = 0;
    
    for(vector<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); ++row){
        
        vector<double> covariate_term_row = covariate_term_map[*id_iter++];
        
        for(int col = 0; col < covariate_terms.size(); ++col){
            covariate_term_matrix(row, col) = covariate_term_row[col];
        }
     
        
    }
    
    for(int col = 0 ; col < covariate_terms.size(); ++col){
        if(!StringCmp(covariate_terms[col].c_str(), "SEX", case_ins)){
            if((covariate_term_matrix.col(col).array() == 2.0).count() != 0){
                for(int row = 0; row < covariate_term_matrix.rows(); row++){
                    if(covariate_term_matrix(row, col) == 2.0)
                        covariate_term_matrix(row, col) = 1.0;
                    else
                        covariate_term_matrix(row, col) = 0.0;
                }
            }
        } else if(covariate_terms[col].find("snp_") != string::npos || covariate_terms[col].find("SNP_") != string::npos){
            continue;
        }else {
            covariate_term_matrix.col(col) = covariate_term_matrix.col(col).array() -covariate_term_matrix.col(col).mean();
         }
        
    }
    Covariate * cov;
    
    for(int col = 0; (cov = Covariate::index(col)); col++){
        CovariateTerm * cov_term;
        for(cov_term = cov->terms(); cov_term; cov_term = cov_term->next){
            int index = 0;
            
            for(vector<string>::iterator cov_iter = covariate_terms.begin(); cov_iter != covariate_terms.end(); cov_iter++){
                if(!StringCmp(cov_term->name, cov_iter->c_str(), case_ins)){
                    break;
                }
                index++;
            }
            
            if(cov_term->exponent == 1){
               covariate_matrix.col(col) = covariate_matrix.col(col).array()*covariate_term_matrix.col(index).array();
            }else{
                covariate_matrix.col(col) = covariate_matrix.col(col).array()*pow(covariate_term_matrix.col(index).array(), cov_term->exponent);
            }
        }
        
        

        
    }
    
   
    
    return covariate_matrix;
    
}

static int load_ids_and_covariates_per_class(Tcl_Interp * interp, const char * phenotype_filename, vector<string> covariate_terms,
                                   vector<string> & ids, unordered_map<string, vector<double> >  &  covariate_map, vector<int> & lines_included, int class_number){
    const char * errmsg = 0;
    
    SolarFile * file = SolarFile::open("sporadic", phenotype_filename, &errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
                                       
    file->setup("id", &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
                                       
    file->setup("class", &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    char ** file_data;
    int line = 0;
    
    while (0 != (file_data = file->get (&errmsg))){
        if(strtol(file_data[1], NULL, 10) == class_number && strlen(file_data[1]) != 0){
            ids.push_back(string(file_data[0]));
            lines_included.push_back(line);
        }
        line++;
    }
    file->rewind(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    
    for(vector<string>::iterator cov_iter = covariate_terms.begin(); cov_iter != covariate_terms.end(); ++cov_iter){
        file->setup(cov_iter->c_str(), &errmsg);
        
        if(errmsg){
            RESULT_LIT(errmsg);
            return TCL_ERROR;
        }
    }
    
    vector<double> covariates(covariate_terms.size());
    vector<string>::iterator id_iter = ids.begin();
    vector<int>::iterator line_iter = lines_included.begin();
    line = 0;
    while (0 != (file_data = file->get (&errmsg))){
        if(line++ == *line_iter){
            vector<double>::iterator cov_iter = covariates.begin();
            bool skip_row = false;
            for(int cov = 0; cov < covariates.size(); ++cov){
                if(!StringCmp(file_data[cov], "M", case_ins)){
                    covariates[cov] = 0.0;
                }else if(!StringCmp(file_data[cov], "F", case_ins)){
                    covariates[cov] = 1.0;
                }else if(StringCmp(file_data[cov], 0, case_ins)){
                    covariates[cov] = strtod(file_data[cov], NULL);
                }else{
                    skip_row = true;
                    break;
                }
                
            }
        
            if(!skip_row){
            
                covariate_map[*id_iter] = covariates;
            }else{
                vector<double> blank_vector;
                covariate_map[*id_iter] = blank_vector;
            }
            ++id_iter;
            ++line_iter;
        }
    }
    
    
    
    
    return TCL_OK;
}


static int load_ids_and_covariates(Tcl_Interp * interp, const char * phenotype_filename, vector<string> covariate_terms,
                                   vector<string> & ids, unordered_map<string, vector<double> >  &  covariate_map){
    const char * errmsg = 0;
   
    SolarFile * file = SolarFile::open("sporadic", phenotype_filename, &errmsg);

    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }

    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->setup("id", &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    char ** file_data;
    while (0 != (file_data = file->get (&errmsg))){
        
        ids.push_back(string(file_data[0]));
        
        
    }
    file->rewind(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    
    
    for(vector<string>::iterator cov_iter = covariate_terms.begin(); cov_iter != covariate_terms.end(); ++cov_iter){
        file->setup(cov_iter->c_str(), &errmsg);
        
        if(errmsg){
            RESULT_LIT(errmsg);
            return TCL_ERROR;
        }
    }
    
    vector<double> covariates(covariate_terms.size());
    vector<string>::iterator id_iter = ids.begin();

    while (0 != (file_data = file->get (&errmsg))){
        vector<double>::iterator cov_iter = covariates.begin();
        bool skip_row = false;
        for(int cov = 0; cov < covariates.size(); ++cov){
            if(!StringCmp(file_data[cov], "M", case_ins)){
                covariates[cov] = 0.0;
            }else if(!StringCmp(file_data[cov], "F", case_ins)){
                covariates[cov] = 1.0;
            }else if(StringCmp(file_data[cov], 0, case_ins)){
                covariates[cov] = strtod(file_data[cov], NULL);
            }else{
                skip_row = true;
                break;
            }
        }
        
        if(!skip_row){
            
            covariate_map[*id_iter] = covariates;
        }else{
            vector<double> blank_vector;
            covariate_map[*id_iter] = blank_vector;
        }
        ++id_iter;
    }
    
    
    
    
    return TCL_OK;
}



extern "C" void cdfnor_ (int* which, double* p, double* q, double* x, double* mean, double* sd, int* status, double* bound);

static inline double compute_inverse_normal(double  pct){
    
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

static void inormalize(vector< pair<double, size_t> > data_in, unordered_map<size_t, size_t> output_map, double * output_data){
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
            output_data[output_map[current_id]] = sum;
        }
        
        sum = 0.0;
        
    }
}



static int sporadic_normalize(Tcl_Interp * interp,  vector<string> ids, vector<string> covariates, unordered_map<string, vector<double> > covariate_map,
                             const char * phenotype_filename, double * array_Y, double * residuals, int n_covariates){
    const char * errmsg = 0;
    vector<string> id_list;
    for(int row = 0; row < ids.size(); row++){
        double value = array_Y[row];
        string id = ids[row];
        if((value == value && (covariate_map[id].size() != 0)) || (value == value && n_covariates == 0)) {
            id_list.push_back(id);
        }
    }
    Eigen::VectorXd residual; 
    if(n_covariates != 0){
        
    
        Eigen::MatrixXd covariate_matrix = create_covariate_matrix(id_list, covariates,
                        covariate_map,  n_covariates);
        Eigen::VectorXd Y(id_list.size());
        vector<double> covariate_vector;
        int index = 0;
        for(int row = 0 ; row < ids.size(); row++){
            double value = array_Y[row];
            string id = ids[row];
            if(value == value && covariate_map[id].size() != 0)
                Y(index++) = value;
        }
                  
        Eigen::VectorXd beta  = covariate_matrix.colPivHouseholderQr().solve(Y);
      
        
        

        residual = Y - covariate_matrix*beta;
        
    }else{
        Eigen::VectorXd Y(id_list.size());
        int index = 0;
        for(int row = 0 ; row < ids.size(); row++){
            double value = array_Y[row];
            string id = ids[row];
            if(value == value)
                Y(index++) = value;
        }
        residual = Y;
    }
    vector< pair<double, size_t> > data_in;
    unordered_map<size_t, size_t> output_map;
    vector<string>::iterator bad_id_iter;
    size_t row = 0;
    for(int index = 0 ; index < ids.size(); index++){
        if(find(id_list.begin(), id_list.end(), ids[index]) != id_list.end()){
            data_in.push_back(pair<double, size_t>(residual(row), row));
            output_map[row++] = index;
        }else{
            residuals[index] = nan("");
        }
    }
    inormalize(data_in, output_map, residuals);
    
    return TCL_OK;
}

static void print_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help sporadic_normalize");
}
vector<string> check_headers_in_phenotype(const char * phenotype_filename, vector<string>  & header_list){
	ifstream phenotype_in(phenotype_filename);
	string header_line;
	getline(phenotype_in, header_line);
        phenotype_in.close();
	vector<string> fields;
	string field;
	for(int i = 0; i < header_line.length(); i++){
		char c  = header_line[i];
               // if(iscntrl(c)) continue;
		if(c != ',' && isgraph(c)){
		   field += c;
		}else{
                   if(field.length() != 0) fields.push_back(field);
		   field.clear();
		}


	}
	if(field.length() != 0) fields.push_back(field);
	
	vector<string> fields_not_found;
	vector<string>::iterator field_iter = header_list.begin();
	while(field_iter != header_list.end()) {
   	     vector<string>::iterator find_iter = find(fields.begin(), fields.end(), *field_iter);
             if(find_iter == fields.end()){
                  fields_not_found.push_back(*field_iter);
		  field_iter = header_list.erase(field_iter);
	     }else{
		  field_iter++;
	     }
	}

     return fields_not_found;

}
extern "C" int sporadicNormalizeCmd(ClientData clientData, Tcl_Interp * interp,
                                    int argc, const char * argv[]){
    
    const char * header_filename = 0;
    
    const char * out_filename = 0;
    string class_list_str;
   // auto start = std::chrono::high_resolution_clock::now();
    for(int arg = 1 ; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--list", case_ins) || !StringCmp(argv[arg], "-list", case_ins)) && arg + 1 != argc){
            header_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--out", case_ins) || !StringCmp(argv[arg], "-out", case_ins) || !StringCmp(argv[arg], "-o", case_ins) ||\
                  !StringCmp(argv[arg], "--o", case_ins)) && arg + 1 != argc){
            out_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--class", case_ins) || !StringCmp(argv[arg], "-class", case_ins) || !StringCmp(argv[arg], "-c", case_ins) ||\
                  !StringCmp(argv[arg], "--c", case_ins)) && arg + 1 != argc){
            class_list_str = string(argv[++arg]);
        }else if(!StringCmp(argv[arg], "help", case_ins)){
            print_help(interp);
            return TCL_OK;
        }else{
            RESULT_LIT("An invalid argument was entered");
            return TCL_ERROR;
        }
    }
    string phenotype_filename =Phenotypes::filenames();
    if(!phenotype_filename.length()){
	RESULT_LIT("No phenotype file is currently loaded");
	return TCL_ERROR;
    }
    if(!out_filename){
        RESULT_LIT("No output file was specified");
        return TCL_ERROR;
    }
    vector<string> headers;
    if(!header_filename){
        if(Trait::Number_Of() == 0){
            RESULT_LIT("No trait was specified or file containing a list of traits");
            return TCL_ERROR;
        }
        
        headers.push_back(string(Trait::Name(0)));
    }else{
        headers = load_headers(header_filename);
        if(!headers.size()){
            RESULT_LIT("No traits could be read in from file specified");
            return TCL_ERROR;
        }
	vector<string> missing_headers = check_headers_in_phenotype( phenotype_filename.c_str(),  headers);
        if(headers.size() == 0){
	   RESULT_LIT("Fields listed in header file were not found in loaded phenotype.");
	   return TCL_ERROR;
	}
	if(missing_headers.size() != 0){
	   cout << "The following fields listed in the header file were not found in the loaded phenotype file: \n";
	   for(vector<string>::iterator field_iter = missing_headers.begin(); field_iter != missing_headers.end(); field_iter++){
		cout << " -" << *field_iter << endl;
	    }
	}
    }

    
    Covariate * c;
    vector<string> covariate_terms;
    vector<string> ids;
    unordered_map<string, vector<double> > covariate_map;
    int n_covariates = 0;
    for (int i = 0;( c = Covariate::index(i)); i++)
    {
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
        n_covariates++;
    }


    
    
    if(class_list_str.length() != 0){
        vector<int> class_list;
        vector<int> lines_included;
       // char * token = strdup(class_list_str.c_str());
       // token = strtok((char*)class_list_str.c_str(), &token, ",");
      //  int class_number;// = strtol(class_list_str.c_str(), &token, 10);
        //class_list.push_back(class_number);
        int string_index = 0;
        string str_number;
        while(string_index != class_list_str.length()){
            if(class_list_str[string_index] != ','){
                str_number += class_list_str[string_index];
            }else{
                class_list.push_back(atoi(str_number.c_str()));
                str_number.clear();
            }
            string_index++;
        }
        if(str_number.length() != 0){
           class_list.push_back(atoi(str_number.c_str()));
        }
        /*
        while(class_number = strtol(token, &token, 10)){
            //class_number = strtol(token, &token, 10);
            class_list.push_back(class_number);
        }*/
        ofstream output_stream(out_filename);
        output_stream << "ID,Class";
        for(vector<string>::iterator header_iter = headers.begin(); header_iter != headers.end(); header_iter++){
            output_stream << "," << *header_iter;
        }
        output_stream << "\n";
        for(vector<int>::iterator class_iter = class_list.begin(); class_iter != class_list.end(); class_iter++){
            ids.clear();
            covariate_map.clear();
            lines_included.clear();
            if(load_ids_and_covariates_per_class(interp, phenotype_filename.c_str(),  covariate_terms, \
                                       ids,  covariate_map, lines_included, *class_iter) == TCL_ERROR) return TCL_ERROR;
            double * residuals = new double[ids.size()*headers.size()];
            
            double * Y = new double[ids.size()*headers.size()];
            if(headers.size() < 1000){
                const char * errmsg = 0;
                SolarFile * file = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
                if(errmsg){
                    RESULT_LIT(errmsg);
                    return TCL_ERROR;
                }
                if(load_traits_per_inclusion(interp, file, headers, lines_included, Y, 0, headers.size(), ids.size()) == TCL_ERROR){
                    return TCL_ERROR;
                }
                for(int column = 0; column < headers.size(); column++){
                    sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                       phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
                }
                
            }else{
#pragma omp parallel
                {
                    int n_threads = omp_get_num_threads();
                        const char * errmsg = 0;
                     int thread_idx = omp_get_thread_num();
                        SolarFile * file  = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
                        int batch_size = ceil(double(headers.size())/n_threads);
                        int end = (1 +thread_idx)*batch_size;
                        if(thread_idx + 1 == n_threads)
                            end = headers.size();
                        load_traits_per_inclusion(interp, file, headers, lines_included, Y + thread_idx*batch_size*ids.size(), thread_idx*batch_size,  end, ids.size());
                   delete file;
           }
                    
                   
                    
#pragma omp parallel for
                    for(int column = 0; column < headers.size(); column++){
                        sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                           phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
                    }
                    
                
 
            }
            delete [] Y;
            for(int row = 0; row < ids.size(); row++){
                output_stream << ids[row] << "," << *class_iter;
                for(int column = 0; column < headers.size(); column++){
                    double value = residuals[column*ids.size() + row];
                    if(value == value)
                        output_stream << "," << value;
                    else
                        output_stream << ",";
                }
                output_stream << endl;
            }
            
            delete [] residuals;
 
        }
        output_stream.close();
        
    }else{
        if(load_ids_and_covariates(interp, phenotype_filename.c_str(),  covariate_terms, \
                                   ids,  covariate_map) == TCL_ERROR) return TCL_ERROR;
        
        /*
         for(int index = 0; index < covariate_terms.size(); index++){
         if(!StringCmp(covariate_terms[index].c_str(), "SEX", case_ins)){
         vector<double> covariate_row;
         for(vector<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); id_iter++){
         string id = *id_iter;
         if(covariate_map[id].size() == 0)
         continue;
         if(covariate_map[id][index] == 2.0){
         covariate_map[id][index]  = 1.0;
         }else{
         covariate_map[id][index] = 0.0;
         }
         }
         break;
         }
         }*/
        double * residuals = new double[ids.size()*headers.size()];
        
        double * Y = new double[ids.size()*headers.size()];
        if(headers.size() < 100){
            const char * errmsg = 0;
            SolarFile * file = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
            if(errmsg){
                RESULT_LIT(errmsg);
                return TCL_ERROR;
            }
            if(load_traits(interp, file, headers, Y, 0, headers.size(), ids.size()) == TCL_ERROR){
                return TCL_ERROR;
            }
            for(int column = 0; column < headers.size(); column++){
                sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                   phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
            }
            
        }else{
         //  SolarFile ** files;
#pragma omp parallel
            {
               int n_threads = omp_get_num_threads();
             
               const char * errmsg = 0;
                SolarFile * file;
//#pragma omp critical
                file  = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
              //  for(int f = 0; f < n_threads; f++){
                     int thread_idx  = omp_get_thread_num();
                   
//                    const char * errmsg = 0;
                  //  files[f]  = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
                    int batch_size = ceil(double(headers.size())/n_threads);
                    int end = (1 +thread_idx)*batch_size;
                    if(thread_idx + 1 == n_threads)
                        end = headers.size();
                    load_traits(interp, file, headers, Y + thread_idx*batch_size*ids.size(), thread_idx*batch_size,  end, ids.size());
              //  }
                 delete file;
             }   
               // delete [] files;
                
#pragma omp parallel for
                for(int column = 0; column < headers.size(); column++){
                    sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                       phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
                }
                
           // }
          //  delete [] files;
            
            
            
        }
        
        delete [] Y;
        
        ofstream file_out(out_filename);
        
        
        file_out << "ID";
        for(vector<string>::iterator trait_iter = headers.begin(); trait_iter != headers.end(); trait_iter++){
            file_out << "," << *trait_iter;
        }
        file_out << endl;
        for(int row = 0; row < ids.size(); row++){
            file_out << ids[row];
            for(int column = 0; column < headers.size(); column++){
                double value = residuals[column*ids.size() + row];
                if(value == value)
                    file_out << "," << value;
                else
                    file_out << ",";
            }
            file_out << endl;
        }
        file_out.close();
        delete [] residuals;
    }
    

    
  //  auto elapsed = std::chrono::high_resolution_clock::now() - start;
    
  //  auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);
   // cout << seconds.count() << endl;
    return TCL_OK;
    
}
static void print_help_rvi(Tcl_Interp * interp){
    Solar_Eval(interp, "help rvi");
}
static vector<string> check_headers_in_phenotype_from_rvi(const char * phenotype_filename, vector<string>  & header_list){
	ifstream phenotype_in(phenotype_filename);
	string header_line;
	getline(phenotype_in, header_line);
        phenotype_in.close();
	vector<string> fields;
	string field;
	for(int i = 0; i < header_line.length(); i++){
		char c  = header_line[i];
               // if(iscntrl(c)) continue;
		if(c != ',' && isgraph(c)){
		   field += c;
		}else{
                   if(field.length() != 0) fields.push_back(field);
		   field.clear();
		}


	}
	if(field.length() != 0) fields.push_back(field);
	
	vector<string> fields_not_found;
	vector<string>::iterator field_iter = header_list.begin();
	while(field_iter != header_list.end()) {
   	     vector<string>::iterator find_iter = find(fields.begin(), fields.end(), *field_iter);
             if(find_iter == fields.end()){
                  fields_not_found.push_back(*field_iter);
		  field_iter = header_list.erase(field_iter);
	     }else{
		  field_iter++;
	     }
	}

     return fields_not_found;

}

static const char * read_effects(vector<string> & effect_names, vector<double> & effect_values, string disorder, const char * modality){
	const char * errmsg = 0;
	const char * effects_filename = 0;
	if(!StringCmp(modality, "cortical", case_ins)){
		effects_filename = cortical_effects_filename;
        }else if(!StringCmp(modality, "subcortical", case_ins)){
		effects_filename = subcortical_effects_filename;
	}else if(!StringCmp(modality, "dtifa", case_ins)){
		effects_filename = dtifa_effects_filename;
	}else{
		return "Unknown modality was selected";
	}
	

	SolarFile * file =  SolarFile::open("rvi", effects_filename, &errmsg);
	if(errmsg){
		string error_output = "Modality file: " + string(effects_filename) + " not found";
		return error_output.c_str();
	}
	file->start_setup(&errmsg);
	if(errmsg) return errmsg;
	file->setup("roi", &errmsg);
	if(errmsg) return errmsg;
	file->setup(disorder.c_str(), &errmsg);
	if(errmsg) return errmsg;
	
	char ** file_data;
	while (0 != (file_data = file->get (&errmsg))){
		effect_names.push_back(string(file_data[0]));
		effect_values.push_back(strtod(file_data[1], NULL));
	} 	
	delete file;
	return 0;
}
extern "C" int rviCmd(ClientData clientData, Tcl_Interp * interp,
                                    int argc, const char * argv[]){
    
    const char * header_filename = 0;
    
    const char * out_filename = 0;
	string disorder;
    const char *  modality;

   // auto start = std::chrono::high_resolution_clock::now();
    for(int arg = 1 ; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--modality", case_ins) || !StringCmp(argv[arg], "-modality", case_ins)) && arg + 1 != argc){
            modality  = argv[++arg];
        }else if((!StringCmp(argv[arg], "--out", case_ins) || !StringCmp(argv[arg], "-out", case_ins) || !StringCmp(argv[arg], "-o", case_ins) ||\
                  !StringCmp(argv[arg], "--o", case_ins)) && arg + 1 != argc){
            out_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--disorder", case_ins) || !StringCmp(argv[arg], "-disorder", case_ins) || !StringCmp(argv[arg], "-d", case_ins) ||\
                  !StringCmp(argv[arg], "--d", case_ins)) && arg + 1 != argc){
            disorder = string(argv[++arg]);
        }else if(!StringCmp(argv[arg], "help", case_ins)){
            print_help_rvi(interp);
            return TCL_OK;
        }else{
            RESULT_LIT("An invalid argument was entered");
            return TCL_ERROR;
        }
    }
	vector<double> effect_values;
	vector<string> effect_names;
	const char * errmsg = 0;
	errmsg = read_effects(effect_names,effect_values, disorder, modality);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}
	
    string phenotype_filename =Phenotypes::filenames();
    if(!phenotype_filename.length()){
	RESULT_LIT("No phenotype file is currently loaded");
	return TCL_ERROR;
    }
    if(!out_filename){
        RESULT_LIT("No output file was specified");
        return TCL_ERROR;
    }
    vector<string> headers = effect_names;

	vector<string> missing_headers = check_headers_in_phenotype_from_rvi( phenotype_filename.c_str(),  headers);
        if(headers.size() == 0){
	   RESULT_LIT("Fields listed in modality file were not found in loaded phenotype.");
	   return TCL_ERROR;
		}
	if(missing_headers.size() != 0){
	   cout << "The following fields listed in the modality file were not found in the loaded phenotype file: \n";
	   for(vector<string>::iterator field_iter = missing_headers.begin(); field_iter != missing_headers.end(); field_iter++){
		cout << " -" << *field_iter << endl;
	    }
	}
    

    Covariate * c;
    vector<string> covariate_terms;
    vector<string> ids;
    unordered_map<string, vector<double> > covariate_map;
    int n_covariates = 0;
    for (int i = 0;( c = Covariate::index(i)); i++)
    {
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
        n_covariates++;
    }

          if(load_ids_and_covariates(interp, phenotype_filename.c_str(),  covariate_terms, \
                                   ids,  covariate_map) == TCL_ERROR) return TCL_ERROR;
        
        /*
         for(int index = 0; index < covariate_terms.size(); index++){
         if(!StringCmp(covariate_terms[index].c_str(), "SEX", case_ins)){
         vector<double> covariate_row;
         for(vector<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); id_iter++){
         string id = *id_iter;
         if(covariate_map[id].size() == 0)
         continue;
         if(covariate_map[id][index] == 2.0){
         covariate_map[id][index]  = 1.0;
         }else{
         covariate_map[id][index] = 0.0;
         }
         }
         break;
         }
         }*/
        double * residuals = new double[ids.size()*headers.size()];
        double * Y = new double[ids.size()*headers.size()];
        SolarFile * file = SolarFile::open("rvi", phenotype_filename.c_str(), &errmsg);
        if(errmsg){
		RESULT_LIT(errmsg);
            	return TCL_ERROR;
         }
         if(load_traits(interp, file, headers, Y, 0, headers.size(), ids.size()) == TCL_ERROR){
           	return TCL_ERROR;
          }
          for(int column = 0; column < headers.size(); column++){
              sporadic_normalize(interp,  ids, covariate_terms, covariate_map,\
			  phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
          }
    
          errmsg = 0;
	file->rewind(&errmsg);
	if(errmsg){
		RESULT_LIT(errmsg);
		return TCL_ERROR;
	}

          //SolarFile * file_two = SolarFile::open("rvi", phenotype_filename.c_str(), &errmsg);
	  file->start_setup(&errmsg);
	file->setup("dx", &errmsg);
	if(errmsg){
	    RESULT_LIT(errmsg);
	    return TCL_ERROR;
	 }
		  char ** file_data;
		  unsigned * controls = new unsigned[ids.size()];
		  unsigned i = 0;
		while (0 != (file_data = file->get (&errmsg))){
			if(!StringCmp(file_data[0], "1", case_ins)){
					controls[i] = 1;
			}else{
				controls[i] = 0;
			}
			i++;
		} 
		delete file;
		double * output_data = new double[ids.size()*headers.size()];
		for(unsigned col = 0; col < headers.size(); col++){
			double mean = 0;
			double sd = 0;
			unsigned n =  0;
			for(unsigned row = 0 ; row < ids.size(); row++){
				double residual = residuals[ids.size()*col + row];
				if(residual == residual && controls[row] == 1){
					mean += residual;
					sd += residual*residual;
					n++;
				}
			}
			mean /= n;
			sd /= (n - 1);
			sd = sqrt(sd);
			//cout << headers[col] << " " << mean << " " << sd << endl;
			for(unsigned row = 0; row < ids.size(); row++){
				double residual = residuals[ids.size()*col + row];
				if(residual == residual){
					output_data[ids.size()*col + row] = (residual - mean)/sd;
				}else{
					output_data[ids.size()*col + row] = nan("");
				}
			}
		}
		unsigned effects_index_map[headers.size()];
		double effects_mean = 0.0;
		double effects_sd = 0.0;
		unsigned n = headers.size();
		for(unsigned i = 0; i < headers.size(); i++){
			vector<string>::iterator find_iter = find(effect_names.begin(), effect_names.end(), headers[i]);
			unsigned dist = distance(effect_names.begin(), find_iter);

			effects_index_map[i] = dist;
			effects_mean += effect_values[effects_index_map[i]];
		}
		effects_mean /= n;
		for(unsigned i = 0 ; i < headers.size(); i++){
			double effect_value =  effect_values[effects_index_map[i]];
			effects_sd += pow(effect_value - effects_mean, 2);
		}
		effects_sd = sqrt(effects_sd);
		ofstream file_out(out_filename);
		file_out << "ID,RVI." << disorder << "." << modality << endl;
		for(unsigned row = 0; row < ids.size(); row++){
			string current_id = ids[row];
			bool skip_id = false;
			double local_mean = 0.0;
			double local_sd = 0.0;
			
			for(unsigned col = 0; col < headers.size(); col++){
				double value = output_data[ids.size()*col + row];
				if(value == value){
					local_mean += value;
					
				}else{
					skip_id = true;
					break;
				}
			}
			if(skip_id) continue;
			local_mean /= n;
			for(unsigned col = 0; col < headers.size(); col++){
				local_sd += pow(output_data[col*ids.size() + row] - local_mean, 2);
			}
			local_sd = sqrt(local_sd);
			double corr = 0.0;
			for(unsigned col = 0; col < headers.size(); col++){
				corr += (output_data[col*ids.size() + row] - local_mean)*(effect_values[effects_index_map[col]] - effects_mean);
			}
			corr /= (local_sd*effects_sd);

			file_out << current_id << "," << corr << endl;

		}
		file_out.close();

	
		delete [] output_data;
        delete [] Y;
        delete [] residuals;

		return TCL_OK;
    
}


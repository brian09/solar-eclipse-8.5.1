#ifndef SOLAR_MLE_SETUP_H
#define SOLAR_MLE_SETUP_H
#include <unordered_map>
#include <vector>
#include "solar.h"
#include <algorithm>
#include <queue>
#include <string>
#include <Eigen/Dense>
#include <utility>
#include <iterator>
#include <exception>
using namespace std;
class Parse_Expression_Error{
private:
 	string error_message;
public:

    Parse_Expression_Error(string n) {
 
    error_message = n;
    };
    string what(){
    	return error_message;
    }
};

class Solar_File_Error{
private:
	string error_message;
public:
    char * filename;
    Solar_File_Error(string n) {
    
    error_message = n;
    };
    string what(){
    	return error_message;
    }
};
class Expression_Eval_Error{
private:
	string error_message;
public:

    Expression_Eval_Error(string n) {
         
        error_message = n;
    };
    
    string what(){
    	return error_message;
    }    
};
class Misc_Error{
private: 
    string error_message;
public:
   
    Misc_Error(string n) {
    	
    	error_message = n;
    };
    
    string what(){
    	return error_message;
    }      
};

class solar_mle_setup{
private:
    friend struct call_default;
    friend struct call_zscore;
    friend struct call_inorm;
    friend double eval_blank_expression(int );
    friend double eval_default_expression(int );
    friend double eval_inorm_expression(int );
    friend double eval_zscore_expression(int );
    void calculate_eigenvectors_and_eigenvalues (double * , double * , double *  ,int );
    void compute_eigen_data();
    void create_output_matrix();
    double compute_inverse_normal(double  );
    int read_phenotype_and_fill_maps();
    vector<double> inormalize(vector< pair<double, size_t> >  );
    int parse_expression(string , int );
    void compute_inorm_and_zscore();
    unordered_map<string, vector<double> > term_map;
    unordered_map<string, vector<double> > inorm_map;
    
    Tcl_Interp * interp;
    vector<double> zscore_means;
    vector<double> zscore_sd;
    
    vector<int> inorm_indices;
    vector<int> zscore_indices;
    
    vector<int> ibdids;
    
    vector<string> inorm_names;
    vector<string> zscore_names;
  
    vector<string> phenotype_names;
    int id_index;
    Context * _Context;
    Expression ** expressions;
    string phenotype_filename;
    bool eigen_data_is_computed = false;
    bool output_matrix_created = false;
    bool use_ped;
    int n_expressions;
    int n_subjects;
    Eigen::MatrixXd Eigenvectors;
    Eigen::VectorXd Eigenvalues;
    Eigen::MatrixXd phi2_matrix;
    Eigen::MatrixXd output_matrix;
    vector<string> expression_names;
    vector<string> ids;
public:
    Eigen::VectorXd get_eigenvalues();
    Eigen::MatrixXd get_eigenvectors();
    Eigen::MatrixXd get_phi2();
    solar_mle_setup(vector<string> , const char * , Tcl_Interp*,  bool );
    Eigen::MatrixXd return_output_matrix();
    ~solar_mle_setup();
    vector<string> get_expression_names();
    vector<string> get_ids();
};

struct  call_default{
   solar_mle_setup * const  loader;
   call_default(solar_mle_setup * i) : loader(i) {}
   double operator()(int data_index) {return loader->term_map[loader->ids[loader->id_index]].at(data_index);}
};

struct  call_zscore{
   solar_mle_setup * const  loader;
   call_zscore(solar_mle_setup * i) : loader(i) {}
   double operator()(int data_index) {return (loader->term_map[loader->ids[loader->id_index]].at(loader->zscore_indices[data_index]) - loader->zscore_means[data_index])/loader->zscore_sd[data_index];}
};


struct  call_inorm{
   solar_mle_setup * const loader;
   call_inorm(solar_mle_setup * i) : loader(i) {}
   double operator()(int data_index) {return  loader->inorm_map[loader->ids[loader->id_index]].at(loader->inorm_indices[data_index]);}
};

#endif

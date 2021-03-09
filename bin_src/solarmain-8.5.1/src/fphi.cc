//
//  fphi.cpp
//  
//
//  Created by Brian Donohue on 8/26/18.
//
#include <chrono>
#include <stdio.h>
#include <cstdlib>
#include "Eigen/Dense"
#include "solar.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "solar_mle_setup.h"
#include "solar-trait-reader.h"
#include "RicVolumeSet.h"
using namespace std;
vector<string> read_trait_list(const char * list_filename);
extern bool loadedPed ();
extern Pedigree *currentPed;
extern "C" void cdfchi_ (int*, double*, double*, double*, double*,
                         int*, double*);
static double chicdf(double chi, double df){
    double p, q, bound;
    int status = 0;
    int which = 1;
    
    
    cdfchi_ (&which, &p, &q, &chi, &df, &status, &bound);
    
    return q/2.0;
}
static  double calculate_fphi_loglik(double variance, Eigen::VectorXd Sigma,  size_t  n_subjects){
    return -0.5*(log(variance)*n_subjects +  log(Sigma.array()).sum() + n_subjects);
}
static double calculate_pvalue(Eigen::VectorXd residual,  double  loglik){

    double null_variance = residual.squaredNorm()/residual.rows();
    
    
    
    double null_loglik = calculate_fphi_loglik(null_variance, Eigen::ArrayXd::Ones(residual.rows()).matrix(), residual.rows());
    
    
    return chicdf(2.0*(loglik - null_loglik), 1);
    
}

static Eigen::VectorXd compute_Score(double SD, Eigen::VectorXd residual,Eigen::VectorXd one_minus_lambda, Eigen::MatrixXd SX, Eigen::VectorXd omega_diagonal){

    Eigen::VectorXd residual_omega = omega_diagonal.cwiseProduct(residual);
    Eigen::VectorXd beta_score = SX.transpose()*residual_omega;

    double e2_score = pow(SD, 2.0)*0.5*(one_minus_lambda.dot(residual_omega.cwiseAbs2())- one_minus_lambda.dot(omega_diagonal));

    double SD_score = (residual.dot(residual_omega) - residual.rows())/SD;
    Eigen::VectorXd score_vector(SX.cols() + 2);
    
    for(int row = 0 ; row < SX.cols(); row++){
        score_vector(row) = beta_score(row);
    }
    score_vector(SX.cols()) = e2_score;
    score_vector(SX.cols() + 1) = SD_score;
    return score_vector;
}

static Eigen::MatrixXd compute_observed_Hessian(double SD, Eigen::VectorXd residual, Eigen::VectorXd one_minus_lambda, Eigen::MatrixXd SX,  Eigen::VectorXd Sigma_inverse){
    Eigen::MatrixXd beta_hessian = SX.transpose()*Sigma_inverse.asDiagonal()*SX;
    
    Eigen::VectorXd one_minus_lambda_squared = one_minus_lambda.cwiseProduct(one_minus_lambda);
    Eigen::VectorXd beta_var_comp_hessian = pow(SD, 2.0)*SX.transpose()*(Sigma_inverse.cwiseAbs2()).cwiseProduct(one_minus_lambda.cwiseProduct(residual));
    Eigen::VectorXd beta_SD_hessian = 2.0*SX.transpose()*residual.cwiseProduct(Sigma_inverse)/SD;
    
    Eigen::MatrixXd hessian(SX.cols() + 2, SX.cols() + 2);
    for(int row = 0 ; row < SX.cols(); row++){
        for(int col = 0 ; col < SX.cols(); col++){
            hessian(row, col) = beta_hessian(row, col);
        }
        hessian(row, SX.cols()) = beta_var_comp_hessian(row);
        hessian(SX.cols(), row) = beta_var_comp_hessian(row);
        hessian(row, SX.cols() + 1) = beta_SD_hessian(row);
        hessian(SX.cols() + 1, row) = beta_SD_hessian(row);
    }
    
    Eigen::VectorXd residual_squared = residual.cwiseProduct(residual);
    double SD_hessian = -pow(SD, -2.0)*(residual.rows() - 3.0*residual_squared.dot(Sigma_inverse));
    double SD_e2_hessian = SD*one_minus_lambda.dot( (residual.cwiseProduct(Sigma_inverse)).cwiseAbs2());
    double e2_hessian = -pow(SD, 4.0)*(0.5*one_minus_lambda_squared.dot(Sigma_inverse.cwiseAbs2()) - one_minus_lambda_squared.dot(Sigma_inverse.cwiseProduct((Sigma_inverse.cwiseProduct(residual)).cwiseAbs2())));
    hessian(SX.cols(), SX.cols()) = e2_hessian;
    hessian(SX.cols() + 1, SX.cols()) = SD_e2_hessian;
    hessian(SX.cols() , SX.cols() + 1) = SD_e2_hessian;
    hessian(SX.cols() + 1, SX.cols() + 1) = SD_hessian;

    
    
    return hessian;
}
static inline double calculate_constraint(const double x){
	return 1.0/(1.0 + exp(-x));
}

static inline double calculate_dconstraint(const double x){
	const double e_x = exp(-x);
	return e_x*pow(e_x + 1, -2);
	
}


static inline double calculate_ddconstraint(const double x){
	const double e_x = exp(x);
	return -(e_x-1.0)*e_x*pow(e_x+1,-3);
	
}

static inline double calculate_dloglik(Eigen::VectorXd lambda_minus_one,Eigen::VectorXd residual_squared, Eigen::VectorXd sigma,const double variance){
	const double part_one = variance*lambda_minus_one.dot(sigma);
	const double part_two = variance*lambda_minus_one.dot(residual_squared.cwiseProduct(sigma.cwiseAbs2()));
	return -0.5*(part_one-part_two);
}
static inline double calculate_ddloglik(Eigen::VectorXd lambda_minus_one,Eigen::VectorXd residual_squared, Eigen::VectorXd sigma, const double variance){
	Eigen::VectorXd lambda_minus_one_squared = lambda_minus_one.cwiseAbs2();
	Eigen::VectorXd sigma_squared = sigma.cwiseAbs2();
	const double part_one = variance*variance*lambda_minus_one_squared.dot(sigma_squared);
	const double part_two = 2.0*variance*variance*lambda_minus_one_squared.dot(residual_squared.cwiseProduct(sigma.cwiseProduct(sigma_squared)));

	return -0.5*(-part_one + part_two);
}
static inline double calculate_ddloglik_with_constraint(const double t,const double dloglik, const double ddloglik){

	return pow(calculate_dconstraint(t), 2)*ddloglik + calculate_ddconstraint(t)*dloglik;
}

static inline double calculate_dloglik_with_constraint(const double t, const double dloglik){

	return calculate_dconstraint(t)*dloglik;
}
static Eigen::ArrayXXd find_max_loglik_2(const int precision, Eigen::VectorXd Y, Eigen::MatrixXd aux, Eigen::MatrixXd X, double & result_loglik, double & result_variance){
    size_t n_subjects = Y.rows();
    Eigen::VectorXd theta(2);
    double parameter_t = 0;
    double h2r = 0.5; 
    theta(0) = 0.5;
    theta(1) = 0.5;
    Eigen::VectorXd Sigma = aux*theta;
    Eigen::MatrixXd Omega = Sigma.cwiseInverse().asDiagonal();
    Eigen::MatrixXd XTOX  = X.transpose()*Omega*X;
    Eigen::ArrayXXd parameters;
    if(XTOX.determinant() == 0){
    	return parameters;
    }
    Eigen::VectorXd beta = XTOX.inverse()*X.transpose()*Omega*Y;
    Eigen::VectorXd residual = Y - X*beta;
    Eigen::VectorXd residual_squared = residual.cwiseAbs2();
    double variance = residual_squared.dot(Omega.diagonal())/Y.rows();
    double loglik = calculate_fphi_loglik(variance, Sigma, n_subjects); 
    Eigen::VectorXd lambda_minus_one = (aux.col(1).array() - 1.0).matrix(); 
    Eigen::VectorXd sigma_inverse_var = (Sigma*variance).cwiseInverse();
    double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma_inverse_var, variance);
    double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma_inverse_var, variance);
    double score = calculate_dloglik_with_constraint( parameter_t,dloglik);
    double hessian = calculate_ddloglik_with_constraint(parameter_t, dloglik, ddloglik);
    double delta = -score/hessian;
    double new_h2r = 0.0;
    if(delta == delta){
    	parameter_t += delta;
    	new_h2r = calculate_constraint(parameter_t);
    }
    const double end = pow(10, -precision);
    int iter = 0;
    while( delta == delta && fabs(new_h2r - h2r) >= end && ++iter < 100){
	h2r = new_h2r;

	theta(0) = 1.0 - h2r;
	theta(1) = h2r;
	Sigma = aux*theta;
	sigma_inverse_var = (Sigma).cwiseInverse();
	Omega = sigma_inverse_var.asDiagonal();
	XTOX = X.transpose()*Omega*X;
   	 if(XTOX.determinant() == 0){
    		return parameters;
   	 }	
	beta = XTOX.inverse()*X.transpose()*Omega*Y;			
	residual = Y - X*beta;
	residual_squared = residual.cwiseAbs2();
	variance = residual_squared.dot(sigma_inverse_var)/Y.rows();
	sigma_inverse_var /= variance;
	loglik =  calculate_fphi_loglik(variance, Sigma, n_subjects); 	
	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma_inverse_var, variance);
	ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma_inverse_var, variance);
	score = calculate_dloglik_with_constraint( parameter_t,dloglik);
	hessian = calculate_ddloglik_with_constraint(parameter_t, dloglik, ddloglik);
	delta = -score/hessian;

	if(delta == delta){
		parameter_t += delta;
		new_h2r = calculate_constraint(parameter_t);
			
	}
		

    }
    if((h2r >= .9 || h2r <= 0.1) && h2r == h2r){
    	double test_h2r;
    	if(h2r >= 0.9) 
    		test_h2r = 1.0;
    	else	
    		test_h2r = 0.0;
    	Eigen::VectorXd test_theta(2);
    	test_theta(0) = 1 - test_h2r;
    	test_theta(1) = test_h2r;
    	Eigen::VectorXd test_sigma = aux*test_theta;
    	Eigen::VectorXd test_sigma_inverse = test_sigma.cwiseInverse();
    	Eigen::MatrixXd test_omega = test_sigma_inverse.asDiagonal();
    	Eigen::MatrixXd test_XTOX = X.transpose()*test_omega*X;
     	if(test_XTOX.determinant() != 0){  	
    		Eigen::VectorXd test_beta = test_XTOX.inverse()*X.transpose()*test_omega*Y;
    		Eigen::VectorXd test_residual = Y - X*test_beta;
    		double test_variance = test_residual.cwiseAbs2().dot(test_sigma_inverse)/n_subjects;
    		double test_loglik = calculate_fphi_loglik(test_variance, test_sigma, n_subjects); 
    		if(test_loglik > loglik){
    			beta = test_beta;
    			theta = test_theta;
    			h2r = test_h2r;
    			variance = test_variance;
    			loglik = test_loglik;
    		}
    	}
    }
    	

    Sigma = variance*aux*theta;
    Eigen::VectorXd omega_diagonal = Sigma.cwiseInverse();
    residual = Y - X*beta;
    Eigen::VectorXd one_minus_lambda = (1.0 - aux.col(1).array()).matrix();
    Eigen::MatrixXd Hessian = compute_observed_Hessian(sqrt(variance), residual, one_minus_lambda, X,  omega_diagonal);

    Eigen::VectorXd Score = compute_Score(sqrt(variance), residual,one_minus_lambda, X,  omega_diagonal);
    parameters.resize(beta.rows() + 3, 3);
    
    for(size_t i = 0; i < beta.rows(); i++){
        parameters(i, 0) =beta(i);
        parameters(i, 2) = Score(i);
    }
    parameters(beta.rows(), 0) = 1.0 - h2r;
    parameters(beta.rows() + 1, 0) = h2r;
    parameters(beta.rows() + 2, 0) = sqrt(variance);
    parameters(beta.rows(), 2) = Score(beta.rows());
    parameters(beta.rows() + 1, 2) = Score(beta.rows());
    parameters(beta.rows() + 2, 2) = Score(beta.rows() + 1);
    if(Hessian.determinant() != 0){
	Eigen::MatrixXd inverse_hessian =  Hessian.inverse();
        Eigen::ArrayXd errors = inverse_hessian.diagonal().cwiseSqrt();
        for(size_t i = 0 ; i < beta.rows();i++){
            parameters(i, 1) = errors(i);
        }
        parameters(beta.rows(), 1) = parameters(beta.rows() + 1, 1) = errors(beta.rows());
        parameters(beta.rows() + 2, 1) = errors(beta.rows() + 1);
    }else{
    	for(int i = 0 ; i < beta.rows(); i++){
    		parameters(i, 1) = nan("");
    	}
    	parameters(beta.rows(), 1) =parameters(beta.rows() + 1, 1) = parameters(beta.rows() + 2, 1) = nan("");
    }
    
    
    result_variance = variance;
    result_loglik = loglik;
    return parameters;
}

static Eigen::ArrayXXd find_max_loglik(const size_t precision, Eigen::VectorXd Y, Eigen::MatrixXd aux, Eigen::MatrixXd X, double & result_loglik, double & result_variance){
    size_t n_subjects = Y.rows();
    Eigen::VectorXd theta(2);
    theta(0) = 0.5;
    theta(1) = 0.5;
    Eigen::VectorXd Sigma = aux*theta;
    Eigen::MatrixXd Omega = Sigma.cwiseInverse().asDiagonal();
    Eigen::MatrixXd XTEX = X.transpose()*Omega*X;
    Eigen::VectorXd beta;
    Eigen::VectorXd residual;
    double variance = nan("");;
    double loglik = nan(""); 
    double h2r = 0.5; 
    Eigen::FullPivLU<Eigen::MatrixXd> solver(XTEX);
    if(solver.isInvertible()){
          beta = solver.inverse()*X.transpose()*Omega*Y;
	  residual = Y - X*beta;
          variance = residual.dot(Omega*residual)/residual.rows();
	  loglik = calculate_fphi_loglik(variance, Sigma, n_subjects);
    }
    for(int decimal = 1; decimal <= precision; decimal++){
        Eigen::VectorXd next_beta = beta;
        Eigen::VectorXd next_theta = theta;
        double next_h2r = h2r;
        double next_variance = variance;
        double next_loglik = loglik;
        double delta = pow(10, -decimal);

        for(int i = -6; i <= 6; i++){
            double test_h2r = h2r + delta*i;
            if(i == 0 || (test_h2r < 0.0 || test_h2r > 1.0)) continue;
            Eigen::VectorXd test_theta(2);
            test_theta(0) = 1.0 - test_h2r;
            test_theta(1) = test_h2r;
            Eigen::VectorXd test_sigma = aux*test_theta;
            Eigen::MatrixXd test_omega = test_sigma.cwiseInverse().asDiagonal();
            Eigen::VectorXd test_beta;
            double test_variance = nan("");
	    double test_loglik = nan("");
	    XTEX = X.transpose()*test_omega*X;
            Eigen::FullPivLU<Eigen::MatrixXd> test_solver(XTEX);		
	    if(test_solver.isInvertible()){
         	 test_beta = test_solver.inverse()*X.transpose()*test_omega*Y;
         	 Eigen::VectorXd test_residual = Y - X*test_beta;
         	 test_variance = test_residual.dot(test_omega*test_residual)/Y.rows();
         	 test_loglik = calculate_fphi_loglik(test_variance, test_sigma, n_subjects);
   	     }
	     

            if(test_loglik == test_loglik && (test_loglik > next_loglik || next_loglik != next_loglik)){
                next_loglik = test_loglik;
                next_h2r = test_h2r;
                next_beta = test_beta;
                next_theta = test_theta;
                next_variance = test_variance;
            }
        }
        h2r = next_h2r;
        loglik = next_loglik;
        variance = next_variance;
        beta = next_beta;
        theta = next_theta;
    }
    Sigma = variance*aux*theta;
    Eigen::VectorXd omega_diagonal = Sigma.cwiseInverse();
    residual = Y - X*beta;
    Eigen::VectorXd one_minus_lambda = (1.0 - aux.col(1).array()).matrix();
    Eigen::MatrixXd Hessian = compute_observed_Hessian(sqrt(variance), residual, one_minus_lambda, X,  omega_diagonal);
   // cout << "hessian: " << endl;
   // cout << Hessian << endl;
    
    Eigen::VectorXd Score = compute_Score(sqrt(variance), residual,one_minus_lambda, X,  omega_diagonal);
    Eigen::ArrayXXd results(beta.rows() + 3, 3);
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(Hessian);
    for(size_t i = 0; i < beta.rows(); i++){
        results(i, 0) =beta(i);
        results(i, 2) = Score(i);
    }
    results(beta.rows(), 0) = 1.0 - h2r;
    results(beta.rows() + 1, 0) = h2r;
    results(beta.rows() + 2, 0) = sqrt(variance);
    results(beta.rows(), 2) = Score(beta.rows());
    results(beta.rows() + 1, 2) = Score(beta.rows());
    results(beta.rows() + 2, 2) = Score(beta.rows() + 1);
    if(qr.isInvertible()){
	Eigen::MatrixXd inverse_hessian =  qr.inverse();
	//cout << "inverse hessian: " << endl;
	//cout << inverse_hessian << endl;
        Eigen::ArrayXd errors = inverse_hessian.diagonal().cwiseSqrt();
        for(size_t i = 0 ; i < beta.rows();i++){
            results(i, 1) = errors(i);
        }
        results(beta.rows(), 1) = results(beta.rows() + 1, 1) = errors(beta.rows());
        results(beta.rows() + 2, 1) = errors(beta.rows() + 1);
    }
    result_variance = variance;
    result_loglik = loglik;
    return results;
}

static int load_fphi_matrices(Tcl_Interp * interp, Eigen::VectorXd  & trait_vector, Eigen::MatrixXd  & covariate_matrix, Eigen::MatrixXd & auxiliary_matrix, Eigen::MatrixXd & eigenvectors,\
                       Eigen::MatrixXd & phi2_matrix, string trait_name){
    
    static vector<string> saved_ids;
    static string saved_phenotype_filename;
    static string saved_pedigree_filename;
    static Eigen::MatrixXd saved_eigenvectors;
    static Eigen::VectorXd saved_eigenvalues;
    static Eigen::MatrixXd saved_phi2;
    
    
    
    
    vector<string> cov_list;
    vector<string> unique_cov_terms;
    int success;
    Covariate * c;
    
    int n_covariates = 0;
    for (int i = 0;( c = Covariate::index(i)); i++)
    {
        char buff[512];
        c->fullname(&buff[0]);
        cov_list.push_back(string(&buff[0]));
        CovariateTerm * cov_term;
        
        for(cov_term = c->terms(); cov_term; cov_term = cov_term->next){
            bool found = false;
            
            for(vector<string>::iterator cov_iter = unique_cov_terms.begin(); cov_iter != unique_cov_terms.end(); cov_iter++){
                if(!StringCmp(cov_term->name, cov_iter->c_str(), case_ins)){
                    found = true;
                    break;
                }
            }
            if(!found){
                unique_cov_terms.push_back(string(cov_term->name));
            }
        }
        n_covariates++;
    }
    
    
    string phenotype_filename = Phenotypes::filenames();
    string pedigree_filename = currentPed->filename();
    
    vector<string> terms_to_be_used;
    terms_to_be_used.push_back(trait_name);
    for(int i = 0; i < unique_cov_terms.size(); i++){
        terms_to_be_used.push_back(unique_cov_terms[i]);
    }
    const char * errmsg = 0;
   /* SolarFile * file_reader =  SolarFile::open("Read FPHI terms", phenotype_filename.c_str(), &errmsg);
    file_reader->start_setup(&errmsg);
    for(int i = 0  ; i < unique_cov_terms.size();i++){
    	bool is_there  =   file_reader->test_name (unique_cov_terms[i].c_str(),&errmsg);
    	if(is_there == false){
    		string error_message = "Field " + unique_cov_terms[i] + " not found in phenotype file";
    		RESULT_BUF(error_message.c_str());
    		delete file_reader;
    		return TCL_ERROR;
    	}
    }
    delete file_reader;*/
    solar_mle_setup * loader;
    Eigen::MatrixXd output_matrix;
    try{
        loader = new solar_mle_setup(terms_to_be_used, phenotype_filename.c_str(),interp, true);
        if(loader->get_ids().size() == 0){
        	RESULT_LIT("No data could be read for the trait and covariates selected");
        	return TCL_ERROR;
        }
        output_matrix = loader->return_output_matrix();
        if(saved_phenotype_filename == phenotype_filename && pedigree_filename == saved_pedigree_filename && loader->get_ids() == saved_ids){
            eigenvectors = saved_eigenvectors;
            auxiliary_matrix = Eigen::ArrayXXd::Ones(loader->get_ids().size(), 2).matrix();
            auxiliary_matrix.col(1) = saved_eigenvalues;
            phi2_matrix = saved_phi2;
        }else{
            eigenvectors = loader->get_eigenvectors();
            auxiliary_matrix = Eigen::ArrayXXd::Ones(loader->get_ids().size(), 2).matrix();
            auxiliary_matrix.col(1) = loader->get_eigenvalues();
            phi2_matrix = loader->get_phi2();
        }
    }catch(Parse_Expression_Error & error){
    	
       
        RESULT_BUF(error.what().c_str());
        return TCL_ERROR;
    }catch(Solar_File_Error & error){
        RESULT_BUF(error.what().c_str());
        return TCL_ERROR;
    }catch(Expression_Eval_Error & error){
        RESULT_BUF(error.what().c_str());
        return TCL_ERROR;
    }catch(Misc_Error & error){
        
        RESULT_BUF(error.what().c_str());
        return TCL_ERROR;
    }catch(Syntax_Error &e){
    
    	RESULT_LIT("Syntax Error occurred in expression");
    	return TCL_ERROR;
    }catch(Undefined_Function &e){
    	RESULT_LIT("Undefined Function Error occurred in expression");
    	return TCL_ERROR;
    }catch(Unresolved_Name &e){
    	RESULT_LIT("Unresolved Name Error occurred in expression");
    	return TCL_ERROR;
    }catch(Undefined_Name &e){
    	RESULT_LIT("Undefined Name Error occurred in expression");
    	return TCL_ERROR;
    }catch(...){
    	RESULT_LIT("Unkown error occurred reading phenotype data");
    	return TCL_ERROR;
    }
    trait_vector = output_matrix.col(0);
   /* ofstream ids_out("fphi_id_list_2.out");
    vector<string> fphi_ids = loader->get_ids();
   // ids_out << fphi_ids[0];
    for(int i = 0; i < fphi_ids.size(); i++){
    	ids_out << fphi_ids[i] << endl;
    }
    ids_out.close();*/
    if(n_covariates == 0) {
        covariate_matrix = Eigen::ArrayXXd::Ones(trait_vector.rows(),  1);
        saved_eigenvectors = eigenvectors;
        saved_eigenvalues = auxiliary_matrix.col(1);
        saved_ids = loader->get_ids();
        saved_pedigree_filename = pedigree_filename;
        saved_phenotype_filename = phenotype_filename;
        saved_phi2 = phi2_matrix;
        delete loader;
        return TCL_OK;
    }
    covariate_matrix = Eigen::ArrayXXd::Ones(trait_vector.rows(), n_covariates + 1);
    Eigen::MatrixXd covariate_term_matrix(trait_vector.rows(), terms_to_be_used.size() - 1);
    for(int col = 0; col < unique_cov_terms.size(); col++){
        
        if(!StringCmp(unique_cov_terms[col].c_str(), "SEX", case_ins)){
            if((output_matrix.col(col + 1).array() == 2.0).count() != 0){
                for(int row = 0 ; row < covariate_term_matrix.rows(); row++){
                    if(output_matrix(row, col + 1) == 2.0){
                        covariate_term_matrix(row, col) = 1.0;
                    }else{
                        covariate_term_matrix(row, col) = 0.0;
                    }
                }
            }else{
                covariate_term_matrix.col(col) = output_matrix.col(col + 1);
            }
            
        }else if(strstr(unique_cov_terms[col].c_str(), "snp_") != NULL || strstr(unique_cov_terms[col].c_str(), "SNP_") != NULL){
            covariate_term_matrix.col(col) = output_matrix.col(col + 1);
            continue;
        } else {
           covariate_term_matrix.col(col) =  output_matrix.col(col + 1).array() - output_matrix.col(col + 1).mean();
        }
    }
    
    Covariate * cov;
    
    for(int col = 0; (cov = Covariate::index(col)); col++){
        CovariateTerm * cov_term;
        for(cov_term = cov->terms(); cov_term; cov_term = cov_term->next){
            int index = 0;
            
            for(vector<string>::iterator cov_iter = unique_cov_terms.begin(); cov_iter != unique_cov_terms.end(); cov_iter++){
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
    saved_eigenvectors = eigenvectors;
    saved_eigenvalues = auxiliary_matrix.col(1);
    saved_ids = loader->get_ids();
    saved_pedigree_filename = pedigree_filename;
    saved_phenotype_filename = phenotype_filename;
    saved_phi2 = phi2_matrix;
    delete loader;
    return TCL_OK;
    
}
static void print_fphi_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help fphi");
}
static void calculate_h2r(Eigen::VectorXd residual, Eigen::MatrixXd Z,  double & h2r , double & loglik, double & SE){

    
    Eigen::VectorXd F = pow(residual.array(), 2).matrix();
    double F_sum = F.array().mean();
    
    double score =  1.0/F_sum * (Z.col(1).array()*((F.array()/F_sum) - 1.0)).sum();
    
    if (score != score) {
        h2r = 0.0;
        SE = 0.0;
        loglik = 0.0;
        return;
    }
    if(score <= 0.0){
        h2r = 0.0;
        SE = 0.0;
        loglik = 0.0;
        return;
    }
    
    Eigen::VectorXd theta = Z.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(F);//Z.colPivHouseholderQr().solve(F);
    
    
    if(theta(0) < 0.0) theta(0) = 0.0;
    if(theta(1) < 0.0) theta(1) = 0.0;
    
    
    Eigen::VectorXd weights = Z*theta;
    
    
    
    Eigen::MatrixXd weight_matrix  = pow(weights.array(), -2).matrix().asDiagonal();

    
    Eigen::VectorXd Sigma = (Z.transpose()*weight_matrix*Z).inverse()*Z.transpose()*weight_matrix*F;
    
    
    if(Sigma(0) < 0.0) Sigma(0) = 0.0;
    if(Sigma(1) < 0.0) Sigma(1) = 0.0;
    h2r = Sigma(1)/(Sigma(1) + Sigma(0));
    weights = Z*Sigma;
    
    if(h2r != h2r) {
        h2r = 0.0;
        SE = 0.0;
        loglik = 0.0;
        return;
    }
    Eigen::VectorXd omega = Z*Sigma;
    

    loglik = -0.5*(log(omega.array()).sum() + residual.dot(residual.cwiseQuotient(omega)));
    
    
    weights = pow(omega.array(), -2).matrix();
    double a = weights.sum();
    
    double b = (Z.col(1).array()*weights.array()).sum();
    
    double c = (pow(Z.col(1).array(), 2)*weights.array()).sum();
    
    double det = a*c - pow(b, 2.0);
    
    double G = Sigma(1)/pow(Sigma(1) + Sigma(0), 2);
    double E = Sigma(0)/pow(Sigma(1) + Sigma(0), 2);
    
    double var =2.0*(pow(G, 2.0)*c + 2.0*(G*E)*b + pow(E, 2.0)*a)/det;
    SE = sqrt(var);

    
}
static vector<string> convert_volume_indices(string trait){
	vector<string> output;
	string index;
	for(unsigned i = 6; i < trait.length(); i++){
		char c = trait[i];
		if(c != '_'){
			index += c;
		}else{
			output.push_back(index);
			index.clear();
		}
	}
	output.push_back(index);
	return output;
}
static const char * run_fast_fphi_trait_list(const char * list_filename, const char * phenotype_filename, string mask_filename){
	RicVolumeSet * mask_volume = 0;
	if(mask_filename.length() != 0) mask_volume = new RicVolumeSet(mask_filename);
	RicVolumeSet * h2r_volume;
	RicVolumeSet * se_volume;
	RicVolumeSet * loglik_volume;
	if(mask_volume != 0){
		h2r_volume = new RicVolumeSet(mask_volume->nx, mask_volume->ny, mask_volume->nz, 1);
		h2r_volume->NIFTIorientation = mask_volume->NIFTIorientation;
		se_volume = new RicVolumeSet(mask_volume->nx, mask_volume->ny, mask_volume->nz, 1);
		se_volume->NIFTIorientation = mask_volume->NIFTIorientation;
		loglik_volume = new RicVolumeSet(mask_volume->nx, mask_volume->ny, mask_volume->nz, 1);
		loglik_volume->NIFTIorientation = mask_volume->NIFTIorientation;
	}
	vector<string> trait_list = read_trait_list(list_filename);
	if(trait_list.size() == 0){
		return "No traits could be read from given list file";
	}
	vector<string> empty_vector;
	Solar_Trait_Reader * reader = new Solar_Trait_Reader(phenotype_filename, trait_list, empty_vector);
	if(mask_volume == 0){
		ofstream output_stream("list-fphi.out");
		output_stream << "Trait,h2r,loglik,SE\n";
		output_stream.close();
	}
	chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
	for(unsigned set = 0; set < reader->get_n_sets(); set++){
		Eigen_Data * eigen_data = reader->get_eigen_data_set(set);
		const unsigned n_subjects = eigen_data->get_n_subjects();
		Eigen::MatrixXd eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>(eigen_data->get_eigenvectors_transposed(), n_subjects, n_subjects);
		Eigen::VectorXd eigenvalues = Eigen::Map<Eigen::VectorXd>(eigen_data->get_eigenvalues(), n_subjects);
		Eigen::MatrixXd aux_matrix = Eigen::ArrayXXd::Ones(n_subjects, 2);
		aux_matrix.col(1) = eigenvalues;
		vector<string> trait_list = eigen_data->get_trait_names();
		vector<double> h2r_list(trait_list.size());
		vector<double> loglik_list(trait_list.size());
		vector<double> SE_list(trait_list.size());
#pragma omp parallel for
		for(unsigned trait = 0; trait < trait_list.size() ; trait++){
			Eigen::VectorXd raw_Y = Eigen::Map<Eigen::VectorXd>(eigen_data->get_phenotype_column(trait), n_subjects);
			Eigen::VectorXd residual = eigenvectors_transposed*(raw_Y.array() - raw_Y.mean()).matrix();
			double h2r, loglik, SE;
			calculate_h2r(residual, aux_matrix, h2r, loglik,SE);
			h2r_list[trait] = h2r;
			loglik_list[trait] = loglik;
			SE_list[trait] = SE;
		}
		if(mask_volume == 0){
			ofstream output_stream("list-fphi.out", std::ofstream::app);
			for(unsigned trait = 0; trait < trait_list.size(); trait++){
				output_stream << trait_list[trait] << "," << h2r_list[trait] << "," << loglik_list[trait] << "," << SE_list[trait] << "\n";
			}
			output_stream.close();
		}else{
			for(unsigned i = 0; i < trait_list.size(); i++){
				vector<string> indices = convert_volume_indices(trait_list[i]);
				unsigned x = stoi(indices[0]);
				unsigned y = stoi(indices[1]);
				unsigned z = stoi(indices[2]);
				h2r_volume->VolSet[0].vox[x][y][z] = h2r_list[i];
				loglik_volume->VolSet[0].vox[x][y][z] = loglik_list[i];
				se_volume->VolSet[0].vox[x][y][z] = SE_list[i];
			}
		}

	}
          chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
     //   cout << time_span.count() << " seconds\n";
	if(mask_volume != 0){
		string base_filename = "-fphi.nii.gz";
		h2r_volume->Write("h2r"+base_filename);
		loglik_volume->Write("loglik" + base_filename);
		se_volume->Write("se"+base_filename);
		delete h2r_volume;
		delete loglik_volume;
		delete mask_volume;
		delete se_volume;
	}

	delete reader;

	return 0;	
			


}
static void format_fphi_output(double h2r, double se, double pvalue, const char * pedigree_filename, const char * phenotype_filename,\
                               const char * trait_name, size_t n_subjects){
    
    cout << "***********************************************************************************\n";
    cout << "*        Fast Permutation Heritability Inference (FPHI) Summary of Results        *\n";
    cout << "***********************************************************************************\n";
    
    cout << "  Pedigree:    " << pedigree_filename << endl;
    cout <<  "  Phenotypes:  " << phenotype_filename << endl;
    string pvalue_message = "(Significant)";
    if(pvalue >=  0.05){
        pvalue_message = "(Insignificant)";
    }
    
    cout <<   "  Trait:   " << trait_name  << " H2r = " << h2r << " SE = " << se <<  " p = " << pvalue << " " << pvalue_message << "  Individuals:  " << n_subjects << endl;
    
}
extern "C" int runfphiCmd(ClientData clientData, Tcl_Interp * interp,
                          int argc, const char * argv[]){
    bool converge = true;
    bool return_sporadic = false;
    bool debug_mode = false;
    bool web_mode = false;
    string mask_filename;
    const char * list_filename = 0;
    float sampling_rate;
    int n_subsamples;
    for(int arg = 1 ;arg < argc ; arg++){
        if(!StringCmp(argv[arg], "help", case_ins) || !StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins)
           || !StringCmp(argv[arg], "h", case_ins) || !StringCmp(argv[arg], "-h", case_ins) || !StringCmp(argv[arg], "--help", case_ins)){
            print_fphi_help(interp);
            return TCL_OK;
        }else if ((!StringCmp(argv[arg], "-fast", case_ins) || !StringCmp(argv[arg], "--f", case_ins) || !StringCmp(argv[arg], "-f", case_ins) ||
                  !StringCmp(argv[arg], "--fast", case_ins)) && arg + 2 < argc ){
            sampling_rate = atof(argv[++arg]);
            n_subsamples = atoi(argv[++arg]);
            converge = false;
        }else if (!StringCmp(argv[arg], "-d", case_ins) || !StringCmp(argv[arg], "--d", case_ins) || !StringCmp(argv[arg], "-debug", case_ins) ||
                  !StringCmp(argv[arg], "--debug", case_ins)){
            
            debug_mode  = true;
        }else if ((!StringCmp(argv[arg], "-list", case_ins) || !StringCmp(argv[arg], "--list", case_ins)) && arg + 1 < argc){
            list_filename = argv[++arg];
        }else if ((!StringCmp(argv[arg], "-mask", case_ins) || !StringCmp(argv[arg], "--mask", case_ins)) && arg + 1 < argc){
            mask_filename = string(argv[++arg]);
        } else{
            RESULT_LIT("Invalid argument enter see help");
            return TCL_ERROR;
        }
    }
    
    if(!loadedPed()){
        RESULT_LIT("No pedigree has been loaded");
        return TCL_ERROR;
    }
    const char * phenotype_filename = 0;
    phenotype_filename = Phenotypes::filenames();
    if(!phenotype_filename){
	RESULT_LIT("No phenotype file has bee loaded");
	return TCL_ERROR;
    }
    const char * pedigree_filename = currentPed->filename();
    size_t subject_count;
    Eigen::VectorXd trait_vector;
    Eigen::MatrixXd cov_matrix;
    Eigen::MatrixXd eigenvectors;
    Eigen::MatrixXd aux;
    if(list_filename){
	try{
	  	load_phi2_matrix(interp);
	}catch(...){
		RESULT_LIT("phi2 matrix could not be loaded");
		return TCL_ERROR;
	}
	const char * error_message = 0;
	error_message = run_fast_fphi_trait_list(list_filename, phenotype_filename, mask_filename);
	if(!error_message){
		RESULT_LIT(error_message);
		return TCL_ERROR;
	}

	return TCL_OK;
    }
    
    int success;
    
    if(Trait::Number_Of() == 0){
        RESULT_LIT("No trait has been selected");
        return TCL_ERROR;
    }
    
    vector<string> trait_list;
    trait_list.push_back(string(Trait::Name(0)));
    Eigen::MatrixXd phi2_matrix;
    success = load_fphi_matrices(interp, trait_vector, cov_matrix, aux, eigenvectors,phi2_matrix, trait_list[0]);
    
    if(success == TCL_ERROR) return success;
    subject_count = aux.rows();
    const char * trait_name = Trait::Name(0);
    
    
    Eigen::VectorXd trait_v = trait_vector;
    if(converge){
        double loglik = 0.0;
        double SE = 0.0;
        double pvalue = 0.0;
        double h2r = 0.0;
        double variance;
        Eigen::VectorXd Y = eigenvectors.transpose()*trait_v;
        Eigen::MatrixXd X = eigenvectors.transpose()*cov_matrix;
        Eigen::VectorXd theta(2);
     //   auto start = std::chrono::high_resolution_clock::now();
        Eigen::ArrayXXd parameters =find_max_loglik_2(8, Y, aux, X, loglik,  variance);
//	auto finish = std::chrono::high_resolution_clock::now();
 //       auto dur = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
        if(parameters.rows() == 0){
        	cout << "Convergence Failure\n";
        	return TCL_OK;
        }
        pvalue = 0.5;
        if(parameters(parameters.rows() - 2, 0) != 0.0){
            Eigen::VectorXd residual = trait_v - cov_matrix*(cov_matrix).colPivHouseholderQr().solve(trait_v);
            pvalue =  calculate_pvalue(residual, loglik);
        }
        cout.precision(9);
        format_fphi_output ( parameters(cov_matrix.cols() + 1, 0), parameters(cov_matrix.cols() + 1, 1),  pvalue,  pedigree_filename,  phenotype_filename,\
                            trait_list[0].c_str(), trait_v.rows());
          cout << "Fully Converged Parameters\n\n";
          cout << setw(10) << "Name" << setw(20) << "Value" << setw(20) << "Std Error" << endl;
        
        cout.precision(9);
        cout << endl << endl;
         Covariate * c;
         cout << setw(10) << "Parameter" << setw(20) << "Fit Value" << setw(20) << "Standard Error" << endl;
         int n_covariates = 0;
         string param_command;
         for (int i = 0;( c = Covariate::index(i)); i++)
         {
         char buff[512];
         cout << setw(10) << c->fullname(&buff[0])  << setw(20) << parameters(i, 0) << setw(20) << parameters(i, 1) << endl;
         string covar_name = "b" + string(c->fullname(&buff[0]));
         param_command = "parameter " + covar_name + " = " + to_string(parameters(i, 0));
         Solar_Eval(interp, param_command.c_str());
         param_command = "parameter " + covar_name + " se " + to_string(parameters(i, 1));
         Solar_Eval(interp, param_command.c_str());
         param_command = "parameter " + covar_name + " score " + to_string(parameters(i, 2));
         Solar_Eval(interp, param_command.c_str());
         n_covariates++;
         }
    
       
        double SD = parameters(cov_matrix.cols() + 2, 0);
        
        
        cout << setw(10) << "mean" << setw(20) << parameters(n_covariates, 0) << setw(20) << parameters(n_covariates, 1) << endl;
        param_command = "parameter mean = " + to_string(parameters(n_covariates, 0));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter mean se " + to_string(parameters(n_covariates, 1));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter mean score " + to_string(parameters(n_covariates, 2));
        Solar_Eval(interp, param_command.c_str());
        
        cout << setw(10) << "e2" << setw(20) << parameters(cov_matrix.cols(), 0) << setw(20) << parameters(cov_matrix.cols(), 1) << endl;
        param_command = "parameter e2 = " + to_string(parameters(n_covariates + 1, 0));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter e2 se " + to_string(parameters(n_covariates + 1, 1));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter e2 score " + to_string(parameters(n_covariates + 1, 2));
        Solar_Eval(interp, param_command.c_str());
        
        cout << setw(10) << "h2r" << setw(20) << parameters(cov_matrix.cols() + 1, 0)<< setw(20) << parameters(cov_matrix.cols() + 1, 1) << endl;
        param_command = "parameter h2r = " + to_string(parameters(n_covariates + 2, 0));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter h2r se " + to_string(parameters(n_covariates +2, 1));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter h2r score " + to_string(parameters(n_covariates + 2, 2));
        Solar_Eval(interp, param_command.c_str());
        
        cout << setw(10) << "sd" << setw(20) << SD << setw(20) << parameters(cov_matrix.cols()+2, 1) << endl << endl;
        param_command = "parameter SD = " + to_string(parameters(n_covariates + 3, 0));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter SD se " + to_string(parameters(n_covariates +3, 1));
        Solar_Eval(interp, param_command.c_str());
        param_command = "parameter SD score " + to_string(parameters(n_covariates + 3, 2));
        Solar_Eval(interp, param_command.c_str());
      //  cout << dur.count() << " microseconds\n";
    }else{	

        Eigen::VectorXd Y = eigenvectors.transpose()*(trait_v.array() - trait_v.mean()).matrix();
        Eigen::VectorXd residual;
        if(cov_matrix.cols() > 1) {
            Eigen::MatrixXd X = eigenvectors.transpose()*cov_matrix;
            Eigen::MatrixXd hat_matrix = Eigen::MatrixXd::Identity(X.rows(), X.rows()) - X*(X.transpose()*X).inverse()*X.transpose();
            //Eigen::VectorXd beta = X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y);
            residual = hat_matrix*Y;// - X*beta;
        }else{
            residual = Y;
        }

        double h2r ;
        double loglik;
        double SE;
        calculate_h2r( residual, aux,  h2r , loglik,  SE);
        double  pvalue = 0.5;
        if(h2r != 0.0) {
            pvalue = calculate_pvalue(residual, loglik);
        }else{
            pvalue = 0.5;
            loglik = calculate_fphi_loglik((residual.squaredNorm()/residual.rows()),Eigen::ArrayXd::Ones(residual.rows()).matrix(), residual.rows());
            
        }
        format_fphi_output ( h2r, SE,  pvalue,  pedigree_filename,  phenotype_filename,\
                            trait_list[0].c_str(), trait_v.rows());
        
    }
    
    
    
    
    return TCL_OK;
    
}

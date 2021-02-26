#include "solar.h"
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <array>
#include <iterator>
#include "plinkio.h"
#include <unordered_map>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <iomanip> 
#include "solar-trait-reader.h"
using namespace std;

void load_phi2_matrix(Tcl_Interp * interp);
static const unsigned GWAS_BATCH_SIZE = 6000;
static const int PERMUTATION_BATCH_SIZE = 1000;
extern "C" void cdfchi_ (int*, double*, double*, double*, double*,
			 int*, double*);
			 
			 
 double chicdf(double chi, double df){
	double p, q, bound;
	int status = 0;
	int which = 1;
	
	
	cdfchi_ (&which, &p, &q, &chi, &df, &status, &bound);
	
	return q;
}
double gwas_chicdf(double chi, double df){
return chicdf(chi, df);
}
extern "C" void symeig_ (int*, double*, double*, double*, double*, int*);
void calculate_eigenvectors (double * phi2, double * eigenvalues, double * eigenvectors ,int n)
{
    double* e =  new double[n];
    memset(e, 0, sizeof(double)*n);
    int * info = new int;
    *info  = 0;
    symeig_(&n, phi2, eigenvalues, e, eigenvectors, info);
    delete [] e;
    delete [] info;
}


static int create_gwas_matrices(Tcl_Interp * interp, pio_file_t * plink_file, int * & plink_index_map, Eigen::VectorXd & trait_vector, Eigen::MatrixXd & phi2,  int * & snp_values, const char * trait_name, const char * phenotype_filename, const char * plink_filename, vector<string> & snp_names){
    const char * errmsg = 0;
    cout << "Reading phenotype file...\n";
    SolarFile * file = SolarFile::open("gwas", phenotype_filename, &errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    int field_count;
    const char ** fields = file->names(&field_count, &errmsg);
    
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
    vector<string> phenotype_ids;
    
    while (0 != (file_data = file->get (&errmsg))){
        
        phenotype_ids.push_back(string(file_data[0]));
        
        
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
    file->setup(trait_name, &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    double * trait_data = new double[phenotype_ids.size()];
    double row_value;
    size_t row_index = 0;

    

    vector<string> missing_trait_ids;
    while (0 != (file_data = file->get (&errmsg))){
        string str_data = file_data[0];
        if(str_data.length() != 0){
            trait_data[row_index] = strtod((const char*)file_data[0], NULL);
        }else{
            missing_trait_ids.push_back(phenotype_ids[row_index]);
            trait_data[row_index] = nan("");
        }
        row_index++;
    }
    int * raw_plink_data;
    
    size_t n_snps;
   // size_t n_snp_rows = phenotype_ids.size();
    vector<string> plink_ids;
    if(plink_file){
      //  pio_file_t * plink_file = new pio_file_t;
       // pio_status_t status = pio_open(plink_file, plink_filename);
	pio_status_t status;
       /* 
        if(status != PIO_OK){
            pio_close(plink_file);
            if(status == P_FAM_IO_ERROR){
                RESULT_LIT("Error in loading .fam file");
                return TCL_ERROR;
            }else if (status == P_BIM_IO_ERROR){
                RESULT_LIT("Error in loading .bim file");
                return TCL_ERROR;
            }else if (status == P_BED_IO_ERROR){
                RESULT_LIT("Error in loading .bed file");
                return TCL_ERROR;
            }
        }
        */
        n_snps  = plink_file->bed_file.header.num_loci;
        size_t n_snp_rows = plink_file->bed_file.header.num_samples;
        pio_sample_t * sample;
        pio_locus_t * locus;
       // int * snp_index_map = new int[n_snp_rows];
        row_index = 0;
        for(size_t i = 0; i < n_snp_rows; i++){
            sample = fam_get_sample(&plink_file->fam_file, i);
            string id  = string(sample->iid);
            vector<string>::iterator find_iter = find(phenotype_ids.begin(), phenotype_ids.end(), id);
            if(find_iter != phenotype_ids.end()){
                vector<string>::iterator trait_iter = find(missing_trait_ids.begin(),missing_trait_ids.end(), id);
                if(trait_iter == missing_trait_ids.end()){
                    plink_ids.push_back(id);
                   // snp_index_map[i] = row_index++;
                }//else{
                 //   snp_index_map[i] = -1;
              //  }
            }//else{
                //snp_index_map[i] = -1;
           // }
        }

        for(size_t i = 0; i < n_snps; i++){
            locus = bim_get_locus(&plink_file->bim_file, i);
            snp_names.push_back(string(locus->name));
        }
      //  raw_plink_data = new int[plink_ids.size()*n_snps];
      //  snp_t * snp_buffer = new snp_t[n_snp_rows];
       // for(size_t snp = 0 ; snp < n_snps; snp++){
       //     pio_next_row(plink_file,snp_buffer);
       //     for(size_t row =0 ;row < n_snp_rows; row++){
          //      int index = snp_index_map[row];
         //       if(index != -1) raw_plink_data[snp*plink_ids.size() +index] = snp_buffer[row];
         //   }
       // }
	//delete [] snp_index_map;
     //   delete [] snp_buffer;
      //  pio_close(plink_file);
     //   delete plink_file;
    }else{
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
        for(size_t field = 0; field < field_count; field++){
            if(strstr(fields[field], "snp_")){
                string snp_name = string(fields[field]);
                snp_name = snp_name.substr(4, snp_name.length() - 4);
                snp_names.push_back(snp_name);
                file->setup(fields[field], &errmsg);
                if(errmsg){
                    RESULT_LIT(errmsg);
                    return TCL_ERROR;
                }
            }
        }
        n_snps = snp_names.size();
        size_t n_snp_rows = phenotype_ids.size();
        int * temp_plink_data = new int[n_snps*phenotype_ids.size()];
        size_t row_index = 0;
        vector<string> missing_plink_ids;
        while (0 != (file_data = file->get (&errmsg))){
            size_t missing_count = 0;
            for(size_t snp_index = 0; snp_index < n_snps; snp_index++){
                string str_data = file_data[snp_index];
		if(str_data.length() != 0){
                    int snp_value = stoi(string(file_data[snp_index]));
                    temp_plink_data[snp_index*phenotype_ids.size() + row_index] = snp_value;
                    if(snp_value == 3) missing_count++;
                }else{
                    temp_plink_data[snp_index*phenotype_ids.size() + row_index] = 3;
                    missing_count++;
                }
            }
            if(missing_count == n_snps){
                missing_plink_ids.push_back(phenotype_ids[row_index]);
            }
            row_index++;
        }
        int * plink_index_map = new int[phenotype_ids.size()];
        int map_index = 0;
        for(size_t row = 0; row < phenotype_ids.size(); row++){
            string id = phenotype_ids[row];
            vector<string>::iterator find_iter = find(missing_trait_ids.begin(), missing_trait_ids.end(), id);
            if(find_iter == missing_trait_ids.end()){
                find_iter = find(missing_plink_ids.begin(), missing_plink_ids.end(), id);
                if(find_iter == missing_plink_ids.end()){
                    plink_index_map[row] = map_index++;
                    plink_ids.push_back(id);
                }else{
                    plink_index_map[row] = -1;
                }
            }else{
                plink_index_map[row] = -1;
            }
        }
        raw_plink_data = new int[n_snps*plink_ids.size()];
        for(size_t row = 0 ;row < phenotype_ids.size(); row++){
            int index = plink_index_map[row];
            if(index != -1){
                for(size_t snp = 0; snp < n_snps ; snp++){
                    raw_plink_data[index + snp*plink_ids.size()] = temp_plink_data[row + snp*phenotype_ids.size()];
                }
            }
        }
        delete [] temp_plink_data;

    }

    vector<string> ids;
    vector<int> ibdids;
    SolarFile * pedindex_file = SolarFile::open("gwas", "pedindex.out", &errmsg);
    
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }

    pedindex_file->start_setup(&errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
    pedindex_file->setup("id", &errmsg);
    if(errmsg){
        RESULT_LIT(errmsg);
        return TCL_ERROR;
    }
	
    vector<string> shared_id_set = plink_ids;
  //  char ** file_data;
    int ibdid = 1;
    vector<string> ids_out;
    while (0 != (file_data = pedindex_file->get (&errmsg))){
        
        string id = string(file_data[0]);
        vector<string>::iterator find_iter = find(shared_id_set.begin(), shared_id_set.end(), id);
        if(find_iter != shared_id_set.end()){
            ids.push_back(id);
         //   shared_id_set.erase(find_iter);
            ibdids.push_back(ibdid);
        }
        
        ibdid++;
    }
    trait_vector.resize(ids.size());
    size_t id_index = 0;
    size_t row = 0;
    for(size_t row_index = 0; row_index < phenotype_ids.size(); row_index++){
        const double trait_value = trait_data[row_index];
        if(trait_value == trait_value){
            string trait_id = phenotype_ids[row_index];
            vector<string>::iterator find_iter = find(ids.begin(), ids.end(), trait_id);
            if(find_iter != ids.end()){
                trait_vector(distance(ids.begin(), find_iter)) = trait_value;
            }
        }
    }
    if(!plink_file){
    int * index_map = new int[plink_ids.size()];
    for(size_t i = 0; i < plink_ids.size(); i++){
        string plink_id = plink_ids[i];
        vector<string>::iterator find_iter = find(ids.begin(), ids.end(), plink_id);
        if(find_iter != ids.end()){
            index_map[i] = distance(ids.begin(), find_iter);
        }else{
            index_map[i] = -1;
        }
    }
    snp_values = new int[ids.size()*n_snps];
    for(size_t row = 0; row < plink_ids.size(); row++){
        int index = index_map[row];
        if(index != -1){
            for(size_t col = 0; col < n_snps; col++){
                snp_values[col*ids.size() + index] = raw_plink_data[col*plink_ids.size() + row];
            }
        }
    }
    delete [] raw_plink_data;
    }else{

	plink_index_map = new int[ids.size()];
	pio_sample_t * sample;
	vector<string> plink_id_list;
	for(unsigned i = 0; i < plink_file->bed_file.header.num_samples; i++){
		sample = fam_get_sample(&plink_file->fam_file, i);
		plink_id_list.push_back(string(sample->iid));
	}
	for(unsigned row = 0; row < ids.size(); row++){
		plink_index_map[row] =  distance(plink_id_list.begin(), find(plink_id_list.begin(), plink_id_list.end(), ids[row]));
	}
    }
 
    
    phi2.resize(ids.size(), ids.size());
    
    Matrix* pm2;
    
    pm2 = Matrix::find("phi2");
    if (!pm2) {
        Solar_Eval(interp, "matrix load phi2.gz phi2");
        pm2 = Matrix::find("phi2");
        if(!pm2){
            RESULT_LIT("phi2 matrix could not be loaded");
            return TCL_ERROR;
        }
    }

    for(int col_idx = 0; col_idx < ids.size(); col_idx++){
        for(int row_idx = col_idx; row_idx < ids.size(); row_idx++) {
            const double phi2_value = pm2->get(ibdids[row_idx], ibdids[col_idx]);
            phi2(row_idx, col_idx) = phi2_value;
            phi2(col_idx, row_idx) = phi2_value;
            
        }
        
    }
   
    return TCL_OK;
}

static int * get_permutation_indices(int * permutation_indices, int n_rows){
	iota (permutation_indices, permutation_indices + n_rows, 0);
	random_shuffle(permutation_indices, permutation_indices + n_rows);
//	return permutation_indices;

}
typedef struct gwas_data{
	double beta;
	double pvalue;
	double SE;
	double chi;
	double SD;
    double h2r;
    double loglik;
	gwas_data & operator = (const gwas_data & var) {
		beta = var.beta;
		pvalue = var.pvalue;
		SE = var.SE;
		chi = var.chi;
		SD= var.SD;
        h2r = var.h2r;
        loglik = var.loglik;
	}
	
}gwas_data;

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

static inline Eigen::VectorXd calculate_sigma(const double t, const double variance, Eigen::MatrixXd aux){
	Eigen::VectorXd theta(2);
	const double h2 = calculate_constraint(t);
	theta(0) = 1.0 - h2;
	theta(1) = h2;
	theta = theta*variance;
	return aux*theta;
}
static inline double calculate_quick_loglik(Eigen::VectorXd Sigma, const unsigned N){
	return -0.5*(-log(Sigma.array()).sum() + N);
}
static inline double calculate_quick_loglik_2(Eigen::VectorXd residual_squared, Eigen::VectorXd sigma){
	return -0.5*(-log(sigma.array()).sum() + residual_squared.dot(sigma));
}
static void gwas_maximize_newton_raphson_method_with_covariates(gwas_data * result, Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::MatrixXd U, gwas_data null_result, const int precision, string & status_string){


	double t = 0;
	double h2 = 0.5;
	Eigen::VectorXd theta(2);
	theta(0) = 0.5;
	theta(1) = 0.5;
	
	Eigen::VectorXd omega = (U*theta).cwiseInverse();

	Eigen::VectorXd beta(X.cols());
	Eigen::MatrixXd XTOX = X.transpose()*omega.asDiagonal()*X;
	if(XTOX.determinant() == 0){
       		result->beta = 0.0;
        	result->chi = 0.0;
        	result->pvalue = 1.0;
        	result->h2r = null_result.h2r;
        	result->loglik = null_result.loglik;
        	result->SD = null_result.SD;
        	result->SE = 0.0;
        	status_string = "Failure";
        	return;
        }		
	beta = XTOX.inverse()*X.transpose()*omega.asDiagonal()*Y;	
	Eigen::VectorXd residual = Y - X*beta;
	Eigen::VectorXd residual_squared = residual.cwiseAbs2();
	double variance = residual_squared.dot(omega)/Y.rows();
	Eigen::VectorXd sigma = omega*pow(variance, -1);
	double loglik = calculate_quick_loglik(sigma, Y.rows());
	Eigen::VectorXd lambda_minus_one = (U.col(1).array() - 1.0).matrix();
	double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
	double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
	double score = calculate_dloglik_with_constraint( t,dloglik);
	double hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
	double delta = -score/hessian;
	double new_h2 = 0;
	//cout << "iteration: 0\n";
	//cout << "delta: " << delta << endl;
	if (delta == delta){
		t += delta;
		new_h2 = calculate_constraint(t);	
		//cout << "new h2r: " << new_h2 << endl;
	}
	const double end = pow(10, -precision);
	int iter = 0;
	while( delta == delta && fabs(new_h2 - h2) >= end && ++iter < 50){
		h2 = new_h2;

		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		XTOX = X.transpose()*omega.asDiagonal()*X;
		if(XTOX.determinant() == 0){
       			result->beta = 0.0;
        		result->chi = 0.0;
        		result->pvalue = 1.0;
        		result->h2r = null_result.h2r;
        		result->loglik = null_result.loglik;
        		result->SD = null_result.SD;
        		result->SE = 0.0;
        		status_string = "Failure";
        		return;
        	}			
		beta = XTOX.inverse()*X.transpose()*omega.asDiagonal()*Y;			
		residual = Y - X*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega*pow(variance, -1);
	 	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
		ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
		score = calculate_dloglik_with_constraint( t,dloglik);
		hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
		delta = -score/hessian;
		//cout << "iteration: " << iter++ << endl;
		//cout << "delta: " << delta << endl;
		if(delta == delta){
			t += delta;
			new_h2 = calculate_constraint(t);
			//cout << "new h2r: " << new_h2 << endl;
		}
		//iter++;

	}

    loglik = calculate_quick_loglik(sigma, Y.rows());
    if((h2 >= .9 || h2 <= 0.1) && h2 == h2){
    	double test_h2r;
    	if(h2 >= 0.9) 
    		test_h2r = 1.0;
    	else	
    		test_h2r = 0.0;
    	Eigen::VectorXd test_theta(2);
    	test_theta(0) = 1 - test_h2r;
    	test_theta(1) = test_h2r;
    	Eigen::VectorXd test_sigma = U*test_theta;
    	Eigen::VectorXd test_sigma_inverse = test_sigma.cwiseInverse();
    	Eigen::MatrixXd test_omega = test_sigma_inverse.asDiagonal();
    	Eigen::MatrixXd test_XTOX = X.transpose()*test_omega*X;
     	if(test_XTOX.determinant() != 0){ 	
    		Eigen::VectorXd test_beta = test_XTOX.inverse()*X.transpose()*test_omega*Y;
    		Eigen::VectorXd test_residual = Y - X*test_beta;
    		double test_variance = test_residual.cwiseAbs2().dot(test_sigma_inverse)/Y.rows();
    		double test_loglik = calculate_quick_loglik(test_sigma_inverse/test_variance, Y.rows()); 
    		if(test_loglik > loglik){
    			beta = test_beta;
    			theta = test_theta;
    			h2 = test_h2r;
    			variance = test_variance;
    			loglik = test_loglik;
    			sigma = test_sigma_inverse/variance;
    		}
    	}
    }

    double chi = 2.0*(loglik - null_result.loglik);

    if(chi >= 0.0 && chi == chi && delta == delta){
        
        //        Eigen::MatrixXd hessian = calculate_REML_hessian(residual, SNP, Sigma, eigenvalues, variance);
        //      Eigen::VectorXd standard_errors = hessian.inverse().diagonal().cwiseSqrt();
        
       // double beta_se = beta_covariance(1,1);
       Eigen::MatrixXd SE_matrix = (X.transpose()*sigma.asDiagonal()*X).inverse();
        result->SE = sqrt(SE_matrix(beta.rows()-1,beta.rows() -1));
        result->beta = beta(beta.rows()-1);
        result->chi = chi;
        result->SD = sqrt(variance);
//        result->SE = 0.0;//beta_se;//standard_errors(0);
        result->pvalue = chicdf(result->chi , 1);
        result->h2r = h2;
        result->loglik = loglik;
        status_string = "Success";
    }else{
        result->beta = 0.0;
        result->chi = 0.0;
        result->pvalue = 1.0;
        result->h2r = null_result.h2r;
        result->loglik = null_result.loglik;
        result->SD = null_result.SD;
        result->SE = 0.0;
        if(iter == 50){
        	status_string = "Iteration Limit Reached";
        }else{
        	status_string = "Failure";
        }

    }
	

}
static int * permutated_indices;
static void gwas_maximize_newton_raphson_method(gwas_data * result, Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::MatrixXd U, gwas_data null_result, const int precision, string & status_string, const unsigned n_permutations = 0){
    	Eigen::VectorXd raw_A = X.col(0).cwiseAbs2();
    	Eigen::VectorXd raw_B = X.col(1).cwiseProduct(X.col(0));
   	Eigen::VectorXd raw_C = X.col(1).cwiseAbs2();
   	Eigen::VectorXd raw_D = Y.cwiseProduct(X.col(0));
    	Eigen::VectorXd raw_E = Y.cwiseProduct(X.col(1));

	double t = 0;
	double h2 = 0.5;
	Eigen::VectorXd theta(2);
	theta(0) = 0.5;
	theta(1) = 0.5;
	double A,B,C,D,E;
	Eigen::VectorXd omega = (U*theta).cwiseInverse();
	A = raw_A.dot(omega);
	B = raw_B.dot(omega);
	C = raw_C.dot(omega);
	D = raw_D.dot(omega);
	E = raw_E.dot(omega);
	Eigen::VectorXd beta(2);
	double beta_denom  = A*C - B*B;
    	beta(0) = C*D - B*E;
    	beta(1) = A*E - B*D;
   	beta /= beta_denom;	
	Eigen::VectorXd residual = Y - X*beta;
	Eigen::VectorXd residual_squared = residual.cwiseAbs2();
	double variance = residual_squared.dot(omega)/Y.rows();
	Eigen::VectorXd sigma = omega*pow(variance, -1);
	Eigen::VectorXd final_sigma = sigma;
	double loglik = calculate_quick_loglik(sigma, Y.rows());
	Eigen::VectorXd lambda_minus_one = (U.col(1).array() - 1.0).matrix();
	double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
	double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
	double score = calculate_dloglik_with_constraint( t,dloglik);
	double hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
	double delta = -score/hessian;
	double new_h2 = 0;
	//cout << "iteration: 0\n";
	//cout << "delta: " << delta << endl;
	if (delta == delta){
		t += delta;
		new_h2 = calculate_constraint(t);	
		//cout << "new h2r: " << new_h2 << endl;
	}
	const double end = pow(10, -precision);
	int iter = 0;
	while( delta == delta && fabs(new_h2 - h2) >= end && ++iter < 50){
		h2 = new_h2;
		final_sigma = sigma;
		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
		C = raw_C.dot(omega);
		D = raw_D.dot(omega);
		E = raw_E.dot(omega);
		beta_denom  = A*C - B*B;
    		beta(0) = C*D - B*E;
    		beta(1) = A*E - B*D;
   		beta /= beta_denom;	
		residual = Y - X*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega*pow(variance, -1);
	 	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
		ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
		score = calculate_dloglik_with_constraint( t,dloglik);
		hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
		delta = -score/hessian;
		//cout << "iteration: " << iter++ << endl;
		//cout << "delta: " << delta << endl;
		if(delta == delta){
			t += delta;
			new_h2 = calculate_constraint(t);
			//cout << "new h2r: " << new_h2 << endl;
		}
		//iter++;

	}

    loglik = calculate_quick_loglik(sigma, Y.rows());
      if((h2 >= .9 || h2 <= 0.1) && h2 == h2){
    	double test_h2r;
    	if(h2 >= 0.9) 
    		test_h2r = 1.0;
    	else	
    		test_h2r = 0.0;
    	Eigen::VectorXd test_theta(2);
    	test_theta(0) = 1 - test_h2r;
    	test_theta(1) = test_h2r;
    	Eigen::VectorXd test_omega =  (U*test_theta).cwiseInverse();
    	double test_A = raw_A.dot(test_omega);
    	double test_B = raw_B.dot(test_omega);
    	double test_C = raw_C.dot(test_omega);
    	double test_D = raw_D.dot(test_omega);
    	double test_E = raw_E.dot(test_omega);
    	Eigen::VectorXd test_beta(2);
	double test_beta_denom  = test_A*test_C - test_B*test_B;
    	test_beta(0) = test_C*test_D - test_B*test_E;
    	test_beta(1) = test_A*test_E - test_B*test_D;
   	test_beta /= test_beta_denom;
    	Eigen::VectorXd test_residual = Y - X*test_beta;  
    	Eigen::VectorXd test_residual_squared =test_residual.cwiseAbs2();
    	double test_variance = test_residual_squared.dot(test_omega)/Y.rows();
    	Eigen::VectorXd test_sigma = test_omega/test_variance;
    	double test_loglik = calculate_quick_loglik(test_sigma, Y.rows());
    	if(test_loglik > loglik){
    		h2 = test_h2r;
    		loglik = test_loglik;
    		variance = test_variance;
    		beta = test_beta;
    		A = test_A;
    		beta_denom = test_beta_denom;
    		final_sigma = test_sigma;
    	}
    }  

    double chi = 2.0*(loglik - null_result.loglik);
  //  Eigen::MatrixXd beta_covariance_matrix = X.transpose()*omega.asDiagonal()*X;
   // double denom = beta_covariance_matrix(0,0)*beta_covariance_matrix(1, 1) - beta_covariance_matrix(1,0)*beta_covariance_matrix(0,1);
    
    //  gwas_data result;
	//cout << "chi : " << chi << " delta: " << delta << " h2: " << h2 << endl;
    if(chi >= 0.0 && chi == chi && delta == delta){
        
        //        Eigen::MatrixXd hessian = calculate_REML_hessian(residual, SNP, Sigma, eigenvalues, variance);
        //      Eigen::VectorXd standard_errors = hessian.inverse().diagonal().cwiseSqrt();
        
       // double beta_se = beta_covariance(1,1);
        result->SE = sqrt(variance*A/beta_denom);
        result->beta = beta(1);
        result->chi = chi;
        result->SD = sqrt(variance);
//        result->SE = 0.0;//beta_se;//standard_errors(0);
        result->pvalue = chicdf(result->chi , 1);
        result->h2r = h2;
        result->loglik = loglik;
        status_string = "Success";
        if(n_permutations != 0){
        	//vector<double> permutated_sigma(final_sigma.rows());
        	//for(int i = 0 ; i < omega.rows(); i++){
        	//	permutated_sigma[i] = final_sigma(i);
        	//}
//        	Eigen::VectorXd mean_column = X.col(0);
//        	Eigen::VectorXd snp_column = X.col(1);
//        	Eigen::VectorXd demeaned_Y = Y - mean_column*(mean_column.transpose()*mean_column).inverse()*mean_column.transpose()*Y;
//        	const double compare_value = pow(demeaned_Y.dot(snp_column.cwiseProduct(final_sigma)), 2)/snp_column.cwiseAbs2().dot(final_sigma);
		const double compare_value = pow(result->beta/result->SE, 2);
        	unsigned pass_count = 1;
        	//Eigen::VectorXd sigma_Y = final_sigma.cwiseProduct(demeaned_Y);
        	Eigen::VectorXd sigma_Y = final_sigma.cwiseProduct(Y);
        	for(int i = 0; i < n_permutations ; i++){
        		//random_shuffle(permutated_sigma.begin(), permutated_sigma.end());
        		Eigen::VectorXd current_sigma(final_sigma.rows());
        		Eigen::VectorXd current_sigma_Y(final_sigma.rows());

        		for(int row = 0; row < current_sigma.rows(); row++){
        			
        			current_sigma(row) = final_sigma(permutated_indices[i*current_sigma.rows() + row]);
        			current_sigma_Y(row) = sigma_Y(permutated_indices[i*current_sigma.rows() + row]);
	
        		}

        		A = raw_A.dot(current_sigma);
			B = raw_B.dot(current_sigma);
			C = raw_C.dot(current_sigma);
			D = current_sigma_Y.dot(X.col(0));
			E = current_sigma_Y.dot(X.col(1));	
			beta_denom  = A*C - B*B;
		
			const double test_value = pow((A*E - B*D)/beta_denom, 2)/(A/beta_denom);        		
        		//const double test_value = pow(current_sigma_Y.dot(snp_column), 2)/snp_column.cwiseAbs2().dot(current_sigma);
        		if(test_value >= compare_value) ++pass_count;
        	}
        	result->pvalue = (double)pass_count/(1.0+n_permutations);
        }
        
    }else{
        result->beta = 0.0;
        result->chi = 0.0;
        result->pvalue = 1.0;
        result->h2r = null_result.h2r;
        result->loglik = null_result.loglik;
        result->SD = null_result.SD;
        result->SE = 0.0;
        if(iter == 50){
        	status_string = "Iteration Limit Reached";
        }else{
        	status_string = "Failure";
        }

    }
	

}

double calculate_GWAS_loglik(Eigen::VectorXd SY2, Eigen::VectorXd lambda, Eigen::VectorXd Sigma){
    // double var_sum  = 0.0;
    // double log_sum = 0.0;
    // Eigen::VectorXd theta(2);
    // theta(0) = 1.0-h2r;
    // theta(1) = h2r;
    // Eigen::VectorXd Sigma  = (lambda*h2r).array() + e2;
    // const double variance_sum  = SY2.cwiseQuotient(variance);
    
    return -0.5*(log(Sigma.array()).sum() + SY2.rows());
    // for(int i = 0 ; i < SY.rows(); i++)
    
}
gwas_data gwas_maximize_newton_raphson_method_with_covariates_null_model(Eigen::VectorXd Y, Eigen::MatrixXd covariate_matrix, Eigen::MatrixXd U, const int precision){
	gwas_data result;
	double t = 0;
	double h2 = 0.5;
	Eigen::VectorXd theta(2);
	theta(0) = 0.5;
	theta(1) = 0.5;
	Eigen::VectorXd omega = (U*theta).cwiseInverse();
	Eigen::VectorXd beta(covariate_matrix.cols());
	Eigen::MatrixXd XTOX = covariate_matrix.transpose()*omega.asDiagonal()*covariate_matrix;
	if(XTOX.determinant() == 0){
      		result.beta = 0.0;
        	result.chi = 0.0;
        	result.pvalue = 1.0;
        	result.h2r = 0.0;
        	result.loglik = 0.0;
        	result.SD = 0.0;
        	result.SE = 0.0;
        	return result;
        }		
	beta = XTOX.inverse()*covariate_matrix.transpose()*omega.asDiagonal()*Y;
	Eigen::VectorXd residual = Y - covariate_matrix*beta;
	Eigen::VectorXd residual_squared = residual.cwiseAbs2();
	double variance = residual_squared.dot(omega)/Y.rows();
	Eigen::VectorXd sigma = omega*pow(variance, -1);
	double loglik = calculate_quick_loglik(sigma, Y.rows());
	Eigen::VectorXd lambda_minus_one = (U.col(1).array() - 1.0).matrix();
	double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
	
	double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
	double score = calculate_dloglik_with_constraint( t,dloglik);
	double hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
	double delta = -score/hessian;
	double new_h2;// = calculate_constraint(t);
	//cout << "iteration: 0\n";
	//cout << "delta: " << delta << endl;
	if (delta == delta){
		t += delta;
		new_h2 = calculate_constraint(t);
		//cout << "new h2r: " << new_h2 << endl;
	}
	const double end = pow(10, -precision);
	int iter = 1;
	while( delta == delta && fabs(new_h2 - h2) >= end){
		h2 = new_h2;

		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		XTOX = covariate_matrix.transpose()*omega.asDiagonal()*covariate_matrix;
		if(XTOX.determinant() == 0){
      			result.beta = 0.0;
        		result.chi = 0.0;
        		result.pvalue = 1.0;
        		result.h2r = 0.0;
        		result.loglik = 0.0;
        		result.SD = 0.0;
        		result.SE = 0.0;
        		return result;
        	}		
		beta = XTOX.inverse()*covariate_matrix.transpose()*omega.asDiagonal()*Y;
		residual = Y - covariate_matrix*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega*pow(variance, -1);
	 	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
		ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
		score = calculate_dloglik_with_constraint( t,dloglik);
		hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
		delta = -score/hessian;
		//cout << "iteration: " << iter++ << endl;
		//cout << "delta: " << delta << endl;
		if(delta == delta){
			t += delta;
			new_h2 = calculate_constraint(t);
			//cout << "new h2r: " << new_h2 << endl;
		}

	}	
 loglik = calculate_quick_loglik(sigma, Y.rows());
    if((h2 >= .9 || h2 <= 0.1) && h2 == h2){
    	double test_h2r;
    	if(h2 >= 0.9) 
    		test_h2r = 1.0;
    	else	
    		test_h2r = 0.0;
    	Eigen::VectorXd test_theta(2);
    	test_theta(0) = 1 - test_h2r;
    	test_theta(1) = test_h2r;
    	Eigen::VectorXd test_sigma = U*test_theta;
    	Eigen::VectorXd test_sigma_inverse = test_sigma.cwiseInverse();
    	Eigen::MatrixXd test_omega = test_sigma_inverse.asDiagonal();
    	Eigen::MatrixXd test_XTOX = covariate_matrix.transpose()*test_omega*covariate_matrix;
     	if(test_XTOX.determinant() != 0){ 	
    		Eigen::VectorXd test_beta = test_XTOX.inverse()*covariate_matrix.transpose()*test_omega*Y;
    		Eigen::VectorXd test_residual = Y - covariate_matrix*test_beta;
    		double test_variance = test_residual.cwiseAbs2().dot(test_sigma_inverse)/Y.rows();
    		double test_loglik = calculate_quick_loglik(test_sigma_inverse/test_variance, Y.rows()); 
    		if(test_loglik > loglik){
    			beta = test_beta;
    			theta = test_theta;
    			h2 = test_h2r;
    			variance = test_variance;
    			loglik = test_loglik;
    			sigma = test_sigma_inverse/variance;
    		}
    	}
    } 
    if(delta == delta){
        result.beta = 0.0;
        result.pvalue = 1.0;
        result.SE = 0.0;
        result.loglik = loglik;
        result.chi = 0.0;
        result.h2r = h2;
        result.SD = sqrt(variance);

    }else{
        result.beta = 0.0;
        result.chi = 0.0;
        result.pvalue = 1.0;
        result.h2r = 0.0;
        result.loglik = 0.0;
        result.SD = 0.0;
        result.SE = 0.0;

    }
    return result;
}
 gwas_data gwas_maximize_newton_raphson_method_null_model(Eigen::VectorXd Y, Eigen::VectorXd mean_column, Eigen::MatrixXd U, const int precision){
	//cout << "null calculation\n";
	gwas_data result;
	Eigen::VectorXd raw_A = Y.cwiseProduct(mean_column);
	Eigen::VectorXd raw_B = mean_column.cwiseAbs2();

	double t = 0;
	double h2 = 0.5;
	Eigen::VectorXd theta(2);
	theta(0) = 0.5;
	theta(1) = 0.5;
	double A,B;
	Eigen::VectorXd omega = (U*theta).cwiseInverse();
	A = raw_A.dot(omega);
	B = raw_B.dot(omega);
	Eigen::VectorXd beta(1);

    	beta(0) = A/B;	
	Eigen::VectorXd residual = Y - mean_column*beta;
	Eigen::VectorXd residual_squared = residual.cwiseAbs2();
	double variance = residual_squared.dot(omega)/Y.rows();
	Eigen::VectorXd sigma = omega*pow(variance, -1);
	double loglik = calculate_quick_loglik(sigma, Y.rows());
	Eigen::VectorXd lambda_minus_one = (U.col(1).array() - 1.0).matrix();
	double dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
	
	double ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
	double score = calculate_dloglik_with_constraint( t,dloglik);
	double hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
	double delta = -score/hessian;
	double new_h2;// = calculate_constraint(t);
	//cout << "iteration: 0\n";
	//cout << "delta: " << delta << endl;
	if (delta == delta){
		t += delta;
		new_h2 = calculate_constraint(t);
		//cout << "new h2r: " << new_h2 << endl;
	}
	const double end = pow(10, -precision);
	int iter = 1;
	while( delta == delta && fabs(new_h2 - h2) >= end){
		h2 = new_h2;

		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
		//Eigen::VectorXd beta(2);
    		beta(0) = A/B;
		residual = Y - mean_column*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega*pow(variance, -1);
	 	dloglik = calculate_dloglik(lambda_minus_one, residual_squared, sigma, variance);
		ddloglik = calculate_ddloglik(lambda_minus_one, residual_squared, sigma, variance);
		score = calculate_dloglik_with_constraint( t,dloglik);
		hessian = calculate_ddloglik_with_constraint(t, dloglik, ddloglik);
		delta = -score/hessian;
		//cout << "iteration: " << iter++ << endl;
		//cout << "delta: " << delta << endl;
		if(delta == delta){
			t += delta;
			new_h2 = calculate_constraint(t);
			//cout << "new h2r: " << new_h2 << endl;
		}

	}
	/*if(delta == delta){
		h2 = new_h2;
		theta(0) = 1.0 - h2;
		theta(1) = h2;
		omega = (U*theta).cwiseInverse();
		A = raw_A.dot(omega);
		B = raw_B.dot(omega);
    		beta(0) = A/B;
		residual = Y - mean_column*beta;
		residual_squared = residual.cwiseAbs2();
		variance = residual_squared.dot(omega)/Y.rows();
		sigma = omega/variance;
	}*/

    loglik = calculate_quick_loglik(sigma, Y.rows());
     if((h2 >= .9 || h2 <= 0.1) && h2 == h2){
    	double test_h2r;
    	if(h2 >= 0.9) 
    		test_h2r = 1.0;
    	else	
    		test_h2r = 0.0;
    	Eigen::VectorXd test_theta(2);
    	test_theta(0) = 1 - test_h2r;
    	test_theta(1) = test_h2r;
    	Eigen::VectorXd test_omega =  (U*test_theta).cwiseInverse();
    	double test_A = raw_A.dot(test_omega);
    	double test_B = raw_B.dot(test_omega);
    	Eigen::VectorXd test_beta(1);
    	test_beta(0) = test_A/test_B;
    	Eigen::VectorXd test_residual = Y - mean_column*test_beta;  
    	Eigen::VectorXd test_residual_squared =test_residual.cwiseAbs2();
    	double test_variance = test_residual_squared.dot(test_omega)/Y.rows();
    	Eigen::VectorXd test_sigma = test_omega/test_variance;
    	double test_loglik = calculate_quick_loglik(test_sigma, Y.rows());
    	if(test_loglik > loglik){
    		h2 = test_h2r;
    		loglik = test_loglik;
    		variance = test_variance;
    	}
    }

    if(delta == delta){
        result.beta = 0.0;
        result.pvalue = 1.0;
        result.SE = 0.0;
        result.loglik = loglik;
        result.chi = 0.0;
        result.h2r = h2;
        result.SD = sqrt(variance);

    }else{
        result.beta = 0.0;
        result.chi = 0.0;
        result.pvalue = 1.0;
        result.h2r = 0.0;
        result.loglik = 0.0;
        result.SD = 0.0;
        result.SE = 0.0;

    }
    return result;
}
/*
double calculate_GWAS_loglik(Eigen::VectorXd SY2, Eigen::VectorXd lambda, Eigen::VectorXd Sigma){
    // double var_sum  = 0.0;
    // double log_sum = 0.0;
    // Eigen::VectorXd theta(2);
    // theta(0) = 1.0-h2r;
    // theta(1) = h2r;
    // Eigen::VectorXd Sigma  = (lambda*h2r).array() + e2;
    // const double variance_sum  = SY2.cwiseQuotient(variance);
    
    return -0.5*(log(Sigma.array()).sum() + SY2.rows());
    // for(int i = 0 ; i < SY.rows(); i++)
    
}*/

gwas_data  compute_null_model_MLE(Eigen::VectorXd Y, Eigen::VectorXd mean_column, Eigen::MatrixXd U, const unsigned precision){
    gwas_data result;
    double SD = 0.0;
    double loglik = 0.0;
    Eigen::VectorXd theta(2);
    theta(0) = 1.0;
    theta(1) = 0.0;
    Eigen::VectorXd Sigma = U*theta;
    Eigen::VectorXd Omega = Sigma.cwiseInverse();
    Eigen::VectorXd Y_mean_column = Y.cwiseProduct(mean_column);
    Eigen::VectorXd mean_column_squared = mean_column.cwiseAbs2();
    double mean = Omega.dot(Y_mean_column)/Omega.dot(mean_column_squared);
    
    Eigen::VectorXd residual = Y  - mean_column*mean;
    Eigen::VectorXd residual2 = residual.cwiseAbs2();
    double variance = residual2.dot(Omega)/residual2.rows();
    theta = theta*variance;
    Sigma = Sigma*variance;
    Eigen::VectorXd eigenvalues = U.col(1);
    double max_loglik = calculate_GWAS_loglik(residual2, eigenvalues, Sigma);
    double h2r = 0.0;
    for(double decimal = 0.1; decimal <= 1.0; decimal += 0.1){
        Eigen::VectorXd test_theta(2);
        test_theta(0) = 1.0 - decimal;
        test_theta(1) = decimal;
        Eigen::VectorXd test_sigma = U*test_theta;
        Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
        double test_mean = test_omega.dot(Y_mean_column)/test_omega.dot(mean_column_squared);
        
        Eigen::VectorXd test_residual = Y - mean_column*mean;
        Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
        double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
        test_sigma = test_sigma*test_variance;
        test_theta = test_theta*test_variance;
        double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
        if(test_loglik > max_loglik && test_loglik == test_loglik){
            max_loglik = test_loglik;
            theta = test_theta;
            Sigma = test_sigma;
            variance = test_variance;
            residual = test_residual;
            h2r = decimal;
            mean = test_mean;
        }
    }
    
    for(int decimal_place = 2; decimal_place <= precision; decimal_place++){
        double scale = pow(10, -decimal_place);
        Eigen::VectorXd next_theta = theta;
        Eigen::VectorXd next_sigma = Sigma;
        Eigen::VectorXd next_residual = residual;
        double next_mean = mean;
        double next_variance = variance;
        double next_h2r= h2r;
        double next_loglik = max_loglik;
        for(int i = -6; i <= 6; i++){
            double value = i*scale;
            double test_h2r = h2r + value;
            if(test_h2r > 1.0 || test_h2r < 0.0 || test_h2r == h2r){
                continue;
            }
            Eigen::VectorXd test_theta(2);
            test_theta(0) = 1.0 - test_h2r;
            test_theta(1) = test_h2r;
            Eigen::VectorXd test_sigma = U*test_theta;
            Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
            double test_mean = test_omega.dot(Y_mean_column)/test_omega.dot(mean_column_squared);
            
            Eigen::VectorXd test_residual = Y - mean_column*test_mean;
            Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
            double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
            test_sigma = test_sigma*test_variance;
            test_theta = test_theta*test_variance;
            double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
            if(test_loglik > next_loglik && test_loglik == test_loglik){
                next_loglik = test_loglik;
                next_theta = test_theta;
                next_sigma = test_sigma;
                next_variance = test_variance;
                next_residual = test_residual;
                next_mean = test_mean;
                next_h2r = test_h2r;
            }
            
            
        }
        
        max_loglik = next_loglik;
        theta = next_theta;
        Sigma = next_sigma;
        variance = next_variance;
        residual = next_residual;
        mean = next_mean;
        
        h2r = next_h2r;
        
    }
    //    cout << "null h2r " << h2r << " null loglik  " << max_loglik << endl;
    SD = 0.0;
    if(max_loglik == max_loglik){
        loglik = max_loglik;
        SD = sqrt(variance);
        result.beta = 0.0;
        result.pvalue = 1.0;
        result.SE = 0.0;
        result.loglik = loglik;
        result.chi = 0.0;
        result.h2r = h2r;
        result.SD = sqrt(variance);

    }else{
        result.beta = 0.0;
        result.pvalue = 0.0;
        result.chi = 0.0;
        result.h2r = 0.0;
        result.SD = 0.0;
        result.SE = 0.0;
        result.loglik = 0.0;
    }
    return result;
    
}

//double calculate_GWAS_loglik(Eigen::VectorXd SY2, Eigen::VectorXd lambda, Eigen::VectorXd Sigma);
static void  MLE_GWAS(gwas_data * result, Eigen::VectorXd Y, Eigen::MatrixXd X, Eigen::MatrixXd U, gwas_data null_result, const unsigned precision){
    
    Eigen::VectorXd raw_A = X.col(0).cwiseAbs2();
    Eigen::VectorXd raw_B = X.col(1).cwiseProduct(X.col(0));
    Eigen::VectorXd raw_C = X.col(1).cwiseAbs2();
    Eigen::VectorXd raw_D = Y.cwiseProduct(X.col(0));
    Eigen::VectorXd raw_E = Y.cwiseProduct(X.col(1));

    Eigen::VectorXd theta(2);
    theta(0) = 1.0;
    theta(1) = 0.0;
    
    Eigen::VectorXd Sigma = U*theta;
    Eigen::VectorXd Omega = Sigma.cwiseInverse();
    double A  = raw_A.dot(Omega);
    double B = raw_B.dot(Omega);
    double C = raw_C.dot(Omega);
    double D = raw_D.dot(Omega);
    double E = raw_E.dot(Omega);
    Eigen::VectorXd beta(2);
    double beta_denom  = A*C - B*B;
    beta(0) = C*D - B*E;
    beta(1) = A*E - B*D;
    beta /= beta_denom;
    Eigen::VectorXd residual = Y - X*beta;
    Eigen::VectorXd residual2 = residual.cwiseAbs2();
    double variance = residual2.dot(Omega)/residual2.rows();
    theta = theta*variance;
    Sigma = Sigma*variance;
    Eigen::VectorXd eigenvalues = U.col(1);
    double max_loglik = calculate_GWAS_loglik(residual2, eigenvalues, Sigma);
    double h2r = 0.0;
    for(double decimal = 0.1; decimal <= 1.0; decimal += 0.1){
        Eigen::VectorXd test_theta(2);
        test_theta(0) = 1.0 - decimal;
        test_theta(1) = decimal;
        Eigen::VectorXd test_sigma = U*test_theta;
        Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
        double test_A = raw_A.dot(test_omega);
        double test_B = raw_B.dot(test_omega);
        double test_C = raw_C.dot(test_omega);
        double test_D = raw_D.dot(test_omega);
        double test_E = raw_E.dot(test_omega);
        double test_beta_denom = test_A*test_C - test_B*test_B;
        Eigen::VectorXd test_beta(2);
        test_beta(0) = (test_C*test_D - test_B*test_E)/test_beta_denom;
        test_beta(1) = (test_A*test_E - test_B*test_D)/test_beta_denom;
        Eigen::VectorXd test_residual = Y - X*test_beta;
        Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
        double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
        test_sigma = test_sigma*test_variance;
        test_theta = test_theta*test_variance;
        double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
        if(test_loglik > max_loglik && test_loglik == test_loglik){
            max_loglik = test_loglik;
            theta = test_theta;
            Sigma = test_sigma;
            variance = test_variance;
            residual = test_residual;
            h2r = decimal;
            beta = test_beta;
        }
    }
    
    for(int decimal_place = 2; decimal_place <= precision; decimal_place++){
        double scale = pow(10, -decimal_place);
        Eigen::VectorXd next_theta = theta;
        Eigen::VectorXd next_sigma = Sigma;
        Eigen::VectorXd next_residual = residual;
        Eigen::VectorXd next_beta = beta;
        double next_variance = variance;
        double next_h2r= h2r;
        double next_loglik = max_loglik;
        for(int i = -5; i <= 6; i++){
            double value = i*scale;
            double test_h2r = h2r + value;
            if(test_h2r > 1.0 || test_h2r < 0.0 || test_h2r == h2r){
                continue;
            }
            Eigen::VectorXd test_theta(2);
            test_theta(0) = 1.0 - test_h2r;
            test_theta(1) = test_h2r;
            Eigen::VectorXd test_sigma = U*test_theta;
            Eigen::VectorXd test_omega = test_sigma.cwiseInverse();
            double test_A = raw_A.dot(test_omega);
            double test_B = raw_B.dot(test_omega);
            double test_C = raw_C.dot(test_omega);
            double test_D = raw_D.dot(test_omega);
            double test_E = raw_E.dot(test_omega);
            double test_beta_denom = test_A*test_C - test_B*test_B;
            Eigen::VectorXd test_beta(2);
            test_beta(0) = (test_C*test_D - test_B*test_E)/test_beta_denom;
            test_beta(1) = (test_A*test_E - test_B*test_D)/test_beta_denom;
            Eigen::VectorXd test_residual = Y - X*test_beta;
            Eigen::VectorXd test_residual2 = test_residual.cwiseAbs2();
            double test_variance = test_residual2.dot(test_omega)/test_residual2.rows();
            test_sigma = test_sigma*test_variance;
            test_theta = test_theta*test_variance;
            double test_loglik = calculate_GWAS_loglik(test_residual2, eigenvalues, test_sigma);
            if(test_loglik > next_loglik && test_loglik == test_loglik){
                next_loglik = test_loglik;
                next_theta = test_theta;
                next_sigma = test_sigma;
                next_variance = test_variance;
                next_residual = test_residual;
                next_h2r = test_h2r;
                next_beta = test_beta;
            }
            
            
        }
        
        max_loglik = next_loglik;
        theta = next_theta;
        Sigma = next_sigma;
        variance = next_variance;
        residual = next_residual;
        h2r = next_h2r;
        beta = next_beta;
        
    }
    double chi = 2.0*(max_loglik - null_result.loglik);
    Eigen::MatrixXd beta_covariance_matrix = X.transpose()*Sigma.cwiseInverse().asDiagonal()*X;
    double denom = beta_covariance_matrix(0,0)*beta_covariance_matrix(1, 1) - beta_covariance_matrix(1,0)*beta_covariance_matrix(0,1);
    
    //  gwas_data result;
    if(chi >= 0.0 && chi == chi && denom > 0.0 && denom == denom){
        
        //        Eigen::MatrixXd hessian = calculate_REML_hessian(residual, SNP, Sigma, eigenvalues, variance);
        //      Eigen::VectorXd standard_errors = hessian.inverse().diagonal().cwiseSqrt();
        
       // double beta_se = beta_covariance(1,1);
        result->SE = sqrt(beta_covariance_matrix(0, 0)/denom);
        result->beta = beta(1);
        result->chi = chi;
        result->SD = sqrt(variance);
//        result->SE = 0.0;//beta_se;//standard_errors(0);
        result->pvalue = chicdf(result->chi , 1);
        result->h2r = h2r;
        result->loglik = max_loglik;
    }else{
        result->beta = 0.0;
        result->chi = 0.0;
        result->pvalue = 1.0;
        result->h2r = null_result.h2r;
        result->loglik = null_result.loglik;
        result->SD = null_result.SD;
        result->SE = 0.0;

    }
    
    
    // return result;
    
    
    
}

static void calculate_pvalue(vector<double>::iterator & pvalue_iterator, Eigen::MatrixXd syP_Sigma_P, Eigen::MatrixXd sigmaP, Eigen::MatrixXd snps_matrix){
	
	Eigen::ArrayXXd Ts_numerator = (snps_matrix*syP_Sigma_P).array();	
	
	Eigen::ArrayXXd Ts_denominator = ((snps_matrix.array()*snps_matrix.array()).matrix()*sigmaP).array();
	
    
	Eigen::ArrayXXd Ts = Ts_numerator*Ts_numerator/Ts_denominator;
	//double compare_value;
	/*double value;
	for(int row = 0 ; row < n_snps;row++){
	pvalue = 0.0;
	compare_value = Ts(0, 0);
	pvalue = (Ts >= compare_value).count()/(double)sigmaP.cols();	
	*/
//#pragma omp parallel for
	for(int row = 0; row < Ts.rows(); row++){
		*pvalue_iterator++ = (Ts.row(row) >= Ts(row, 0)).count()/double(Ts.cols());
	}
}	
static void print_gwas_help(Tcl_Interp * interp){
	Solar_Eval(interp, "help gwas");
}
static void calculate_eigen_data_null(Eigen::VectorXd & Y,  Eigen::MatrixXd & eigenvectors_transposed,  Eigen::MatrixXd & U,  Eigen::MatrixXd  phi2){
    const int n_subjects = Y.rows();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(phi2);
    
    Eigen::MatrixXd temp_eigenvectors = es.eigenvectors();
    Eigen::VectorXd temp_eigenvalues = es.eigenvalues();
    const int n_subjects_reml = n_subjects - (temp_eigenvalues.array() == 0.0).count();
    gwas_data result;

        eigenvectors_transposed = temp_eigenvectors.transpose();
        U.resize(n_subjects, 2);
        U.col(0) = Eigen::ArrayXd::Ones(n_subjects).matrix();
        U.col(1) = temp_eigenvalues;
        Y = eigenvectors_transposed * Y;//(Y.array() - Y.mean()).matrix();

    
}

static vector<gwas_data> GWAS_MLE_fix_missing_run_with_covariates(gwas_data null_result, Eigen::VectorXd default_Y, Eigen::MatrixXd default_covariate_matrix,\
						Eigen::MatrixXd snp_matrix,Eigen::MatrixXd default_U, \
                                               const int n_snps,  const unsigned precision, vector<string> & status_vector){
    // const int n_subjects = Y.rows();

    vector<gwas_data> results(n_snps);
#pragma omp parallel for schedule(dynamic)
    for(int iteration = 0; iteration < n_snps; iteration++){
	//cout << "snp: " << iteration << endl;
	Eigen::MatrixXd local_covariate_matrix = default_covariate_matrix;
	local_covariate_matrix.col(local_covariate_matrix.cols() - 1) = snp_matrix.col(iteration);
	

       /* Eigen::MatrixXd local_X(n_subjects , 2);
        local_X.col(0) = default_mean;
        int local_n_subjects = n_subjects;
        double SNP_mean = 0.0;
        for(int row = 0 ; row < n_subjects; row++){
            local_X(row, 1) = SNP_data[iteration*n_subjects + row];
            if(local_X(row, 1) != 3.0){
                SNP_mean += local_X(row, 1);
            }else{
                local_n_subjects--;
            }
        }
        SNP_mean /= local_n_subjects;
        for(int row = 0; row < n_subjects; row++){
            if(local_X(row, 1) != 3.0){
                local_X(row, 1) -= SNP_mean;
            }else{
                local_X(row, 1) = 0.0;
            }
        }*/
	//cout << "multiplying evectors times snp\n";
       // local_X.col(1) = default_eigenvectors_transposed*local_X.col(1);
        gwas_data result;
	//cout << "starting newton raphson\n";
gwas_maximize_newton_raphson_method_with_covariates(&result, default_Y, local_covariate_matrix, default_U, null_result, precision, status_vector[iteration]);
        //MLE_GWAS(&result,default_Y, local_X, default_U, null_result, precision);
       // if(result.chi == 0.0) result.SD = default_SD;
        results[iteration] = result;
        
    }
    
    return results;
}
static vector<gwas_data> GWAS_MLE_fix_missing_run(gwas_data null_result, Eigen::VectorXd default_Y, Eigen::VectorXd default_mean,\
						Eigen::MatrixXd snp_matrix,Eigen::MatrixXd default_U, \
                                               const int n_snps,  const unsigned precision, vector<string> & status_vector, const unsigned n_permutations = 0){
    // const int n_subjects = Y.rows();

    vector<gwas_data> results(n_snps);
#pragma omp parallel for schedule(dynamic)
    for(int iteration = 0; iteration < n_snps; iteration++){
	//cout << "snp: " << iteration << endl;
	Eigen::MatrixXd local_X(default_Y.rows(), 2);
	local_X.col(0) = default_mean;
	local_X.col(1) = snp_matrix.col(iteration);
       /* Eigen::MatrixXd local_X(n_subjects , 2);
        local_X.col(0) = default_mean;
        int local_n_subjects = n_subjects;
        double SNP_mean = 0.0;
        for(int row = 0 ; row < n_subjects; row++){
            local_X(row, 1) = SNP_data[iteration*n_subjects + row];
            if(local_X(row, 1) != 3.0){
                SNP_mean += local_X(row, 1);
            }else{
                local_n_subjects--;
            }
        }
        SNP_mean /= local_n_subjects;
        for(int row = 0; row < n_subjects; row++){
            if(local_X(row, 1) != 3.0){
                local_X(row, 1) -= SNP_mean;
            }else{
                local_X(row, 1) = 0.0;
            }
        }*/
	//cout << "multiplying evectors times snp\n";
       // local_X.col(1) = default_eigenvectors_transposed*local_X.col(1);
        gwas_data result;
	//cout << "starting newton raphson\n";
gwas_maximize_newton_raphson_method(&result, default_Y, local_X, default_U, null_result, precision, status_vector[iteration], n_permutations);
        //MLE_GWAS(&result,default_Y, local_X, default_U, null_result, precision);
       // if(result.chi == 0.0) result.SD = default_SD;
        results[iteration] = result;
        
    }
    
    return results;
}


static vector<gwas_data> GWAS_MLE_run(gwas_data default_null_result,Eigen::VectorXd Y, Eigen::VectorXd default_Y, Eigen::VectorXd default_mean,\
				Eigen::MatrixXd default_U, Eigen::MatrixXd default_eigenvectors_transposed, Eigen::MatrixXd default_phi2,\
				int * SNP_data, const int n_snps, const int n_subjects, const unsigned precision, vector<string> & status_vector){
    //const int n_subjects = Y.rows();

    vector<gwas_data> results(n_snps);
#pragma omp parallel for schedule(dynamic)
    for(int iteration = 0; iteration < n_snps; iteration++){
        Eigen::VectorXd local_Y;
        Eigen::MatrixXd local_U;
        Eigen::MatrixXd local_X;

        Eigen::MatrixXd local_eigenvectors_transposed;
        gwas_data local_null_result = default_null_result;
        double local_SD;
        Eigen::MatrixXd local_phi2;
        Eigen::VectorXd local_SNP(n_subjects);
        int local_n_subjects = n_subjects;
        for(int row = 0 ; row < n_subjects; row++){
            local_SNP(row) = SNP_data[iteration*n_subjects + row];
            if(local_SNP(row) == 3.0) local_n_subjects--;
        }
        if(local_n_subjects == 0){
            gwas_data result;
            result.beta = 0.0;
            result.SD = 0.0;
            result.SE = 0.0;
            result.chi = 0.0;
            result.pvalue = 0.0;
            results[iteration] = result;
        }else{
        if(local_n_subjects == n_subjects){
            local_U = default_U;
            local_eigenvectors_transposed = default_eigenvectors_transposed;
            local_Y = default_Y;
            local_SNP = local_eigenvectors_transposed*local_SNP;
            Eigen::MatrixXd temp_X(n_subjects, 2);
            temp_X.col(0) = default_mean;
            temp_X.col(1) = local_SNP;
            local_X = temp_X;
            local_SD = default_null_result.SD;
	    local_null_result = default_null_result;
        }else{
            Eigen::VectorXd temp_local_Y(local_n_subjects);
            Eigen::MatrixXd local_phi2(local_n_subjects,local_n_subjects);
            local_X = Eigen::ArrayXXd::Ones(local_n_subjects, 2).matrix();
            
            int local_row = 0;
            for(int row = 0; row < n_subjects; row++){
                if(local_SNP(row) != 3.0){
                    temp_local_Y(local_row) = Y(row);
                    local_X(local_row, 1) = local_SNP(row);
                    int local_col = local_row;
                    for(int col = row; col < n_subjects; col++){
                        if( local_SNP(col) != 3.0){
                            local_phi2(local_row, local_col) = default_phi2(row,col);
                            local_phi2(local_col++, local_row) = default_phi2(col, row);
                        }
                    }
                    local_row++;
                }
            }
            local_Y = temp_local_Y;
            calculate_eigen_data_null(local_Y,  local_eigenvectors_transposed, local_U,   local_phi2);
            
            local_X = local_eigenvectors_transposed*local_X;
            Eigen::VectorXd local_mean = local_X.col(0);

            local_null_result = gwas_maximize_newton_raphson_method_null_model(local_Y, local_mean, local_U, precision);
            
            
            
        }
        gwas_data result;
        if(local_null_result.loglik == local_null_result.loglik){
            gwas_maximize_newton_raphson_method(&result, local_Y, local_X, local_U, local_null_result, precision, status_vector[iteration]);
        }else{
            result.chi  = 0.0;
            result.beta = 0.0;
            result.loglik = 0.0;
            result.SD = 0.0;
            result.SE = 0.0;
            result.h2r = 0.0;
            result.pvalue = 0.0;
        }
        results[iteration] = result;
        }
    }
    
    return results;
    
}
 Eigen::VectorXd calculate_theta(Eigen::VectorXd Y, Eigen::MatrixXd aux){
    Eigen::VectorXd theta(2);
    
    theta(0) = 0.0;
    theta(1) = 0.0;
    
    Eigen::VectorXd F = Y.cwiseAbs2();
    
    const double F_mean = F.mean();
    
    const double score = aux.col(1).dot(((F/F_mean).array() - 1.0).matrix())/F_mean;
    
    if(score <= 0.0) return theta;
    
    theta = aux.colPivHouseholderQr().solve(F);
    if(theta(0) <= 0.0 && theta(1) > 0.0){
    	 theta(0) = 0.0;
    	 Eigen::VectorXd single_value = aux.col(1).colPivHouseholderQr().solve(F);
    	 if(single_value(0) > 0.0){
    	 	theta(1) = single_value(0);
    	 }else{
    	 	theta(1) = 0.0;
    	 }
    }else if (theta(0) > 0.0 && theta(1) <= 0.0){
    	theta(1) = 0.0;
    	Eigen::VectorXd single_value = aux.col(0).colPivHouseholderQr().solve(F);
    	if(single_value(0) > 0.0){
    		theta(0) = single_value(0);
    	}else{
    		theta(0) = 0.0;
    	}
    }   
    
    if(theta(0) < 0.0) theta(0) = 0.0;
    
    if(theta(1) < 0.0) theta(1) = 0.0;
    
    if(theta(0) == 0.0 && theta(1) == 0.0) return theta;
    
    Eigen::VectorXd Sigma = aux*theta;
    
    Eigen::MatrixXd Omega = Sigma.cwiseAbs2().cwiseInverse().asDiagonal();
    
    theta = (Omega*aux).colPivHouseholderQr().solve(Omega*F);
    
    if(theta(0) <= 0.0 && theta(1) > 0.0){
    	 theta(0) = 0.0;
    	 Eigen::VectorXd single_value = (Omega*aux.col(1)).colPivHouseholderQr().solve(Omega*F);
    	 if(single_value(0) > 0.0){
    	 	theta(1) = single_value(0);
    	 }else{
    	 	theta(1) = 0.0;
    	 }
    }else if (theta(0) > 0.0 && theta(1) <= 0.0){
    	theta(1) = 0.0;
    	Eigen::VectorXd single_value = (Omega*aux.col(0)).colPivHouseholderQr().solve(Omega*F);
    	if(single_value(0) > 0.0){
    		theta(0) = single_value(0);
    	}else{
    		theta(0) = 0.0;
    	}
    }else if (theta(0) <= 0.0 && theta(1) <= 0.0){
    	theta(0) = 0.0;
    	theta(1) = 0.0;
    }  	

    
    
    return theta;
    
    
}
static void calculate_pvalues_permutation_method(vector<double>::iterator & pvalue_iterator,Eigen::MatrixXd Sigma_Y_permutated, Eigen::MatrixXd Sigma_permutated, Eigen::MatrixXd SNP_data){
    Eigen::ArrayXXd Ts_numerator = (SNP_data*Sigma_Y_permutated).cwiseAbs2().array();
    
    Eigen::ArrayXXd Ts_denominator = (SNP_data.cwiseAbs2()*Sigma_permutated).array();
    Eigen::ArrayXXd Ts = Ts_numerator/Ts_denominator;
   
    for(int row = 0; row < Ts.rows(); row++){
        *pvalue_iterator++ = (Ts.row(row) >= Ts(row, 0)).count()/double(Ts.cols());
    }
    
   
    
}


static void compute_pvalue_batch(vector<double>::iterator & pvalue_iterator, Eigen::MatrixXd Sigma_Y_permutated, Eigen::MatrixXd Sigma_permutated, Eigen::MatrixXd eigenvectors_transposed,\
                                 int * SNP_data, const int batch_size,const int n_subjects, const int n_permutations){
    Eigen::MatrixXd SNP_matrix(n_subjects, batch_size);
    
    for(int col = 0 ; col < batch_size; col++){
        double SNP_mean = 0.0;
        int subject_count = 0;
        for(int subject = 0; subject < n_subjects; subject++){
            int value = SNP_data[col*n_subjects + subject];
            if(value != 3){
                SNP_mean+= value;
                subject_count++;
            }
            SNP_matrix(subject, col) = value;
        }
        SNP_mean /= subject_count;
        for(int subject = 0; subject < n_subjects; subject++){
            if(SNP_matrix(subject,col) != 3.0){
                SNP_matrix(subject,col) -= SNP_mean;
            }else{
                SNP_matrix(subject, col) = 0.0;
            }
        }
    }
    SNP_matrix = eigenvectors_transposed*SNP_matrix;
    SNP_matrix.transposeInPlace();
    calculate_pvalues_permutation_method(pvalue_iterator, Sigma_Y_permutated, Sigma_permutated, SNP_matrix);
    
    
}

static void compute_permutation_data(Eigen::VectorXd Sigma_Y, Eigen::VectorXd Sigma, Eigen::MatrixXd & Sigma_Y_permutated, \
                                     Eigen::MatrixXd & Sigma_permutated, const int n_subjects, const int n_permutations){
    vector<int> indices(n_subjects);
    iota(indices.begin(), indices.end(), 0);
    
    Sigma_Y_permutated.col(0) = Sigma_Y;
    Sigma_permutated.col(0) = Sigma;
    
    
    for(int col = 1 ; col < n_permutations + 1; col++){
        
        random_shuffle(indices.begin(), indices.end());
        
        for(int row = 0; row < n_subjects; row++){
            Sigma_Y_permutated(row, col) = Sigma_Y(indices[row]);
            Sigma_permutated(row, col) = Sigma(indices[row]);
        }
        
    }
    
}

static vector<double> compute_pvalues_permutation(Eigen::VectorXd Y, Eigen::VectorXd Sigma, Eigen::MatrixXd eigenvectors_transposed, int * snp_data, \
                                                    const int n_snps, const int n_subjects, const  int n_permutations){
    
    Eigen::VectorXd Sigma_Y = Y.cwiseProduct(Sigma);
    
    Eigen::MatrixXd Sigma_Y_permutated(n_subjects, n_permutations + 1);
    Eigen::MatrixXd Sigma_permutated(n_subjects, n_permutations + 1);
    
    compute_permutation_data(Sigma_Y, Sigma, Sigma_Y_permutated, \
                             Sigma_permutated,  n_subjects,  n_permutations);
    vector<double> pvalues(n_snps);
    if(n_snps <= PERMUTATION_BATCH_SIZE){
        vector<double>::iterator pvalue_iterator = pvalues.begin();
        compute_pvalue_batch(pvalue_iterator, Sigma_Y_permutated, Sigma_permutated, eigenvectors_transposed, \
                             snp_data, n_snps, n_subjects, n_permutations);
    }else{
        const int batch_size = PERMUTATION_BATCH_SIZE;
        const int n_batches = ceil(double(n_snps)/batch_size);
        int remainder_snps = n_snps % batch_size;
        if(remainder_snps == 0)
            remainder_snps = batch_size;

#pragma omp parallel for
        for(int batch = 0; batch < n_batches; batch++){
            int  thread_batch_size = batch_size;
            if(batch == n_batches - 1)
                thread_batch_size = remainder_snps;
             vector<double>::iterator pvalue_iterator = pvalues.begin() + batch_size*batch;
            compute_pvalue_batch(pvalue_iterator, Sigma_Y_permutated, Sigma_permutated, eigenvectors_transposed, \
                                 snp_data + n_subjects*batch*batch_size, thread_batch_size, n_subjects, n_permutations);
        }
    }
    

    
    return pvalues;
    
}


vector<string> read_trait_list(const char * list_filename){
	ifstream input_stream(list_filename);
	vector<string> output;
	if(!input_stream.is_open()) return output;
	string trait;
	while(input_stream >> trait){
		output.push_back(trait);
	}

	return output;

}

static const char * load_covariate_terms(const char * phenotype_filename, Eigen::MatrixXd & covariate_term_matrix, vector<string> covariate_terms, vector<string> & ids){
	const char * errmsg = 0;
	SolarFile * file_reader =  SolarFile::open("Read GWAS Covariate Terms", phenotype_filename, &errmsg);
	if(errmsg) return errmsg;
	file_reader->start_setup(&errmsg);
	if(errmsg) return errmsg;
	file_reader->setup("id", &errmsg);
	if(errmsg) return errmsg;
	for(int i = 0; i < covariate_terms.size(); i++){
		file_reader->setup(covariate_terms[i].c_str(), &errmsg);
		if(errmsg) return errmsg;
	}
	
	vector< vector<double> > file_values;
	vector<string> temp_ids;
	char ** file_data;
    while (0 != (file_data = file_reader->get (&errmsg))){
        
        string id = file_data[0];
        
        bool skip_row = false;
        vector<double> row_values(covariate_terms.size());
        for(int i = 1; i < covariate_terms.size() + 1; i++){
        	if(StringCmp(file_data[i], 0, case_ins)){
        		if(!StringCmp(file_data[i], "m", case_ins)){
        			row_values[i - 1] = 0;
        		}else if (!StringCmp(file_data[i], "f", case_ins)){
        			row_values[i - 1] = 1;
        		}else{
        			row_values[i - 1] = atof (file_data[i]);
        		}
        	}else{
        		skip_row = true;
        		break;
        	}
        }
        
        if(!skip_row){
        	temp_ids.push_back(id);
        	file_values.push_back(row_values);
        }
    }
    
    delete file_reader;
    
    covariate_term_matrix.resize(temp_ids.size(), covariate_terms.size());
    for(int row = 0; row < covariate_term_matrix.rows(); row++){
    	vector<double> row_values = file_values[row];
    	for(int col = 0; col < covariate_term_matrix.cols(); col++){
    		covariate_term_matrix(row, col) = row_values[col];
    	}
    }
    for(int i = 0; i < covariate_terms.size();i++){
    	if(!StringCmp(covariate_terms[i].c_str(), "sex", case_ins)){
    		if((covariate_term_matrix.col(i).array() == 2.0).count() != 0){
    			for(int row = 0 ; row < covariate_term_matrix.rows(); row++){
    				if(covariate_term_matrix(row,i) == 2)
    					covariate_term_matrix(row,i) = 1;
    				else
    					covariate_term_matrix(row,i) = 0;
    			}
    		}
    		break;
    	}
    }
    ids = temp_ids;
    return 0;
}
static Eigen::MatrixXd create_covariate_matrix(vector<string> ids, vector<string> covariate_ids, vector<string> covariate_terms, Eigen::MatrixXd raw_covariate_term_matrix, const int n_covariates){
	
	Eigen::MatrixXd ordered_covariate_term_matrix(ids.size(), covariate_terms.size());
	for(int i = 0; i < ids.size(); i++){
		vector<string>::iterator find_iter = find(covariate_ids.begin(), covariate_ids.end(), ids[i]);
		if(find_iter !=  covariate_ids.end()){
			const int row_index = distance(covariate_ids.begin(), find_iter);
			for(int col = 0 ; col < covariate_terms.size(); col++){
				ordered_covariate_term_matrix(i, col) = raw_covariate_term_matrix(row_index, col);
			}
		}else{
			Eigen::MatrixXd null_matrix;
			return null_matrix;
		}
	}
	
	for(int i = 0; i < covariate_terms.size(); i++){
		if(StringCmp(covariate_terms[i].c_str(), "sex", case_ins)){
			ordered_covariate_term_matrix.col(i) = ordered_covariate_term_matrix.col(i).array() - ordered_covariate_term_matrix.col(i).mean();
		}
	}
	
    Eigen::MatrixXd covariate_matrix = Eigen::MatrixXd::Ones(ids.size(), n_covariates);	
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
                covariate_matrix.col(col) = covariate_matrix.col(col).array()*ordered_covariate_term_matrix.col(index).array();
            }else{
                covariate_matrix.col(col) = covariate_matrix.col(col).array()*pow(ordered_covariate_term_matrix.col(index).array(), cov_term->exponent);
            }
        }
        
    }
    
    return covariate_matrix;	


}
typedef struct gwas_screen_data{
	double beta_se;
	double beta;
	double chi;

	gwas_screen_data & operator = (const gwas_screen_data & var) {
		beta = var.beta;
	
		beta_se = var.beta_se;
		chi = var.chi;

	}
	
}gwas_screen_data;
static vector<gwas_screen_data> calculate_gwas_screen(Eigen::VectorXd Y, Eigen::MatrixXd SNP_matrix, Eigen::VectorXd Sigma, const unsigned n_snps){
    vector<gwas_screen_data> results(n_snps);
#pragma omp parallel for schedule(dynamic)
    for(int iteration = 0; iteration < n_snps; iteration++){
    	double denom = Sigma.dot(SNP_matrix.col(iteration).cwiseAbs2());
    	double numerator = SNP_matrix.col(iteration).dot(Y.cwiseProduct(Sigma));
    	results[iteration].beta_se = pow(denom, -0.5);
    	results[iteration].beta = numerator/denom;
    	results[iteration].chi = numerator*numerator/denom;
    }
    
    return results;


}


static const char * run_gwas_screen_list(const char * phenotype_filename, const char * list_filename, const char * evd_data_filename, const char * plink_filename, const bool verbose, unsigned batch_size = GWAS_BATCH_SIZE){

	vector<string> trait_list;
	if(list_filename) {
		trait_list = read_trait_list(list_filename);
	}else{
		string trait_name = string(Trait::Name(0));
		trait_list.push_back(trait_name);
	}
	if(trait_list.size() == 0){
		return "No traits read from list file";
	}
	
	pio_file_t * plink_file = new pio_file_t;
	if (pio_open(plink_file, plink_filename) != PIO_OK){
		return "Error opening plink file";
	}
	const unsigned num_plink_samples = pio_num_samples(plink_file);
	const unsigned n_snps = pio_num_loci(plink_file);
	//unsigned batch_size = GWAS_BATCH_SIZE;
	if(batch_size >= n_snps){
		batch_size = n_snps;
	}
	unsigned iterations = ceil(n_snps/batch_size);
	unsigned final_batch_size = n_snps % batch_size;
	if(final_batch_size == 0) final_batch_size = batch_size;
	if(n_snps >= GWAS_BATCH_SIZE){
		batch_size = GWAS_BATCH_SIZE;
		iterations = ceil(n_snps/GWAS_BATCH_SIZE);
		final_batch_size = n_snps % batch_size;
		if(final_batch_size == 0 ) final_batch_size = batch_size;
	}
	pio_sample_t * sample;
	vector<string> plink_ids;
	for(unsigned i = 0; i < num_plink_samples; i++){
		sample = pio_get_sample(plink_file, i);
		plink_ids.push_back(string(sample->iid));
	}
	
	vector<string> snp_names;
	for(unsigned snp = 0; snp < n_snps; snp++){
		pio_locus_t * locus = pio_get_locus(plink_file, snp);
		snp_names.push_back(string(locus->name));
	}
	//snp_t * snp_buffer = new snp_t[num_plink_samples];
	Solar_Trait_Reader * trait_reader;
	if(evd_data_filename){
		try{
			trait_reader = new Solar_Trait_Reader(phenotype_filename,evd_data_filename, trait_list);
		}catch(Solar_Trait_Reader_Exception & e){
			const char * error_message = e.what();
			return error_message;
		}catch(...){
			return "Unknown error occurred reading phenotype or pedigree data";
		}
	}else{

		try{
			trait_reader = new Solar_Trait_Reader(phenotype_filename, trait_list, plink_ids);
		}catch(Solar_Trait_Reader_Exception & e){
			return e.what();
		}catch(...){
			return "Unknown error occurred reading phenotype or pedigree data";
		}
	}
	for(unsigned set = 0; set < trait_reader->get_n_sets(); set++){
		Eigen_Data * eigen_data = trait_reader->get_eigen_data_set(set);
		vector<string> ids = eigen_data->get_ids();
		unsigned plink_index_map[ids.size()];
		for(unsigned index = 0; index < ids.size(); index++){
			string id = ids[index];
			plink_index_map[index] = distance(plink_ids.begin(), find(plink_ids.begin(), plink_ids.end(), id));
		}
		Eigen::VectorXd eigenvalues = Eigen::Map<Eigen::VectorXd>(eigen_data->get_eigenvalues(), ids.size());
		Eigen::MatrixXd eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>(eigen_data->get_eigenvectors_transposed(), ids.size(), ids.size());
		//eigenvectors_transposed = Eigen::MatrixXd::Identity(eigenvalues.rows(), eigenvalues.rows());
		//eigenvalues = Eigen::VectorXd::Ones(eigenvectors_transposed.rows());

		int * snp_data;
				
		for(unsigned trait = 0; trait < eigen_data->get_n_phenotypes(); trait++){
			string trait_name = eigen_data->get_trait_name(trait);
			string output_filename = trait_name + "-screen-gwas.out";
			ofstream output_stream(output_filename.c_str());
			output_stream << "SNP,p-value,beta,beta_se\n";

			Eigen::VectorXd trait_vector = Eigen::Map<Eigen::VectorXd>(eigen_data->get_phenotype_column(trait), ids.size());
			pio_reset_row(plink_file);
			unsigned current_batch_size = batch_size;
			
    		
			Eigen::VectorXd default_Y;// = trait_vector;
    			
    			Eigen::MatrixXd default_U;
   		


			default_Y = eigenvectors_transposed*(trait_vector.array() - trait_vector.mean()).matrix();
			default_U = Eigen::MatrixXd::Ones(eigenvalues.rows(),2 );
			default_U.col(1) = eigenvalues;				
			Eigen::VectorXd theta = calculate_theta(default_Y, default_U);
			Eigen::VectorXd Sigma = (default_U*theta).cwiseInverse();
			//}else{
			//	default_Y = trait_vector;
			//}


    			//gwas_data default_null_result = compute_null_model_MLE(default_Y, default_mean, default_U, precision);
    			unsigned n_snps_computed = 0;
    			vector<gwas_screen_data> results;
    	    int max_output_width = 0;
			for(unsigned iteration = 0; iteration < iterations; iteration++){
				if(iteration == iterations - 1) current_batch_size = final_batch_size;
		
				Eigen::MatrixXd snp_matrix(ids.size(), current_batch_size);
					
				unsigned snp_index = 0;
				#pragma omp parallel for
				 for(unsigned snp = 0; snp < current_batch_size; snp++){
					snp_t snp_buffer[num_plink_samples];
					unsigned local_snp_index;
				 	#pragma omp critical
					{
						pio_next_row(plink_file, &snp_buffer[0]);
						local_snp_index = snp_index;
						snp_index++;
					}
					double mean = 0.0;
					unsigned current_n_subjects = ids.size();
					for(unsigned id = 0; id < ids.size(); id++){
						const double value = snp_buffer[plink_index_map[id]];
						if(value != 3){
							mean += value;
						}else{
							current_n_subjects--;
						}
						snp_matrix(id, local_snp_index) = value;
							//snp_data[ids.size()*snp + id] = snp_buffer[plink_index_map[id]];
					}
						mean /= current_n_subjects;
					for(unsigned id = 0; id < ids.size() ;id++){
						if(snp_matrix(id,local_snp_index) != 3) 
							snp_matrix(id,local_snp_index) -= mean;
						else
							snp_matrix(id, local_snp_index) = 0;
					}
				}
					snp_matrix = eigenvectors_transposed*snp_matrix;
				results = calculate_gwas_screen(default_Y, snp_matrix,  Sigma, current_batch_size);

				for(unsigned snp = 0; snp < current_batch_size; snp++){
					double pvalue = 1.0;
					if(results[snp].chi > 0.0){
						 pvalue =  chicdf(results[snp].chi, 1); 
						output_stream << snp_names[iteration*batch_size + snp] << "," << pvalue << "," << \
						results[snp].beta << "," << results[snp].beta_se << "\n";
					}else{
						output_stream << snp_names[iteration*batch_size + snp] << "," << 1.0 << "," << \
						0.0 << "," << 0.0 << "\n";
					}
				}
				if(verbose){
					n_snps_computed += current_batch_size;
					std::string output_str = "Trait: " + trait_name + " SNPs Computed: " + to_string(n_snps_computed) + " Percent Complete " + to_string(floor(100.0*n_snps_computed/n_snps)) + "%\r";	
					if(max_output_width < output_str.length()) max_output_width = output_str.length();	
					//std::cout << "Trait: " << trait_name <<  " SNPs Computed: " << n_snps_computed << " SNPs Left: " << n_snps - n_snps_computed << " Percent Complete " << 100*n_snps_computed/n_snps << "% \r";
					std::cout << std::setw(max_output_width) << output_str << std::flush;
					//std::cout.flush();	
				}

			}
			if(verbose){
				std::cout.flush();
				std::cout << "\n";
				std::cout << "Trait: " << trait_name << " is finished GWAS Screen computation\n";
			}
			
			output_stream.close();
		}
		
		

	}		
	
}
static const char * run_gwas_list(const char * phenotype_filename, const char * list_filename, const char * evd_data_filename, const char * plink_filename, const bool fix_missing, const unsigned precision, const bool verbose, const bool use_covariates, unsigned n_permutations = 0, unsigned batch_size = GWAS_BATCH_SIZE){
	vector<string> trait_list;// = read_trait_list(list_filename);
	if(list_filename) {
		trait_list = read_trait_list(list_filename);
	}else{
		string trait_name = string(Trait::Name(0));
		trait_list.push_back(trait_name);
	}
	int n_covariates = 0;
	vector<string> covariate_terms;
	Eigen::MatrixXd raw_covariate_term_matrix;
	vector<string> covariate_term_ids;
	if(use_covariates){
	    Covariate * c;
    
   	  
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
    	const char * error_message = 0;
    	error_message = load_covariate_terms(phenotype_filename, raw_covariate_term_matrix, covariate_terms, covariate_term_ids);
    		
    	if(error_message) return error_message;
    	
    }
	
	if(trait_list.size() == 0){
		return "No traits read from list file";
	}
	
	pio_file_t * plink_file = new pio_file_t;
	if (pio_open(plink_file, plink_filename) != PIO_OK){
		return "Error opening plink file";
	}
	const unsigned num_plink_samples = pio_num_samples(plink_file);
	const unsigned n_snps = pio_num_loci(plink_file);
	//unsigned batch_size = GWAS_BATCH_SIZE;
	if(batch_size >= n_snps){
		batch_size = n_snps;
	}
	unsigned iterations = ceil(n_snps/batch_size);
	unsigned final_batch_size = n_snps % batch_size;
	if(final_batch_size == 0) final_batch_size = batch_size;
	if(n_snps >= GWAS_BATCH_SIZE){
		batch_size = GWAS_BATCH_SIZE;
		iterations = ceil(n_snps/GWAS_BATCH_SIZE);
		final_batch_size = n_snps % batch_size;
		if(final_batch_size == 0 ) final_batch_size = batch_size;
	}
	pio_sample_t * sample;
	vector<string> plink_ids;
	for(unsigned i = 0; i < num_plink_samples; i++){
		sample = pio_get_sample(plink_file, i);
		plink_ids.push_back(string(sample->iid));
	}
	//std::cout << "plink id count: " << plink_ids.size() << std::endl;
	vector<string> snp_names;
	for(unsigned snp = 0; snp < n_snps; snp++){
		pio_locus_t * locus = pio_get_locus(plink_file, snp);
		snp_names.push_back(string(locus->name));
	}
	//snp_t * snp_buffer = new snp_t[num_plink_samples];
	Solar_Trait_Reader * trait_reader;
	if(evd_data_filename){
		try{
			trait_reader = new Solar_Trait_Reader(phenotype_filename,evd_data_filename, trait_list);
		}catch(Solar_Trait_Reader_Exception & e){
			const char * error_message = e.what();
			return error_message;
		}catch(...){
			return "Unknown error occurred reading phenotype or pedigree data";
		}
	}else{
		vector<string> id_include_list;
		if(use_covariates){
			for(int i = 0; i < plink_ids.size(); i++){
				vector<string>::iterator find_iter = find(covariate_term_ids.begin(), covariate_term_ids.end(), plink_ids[i]);
				if(find_iter !=  covariate_term_ids.end()){
					id_include_list.push_back(plink_ids[i]);
				}
			}
		}else{
			id_include_list = plink_ids;
		}
		try{
			trait_reader = new Solar_Trait_Reader(phenotype_filename, trait_list, id_include_list);
		}catch(Solar_Trait_Reader_Exception & e){
			return e.what();
		}catch(...){
			return "Unknown error occurred reading phenotype or pedigree data";
		}
	}
	if(trait_reader->get_n_sets() == 0){
		return "No viable data could be read";
	}
     //   int * iteration_count = new int[batch_size];

	for(unsigned set = 0; set < trait_reader->get_n_sets(); set++){
		Eigen_Data * eigen_data = trait_reader->get_eigen_data_set(set);
		vector<string> ids = eigen_data->get_ids();
		unsigned plink_index_map[ids.size()];
		for(unsigned index = 0; index < ids.size(); index++){
			string id = ids[index];
			plink_index_map[index] = distance(plink_ids.begin(), find(plink_ids.begin(), plink_ids.end(), id));
		}
		Eigen::VectorXd eigenvalues = Eigen::Map<Eigen::VectorXd>(eigen_data->get_eigenvalues(), ids.size());
		Eigen::MatrixXd eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>(eigen_data->get_eigenvectors_transposed(), ids.size(), ids.size());
		//eigenvectors_transposed = Eigen::MatrixXd::Identity(eigenvalues.rows(), eigenvalues.rows());
		//eigenvalues = Eigen::VectorXd::Ones(eigenvectors_transposed.rows());
		Eigen::MatrixXd phi2;// = eigenvectors_transposed.transpose()*eigenvalues.asDiagonal()*eigenvectors_transposed;
		if(!fix_missing){
			phi2 = eigenvectors_transposed.transpose()*eigenvalues.asDiagonal()*eigenvectors_transposed;
		}
		int * snp_data;
		if(!fix_missing) snp_data = new int[batch_size*ids.size()];
		Eigen::MatrixXd covariate_matrix;
		Eigen::MatrixXd default_covariate_matrix;
		if(use_covariates){
	 		covariate_matrix = create_covariate_matrix(ids, covariate_term_ids, covariate_terms,  raw_covariate_term_matrix, n_covariates);
	 		
	 		if(covariate_matrix.rows() == 0){
	 			return "Failure loading covariates";
	 		}
	 		Eigen::MatrixXd temp_covariate_matrix = Eigen::MatrixXd::Ones(covariate_matrix.rows(), covariate_matrix.cols() +1);	
			default_covariate_matrix.resize(covariate_matrix.rows(), covariate_matrix.cols() + 2);
			for(int col = 0; col < covariate_matrix.cols(); col++){
				temp_covariate_matrix.col(col) = covariate_matrix.col(col);
						
			}
			covariate_matrix = eigenvectors_transposed*temp_covariate_matrix;
			for(int col = 0; col < covariate_matrix.cols(); col++){
				default_covariate_matrix.col(col) = covariate_matrix.col(col);
						
			}		 		
	 	}
 		if(n_permutations){
     			permutated_indices = new int[n_permutations*ids.size()];
     			vector<int> indices(ids.size());
     			for(int i = 0; i < indices.size() ; i++){
     				indices[i] = i;
     			}
     			for(int p = 0; p < n_permutations; p++){
     				random_shuffle(indices.begin(), indices.end());
     				for(int row = 0; row < ids.size(); row++){
     					permutated_indices[p*ids.size() + row] = indices[row];
     				}
     			}
     		}		
		for(unsigned trait = 0; trait < eigen_data->get_n_phenotypes(); trait++){
			string trait_name = eigen_data->get_trait_name(trait);
			string output_filename = trait_name + "-gwas.out";
			ofstream output_stream(output_filename.c_str());
			output_stream << "SNP,h2r,loglik,SD,beta_snp,beta_snp_se,chi2,p-value,Status\n";

			Eigen::VectorXd trait_vector = Eigen::Map<Eigen::VectorXd>(eigen_data->get_phenotype_column(trait), ids.size());
			vector<gwas_data> results;
			pio_reset_row(plink_file);
			unsigned current_batch_size = batch_size;
			
    		
			Eigen::VectorXd default_Y;// = trait_vector;
    			Eigen::VectorXd mean = Eigen::ArrayXd::Ones(trait_vector.rows()).matrix();
    			Eigen::MatrixXd default_U;
   		       // Eigen::MatrixXd default_eigenvectors_transposed;
   		        Eigen::VectorXd default_mean;
   			Eigen::MatrixXd default_phi2;// = phi2;
			//Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(ids.size(),ids.size());
			//Eigen::MatrixXd hat = identity - mean*(mean.transpose()*mean).inverse()*mean.transpose();
			//Eigen::MatrixXd Xphi2X = phi2;
			//Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(phi2);
			//cout << Xphi2X.rows() << endl;
			//Eigen::MatrixXd all_eigenvectors = es.eigenvectors();
			//cout << all_eigenvectors.cols() << endl;
			//Eigen::VectorXd all_eigenvalues = es.eigenvalues();
			//Eigen::VectorXd reml_eigenvalues(all_eigenvalues.rows());
			//Eigen::MatrixXd reml_eigenvectors(all_eigenvalues.rows(), all_eigenvalues.rows());
			//cout << all_eigenvalues(all_eigenvalues.rows() - 1) << endl;
			//for(unsigned row = 0 ; row < all_eigenvalues.rows(); row++){
			//	reml_eigenvalues(row) = all_eigenvalues(row );
			//	reml_eigenvectors.col(row) = all_eigenvectors.col(row);
			//}
			//default_eigenvectors_transposed = eigenvectors_transposed;//es.eigenvectors().transpose();
			
			//default_U = Eigen::MatrixXd::Ones(eigenvectors_transposed.rows(),2 );
			//default_U.col(1) = eigenvalues;//es.eigenvalues();

			//default_Y = default_eigenvectors_transposed*(default_Y.array() - default_Y.mean()).matrix();
      			//calculate_eigen_data_null(default_Y, default_eigenvectors_transposed, default_U, default_phi2);
   			//Eigen::VectorXd default_mean; = eigenvectors_transposed*mean;
			//if(fix_missing){
				if(!use_covariates){
					default_Y = eigenvectors_transposed*trait_vector;
					default_U = Eigen::MatrixXd::Ones(eigenvectors_transposed.rows(),2 );
					default_U.col(1) = eigenvalues;
					default_mean = eigenvectors_transposed*mean;
				}else {
					default_Y = eigenvectors_transposed*trait_vector;
					default_U = Eigen::MatrixXd::Ones(eigenvalues.rows(),2 );
					default_U.col(1) = eigenvalues;	
				}

			//}else{
			//	default_Y = trait_vector;
			//}

    			gwas_data default_null_result;
    			if(!use_covariates) default_null_result = gwas_maximize_newton_raphson_method_null_model(default_Y,default_mean, default_U,  precision);
			if(use_covariates) default_null_result = gwas_maximize_newton_raphson_method_with_covariates_null_model(default_Y,covariate_matrix, default_U,  precision);
    			//gwas_data default_null_result = compute_null_model_MLE(default_Y, default_mean, default_U, precision);
    			unsigned n_snps_computed = 0;
    			//std::cout << "n_subjects : " << ids.size() << std::endl;
    		int max_output_width = 0;
			for(unsigned iteration = 0; iteration < iterations; iteration++){
				if(iteration == iterations - 1) current_batch_size = final_batch_size;
				vector<string> status_vector(current_batch_size);
				if(fix_missing){
					Eigen::MatrixXd snp_matrix(ids.size(), current_batch_size);
					
					unsigned snp_index = 0;
					#pragma omp parallel for
				 	for(unsigned snp = 0; snp < current_batch_size; snp++){
						snp_t snp_buffer[num_plink_samples];
						unsigned local_snp_index;
					#pragma omp critical
						{
							pio_next_row(plink_file, &snp_buffer[0]);
							local_snp_index = snp_index;
							snp_index++;
						}
						double mean = 0.0;
						unsigned current_n_subjects = ids.size();
						for(unsigned id = 0; id < ids.size(); id++){
							const double value = snp_buffer[plink_index_map[id]];
							if(value != 3){
								mean += value;
							}else{
								current_n_subjects--;
							}
							snp_matrix(id, local_snp_index) = value;
							//snp_data[ids.size()*snp + id] = snp_buffer[plink_index_map[id]];
						}
						mean /= current_n_subjects;
						for(unsigned id = 0; id < ids.size() ;id++){
							if(snp_matrix(id,local_snp_index) != 3) 
								snp_matrix(id,local_snp_index) -= mean;
							else
								snp_matrix(id, local_snp_index) = 0;
						}
					}
					snp_matrix = eigenvectors_transposed*snp_matrix;
				
					
					if(!use_covariates){ 
						results = GWAS_MLE_fix_missing_run(default_null_result, default_Y, default_mean,\
						snp_matrix, default_U, \
                                               current_batch_size, precision, status_vector, n_permutations);	
                                        }else{
						results = GWAS_MLE_fix_missing_run_with_covariates(default_null_result, default_Y, default_covariate_matrix,\
						snp_matrix, default_U, \
                                               current_batch_size, precision, status_vector);	
                                        }
				}else{
					snp_t * snp_buffer = new snp_t[num_plink_samples];
				 	for(unsigned snp = 0; snp < current_batch_size; snp++){
						pio_next_row(plink_file, snp_buffer);
						for(unsigned id = 0; id < ids.size(); id++){
							snp_data[ids.size()*snp + id] = snp_buffer[plink_index_map[id]];
						}
					}
					delete [] snp_buffer;
					results = GWAS_MLE_run(default_null_result,trait_vector,  default_Y, default_mean,\
				 default_U,  eigenvectors_transposed, phi2,\
				snp_data, current_batch_size, ids.size(), precision,status_vector);
				}

				for(unsigned snp = 0; snp < current_batch_size; snp++){
					string status_string;
					output_stream << snp_names[iteration*batch_size + snp] << "," << results[snp].h2r << "," << \
					results[snp].loglik << "," << results[snp].SD << "," << results[snp].beta \
					 << "," << results[snp].SE << "," << results[snp].chi << "," << results[snp].pvalue << "," << status_vector[snp] << "\n";

			
				}
				if(verbose){
					n_snps_computed += current_batch_size;		
					std::string output_str = "Trait: " + trait_name + " SNPs Computed: " + to_string(n_snps_computed) + " Percent Complete " + to_string(floor(100.0*n_snps_computed/n_snps)) + "%\r";	
					if(max_output_width < output_str.length()) max_output_width = output_str.length();	
					//std::cout << "Trait: " << trait_name <<  " SNPs Computed: " << n_snps_computed << " SNPs Left: " << n_snps - n_snps_computed << " Percent Complete " << 100*n_snps_computed/n_snps << "% \r";
					std::cout <<  setw(max_output_width) << output_str << std::flush;	
				}

			}
			if(verbose){
				std::cout.flush();
				std::cout << "\n";
				std::cout << "Trait: " << trait_name << " is finished GWAS computation\n";
			}
			
			output_stream.close();
		}
		if(n_permutations) delete [] permutated_indices;
		if(!fix_missing) delete [] snp_data;

	}
	//delete [] iteration_count;
	delete trait_reader;

	return 0;
}
static Eigen::VectorXd create_snp_vector(vector<string> ids, vector<string> snp_ids, vector<double> snp_values){
	Eigen::VectorXd snp_vector(ids.size());
	for(unsigned i = 0; i < ids.size(); i++){
		string current_id = ids[i];
		for(unsigned j = 0; j < snp_ids.size(); j++){
			if(current_id == snp_ids[j]){
				snp_vector[i] = snp_values[j];
				break;
			}

		}

	}

	return snp_vector;
}
static vector<gwas_data> calculate_single_snp_gwas(Eigen::MatrixXd trait_matrix,Eigen::VectorXd mean_column,Eigen::VectorXd snp_vector, Eigen::MatrixXd aux, const unsigned precision){
	vector<gwas_data> results(trait_matrix.cols());
	Eigen::MatrixXd X(trait_matrix.rows(), 2);
	X.col(0) = mean_column;
	X.col(1) = snp_vector;
#pragma omp parallel for
	for(unsigned index = 0; index < trait_matrix.cols(); index++){
		Eigen::VectorXd Y = trait_matrix.col(index);
		gwas_data  null_data = compute_null_model_MLE(Y,  mean_column, aux, precision);
		gwas_data result;
       		if(null_data.loglik == null_data.loglik){
           	    MLE_GWAS(&result, Y, X, aux, null_data, precision);
       		}else{
          	   result.chi  = 0.0;
            	   result.beta = 0.0;
           	   result.loglik = 0.0;
            	   result.SD = 0.0;
                   result.SE = 0.0;
                   result.h2r = 0.0;
                   result.pvalue = 0.0;
      		}
		results[index] = result;
		
	}
	return results;

}
static const char * run_single_snp_gwas(const char * phenotype_filename, const char * snp_name, const char * list_filename, const unsigned precision){
	vector<string> trait_list;// = read_trait_list(list_filename);

	trait_list = read_trait_list(list_filename);
	
	
	if(trait_list.size() == 0){
		return "No traits read from list file";
	}
    const char * errmsg = 0;	
    SolarFile * file = SolarFile::open("gwas", phenotype_filename, &errmsg);
    
    if(errmsg){

        return errmsg;
    }

    
    file->start_setup(&errmsg);
    if(errmsg){

        return errmsg;
    }
    file->setup("id", &errmsg);
    if(errmsg){

        return errmsg;
    }

    file->setup(snp_name, &errmsg);
    if(errmsg){
	return errmsg;
    }
    char ** file_data;
    vector<string> phenotype_ids;
    vector<double> snp_values;
    
     while (0 != (file_data = file->get (&errmsg))){
        
       	
       string snp_value_str = file_data[1];
       if(snp_value_str.length() != 0){
	  if(snp_value_str[0] != '3'){
		snp_values.push_back(stod(snp_value_str));
		phenotype_ids.push_back(string(file_data[0]));
	  }
	}
     }
     Solar_Trait_Reader * trait_reader = new Solar_Trait_Reader(phenotype_filename, trait_list, phenotype_ids);
     string output_filename = string(snp_name) + "-gwas.out";
     ofstream output_stream(output_filename.c_str());
     output_stream << "trait,h2r,loglik,SD,beta_snp,beta_snp_se,chi2,p-value\n";
     unsigned total_index = 0; 
     for(unsigned set = 0; set < trait_reader->get_n_sets(); set++){
	Eigen_Data * data_set = trait_reader->get_eigen_data_set(set);
	vector<string> current_ids = data_set->get_ids();
 	Eigen::MatrixXd eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>(data_set->get_eigenvectors_transposed(), current_ids.size(), current_ids.size());
	Eigen::VectorXd eigenvalues = Eigen::Map<Eigen::VectorXd>(data_set->get_eigenvalues(), current_ids.size());
	Eigen::VectorXd snp_vector = create_snp_vector(current_ids, phenotype_ids, snp_values);
	snp_vector = eigenvectors_transposed*(snp_vector.array() - snp_vector.mean()).matrix();
	Eigen::MatrixXd trait_matrix = Eigen::Map<Eigen::MatrixXd>(data_set->get_phenotype_buffer(), data_set->get_n_subjects(), data_set->get_n_phenotypes());
	trait_matrix = eigenvectors_transposed*trait_matrix;
	Eigen::VectorXd mean_column = eigenvectors_transposed*Eigen::ArrayXd::Ones(trait_matrix.rows()).matrix();
	Eigen::MatrixXd aux = Eigen::ArrayXXd::Ones(trait_matrix.rows(), 2);
	aux.col(1) = eigenvalues;
	vector<gwas_data> results = calculate_single_snp_gwas(trait_matrix, mean_column, snp_vector, aux, precision);
	for(unsigned trait_index = 0; trait_index < trait_matrix.cols(); trait_index++){
		gwas_data current_result = results[trait_index];
		output_stream << trait_list[total_index + trait_index] << "," << current_result.h2r << "," << current_result.loglik << ","\
		<< current_result.SD << "," << current_result.beta << "," << current_result.SE << "," << current_result.chi << "," << \
		current_result.pvalue << "\n";
	}
	total_index += trait_matrix.cols();
     }  
     output_stream.close();
     delete trait_reader;
     return 0;

}

extern "C" int gwaCmd(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[]){
//	int n_permutations = 0;
	const char *  plink_filename = 0;
	
    bool correct_missing = false;
    bool verbose = false;
    bool use_covariates= false;
    bool use_screen_option = false;
	const char * list_filename = 0;
    const char * single_snp_name = 0;
    const char * evd_data_filename = 0;
    unsigned precision = 8;
    unsigned batch_size = 0;
    unsigned n_permutations = 0;
    for(unsigned arg = 1; arg < argc; arg++){
	if(!StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins) || !StringCmp(argv[arg], "help", case_ins)){
		print_gwas_help(interp);
		return TCL_OK;
	}else if((!StringCmp(argv[arg], "-plink", case_ins) || !StringCmp(argv[arg], "--plink", case_ins)) && arg + 1 < argc){
		plink_filename = argv[++arg];
	}else if((!StringCmp(argv[arg], "-batch_size", case_ins) || !StringCmp(argv[arg], "--batch_size", case_ins)) && arg + 1 < argc){
		batch_size = atoi(argv[++arg]);
	}else if((!StringCmp(argv[arg], "-np", case_ins) || !StringCmp(argv[arg], "--np", case_ins)) && arg + 1 < argc){
		n_permutations = atoi(argv[++arg]);
			
	}else if (!StringCmp(argv[arg], "-fix", case_ins) || !StringCmp(argv[arg], "--fix", case_ins) || !StringCmp(argv[arg], "-f", case_ins)){
		correct_missing = true;
	}else if (!StringCmp(argv[arg], "-use_covs", case_ins) || !StringCmp(argv[arg], "--use_covs", case_ins)){
		use_covariates = true;
	}else if (!StringCmp(argv[arg], "-verbose", case_ins) || !StringCmp(argv[arg], "--verbose", case_ins) || !StringCmp(argv[arg], "-v", case_ins)){
		verbose = true;
	}else if (!StringCmp(argv[arg], "-screen", case_ins) || !StringCmp(argv[arg], "--screen", case_ins) || !StringCmp(argv[arg], "-s", case_ins)){
		use_screen_option = true;
	}else if((!StringCmp(argv[arg], "-precision", case_ins) || !StringCmp(argv[arg], "--precision", case_ins)) && arg + 1 < argc){
		precision = atoi(argv[++arg]);
	}else if((!StringCmp(argv[arg], "-evd_data", case_ins) || !StringCmp(argv[arg], "--evd_data", case_ins)) && arg + 1 < argc){
		evd_data_filename = argv[++arg];
	}else if((!StringCmp(argv[arg], "-list", case_ins) || !StringCmp(argv[arg], "--list", case_ins)) && arg + 1 < argc){
		list_filename = argv[++arg];
	}else if((!StringCmp(argv[arg], "-single-snp", case_ins) || !StringCmp(argv[arg], "--single-snp", case_ins)) && arg + 1 < argc){
		single_snp_name = argv[++arg];
	}else{
		RESULT_LIT("Invalid argument entered");
		return TCL_ERROR;
	}
   }
    ///string phen_file = Phenotypes::filenames();
    // research_function(phen_file.c_str(), plink_filename, batch_size, precision);
    // return TCL_OK;
   // research_function(phen_file.c_str(), plink_filename, batch_size, precision);
   // return TCL_OK;
    if(use_covariates && (!correct_missing || (!correct_missing && evd_data_filename) || single_snp_name || use_screen_option)){
    	RESULT_LIT("Covariates can only be used with fix missing option and cannot be used with screen option or single snp computation");
    	return TCL_ERROR;
    }
    if(evd_data_filename && !correct_missing){
    	RESULT_LIT("Fix missing option -f must be used with -evd_data option");
    	return TCL_ERROR;
    }
    if((precision < 1 || precision > 9)){
        RESULT_LIT("Precision must be between 1 and 9");
        return TCL_ERROR;
    }
	Eigen::VectorXd trait_vector;

	int success;
	string phenotype_filename = Phenotypes::filenames(); 
	if(phenotype_filename.length() == 0){
		RESULT_LIT("No phenotype file is currently loaded");
		return TCL_ERROR;
	}
	if(!plink_filename && !single_snp_name){
		RESULT_LIT("No plink file was specified");
		return TCL_ERROR;
	}				
	
	if(!list_filename && Trait::Number_Of() == 0){
		RESULT_LIT("No trait selected with trait command or specified with list file");
		return TCL_ERROR;
	}
	if(!evd_data_filename){
		try{
			load_phi2_matrix(interp);
		}catch(...){
			RESULT_LIT("phi2 matrix could not be loaded.  Check to see if pedigree has been properly loaded.");
			return TCL_ERROR;
		}
	}
	const char * error = 0;
	if(single_snp_name){
		if(!list_filename){
			RESULT_LIT("Single SNP gwas requires a trait list specified by --list <list file name>");
			return TCL_ERROR;
		}
		error = run_single_snp_gwas(phenotype_filename.c_str(), single_snp_name, list_filename, precision);

	
	}else if(!use_screen_option){
		if(batch_size == 0){
			error = run_gwas_list(phenotype_filename.c_str(), list_filename,evd_data_filename,\
				 plink_filename, correct_missing, precision, verbose, use_covariates, n_permutations);	
		}else{
			error = run_gwas_list(phenotype_filename.c_str(), list_filename,evd_data_filename,\
				 plink_filename, correct_missing, precision,verbose, use_covariates, n_permutations ,batch_size);
		}
	}else{
		if(batch_size == 0){
			error = run_gwas_screen_list(phenotype_filename.c_str(),  list_filename,  evd_data_filename,  plink_filename,  verbose);
		}else{
			error = run_gwas_screen_list(phenotype_filename.c_str(),  list_filename,  evd_data_filename,  plink_filename,  verbose, batch_size);
		}
	}
		
		
	if(error){
		//cout << error << endl;
		RESULT_BUF(error);
		return TCL_ERROR;
	}
										 
											
    return TCL_OK;
    
}

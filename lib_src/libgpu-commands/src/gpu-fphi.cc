#include "gpu-exception.h"
#include "gpu-fphi-settings.h"
#include "gpu-fphi-variables.h"
#include "solar-trait-reader.h"
#include <iostream>
#include <fstream>
#include "RicVolumeSet.h"
#include <string>
#include <vector>
#include "solar.h"
#include <Eigen/Dense>
using namespace std;
extern double gwas_chicdf(double, double);
static int get_total_dynamic_memory(const int pitch){
	int memory = 0; 
	memory += GPU_FPHI_Stream_Variables::Adjustable_GPU_Memory_Cost(pitch);
	memory += GPU_FPHI_Results::Adjustable_GPU_Memory_Cost();
	return memory;	
		
}
static const char * calculate_pitch_blockSize_and_batch_size(vector<int> gpu_id_list, int n_streams,const bool use_covariates, int & pitch,int & blockSize, int  & batch_size,const int defined_blockSize,\
							 const int defined_batch_size, const int n_subjects, const int n_traits){

    int blockSize_list[6] = {32, 64, 128, 256, 512, 1024};
    int remainder_list[6];
    
    int scale_list[6];
  //  cout << "n_subjects " << n_subjects << endl;
    pitch = ceil(n_subjects/32.f)*32;
    for(int index = 0; index < 6; index++){
	scale_list[index] = ceil((double)pitch/blockSize_list[index]);


	
	
	remainder_list[index] = scale_list[index]*blockSize_list[index] - pitch;
	
    }
    int selection_index = 0;
    for(int index = 1; index < 6; index++){
	if(remainder_list[selection_index] > remainder_list[index]) selection_index = index;
    }
  
    if(scale_list[selection_index] >= 20 && selection_index < 5){
	int new_selection_index;// = selection_index + 1;
	for(new_selection_index = selection_index + 1; new_selection_index < 6; new_selection_index++){
		if(scale_list[new_selection_index] < 20) {
			selection_index = new_selection_index;
			break;
		}
	}
    }

    if(defined_blockSize != -1){
	for(int index = 0; index < 6; index++){
		if(blockSize_list[index] == defined_blockSize){
			selection_index = index;
			break;
		}

	}
    }
    for(int index = 0; index < 6; index++){
    	if(pitch <= blockSize_list[index]){
    		if(blockSize_list[selection_index] > blockSize_list[index]){
    			selection_index = index;
    			break;
    		}
    	}
    }
    blockSize = blockSize_list[selection_index];
 
    
    size_t minimum_gpu_memory = 0;
    
    for(size_t idx = 0; idx < gpu_id_list.size(); idx++){
    	try{
        	cudaErrorCheck(cudaSetDevice(gpu_id_list[idx]));
        }catch(GPU_Exception & e){
        	return e.what();
        }
        size_t total_mem;
        size_t free_mem;
        try{
        	cudaErrorCheck(cudaMemGetInfo(&free_mem,&total_mem));
        }catch(GPU_Exception & e){
        	return e.what();
        }
        if(minimum_gpu_memory == 0) minimum_gpu_memory = free_mem;
        if(minimum_gpu_memory > free_mem) minimum_gpu_memory = free_mem;
        
    }
    
   // const size_t multiple = ceil(n_subjects/32.f);
   // pitch = multiple*32;
    const int static_memory = GPU_FPHI_Shared_Variables::Static_GPU_Memory_Cost(pitch, use_covariates);
    const int dynamic_memory = get_total_dynamic_memory(pitch);

    
    int temp_batch_size = floor(double(minimum_gpu_memory - static_memory)/double(n_streams*dynamic_memory)) - 100;
    if(temp_batch_size > MAX_BATCH_SIZE)
        temp_batch_size = MAX_BATCH_SIZE;
    if(defined_batch_size != -1){
    	 if(defined_batch_size > temp_batch_size){
    	 	std::cout << "Defined batch size is too large for memory considerations defaulting to smallest workable batch size\n";
    	 }else{
    	 	temp_batch_size = defined_batch_size;
    	 }
    }
    if(temp_batch_size*n_streams*gpu_id_list.size() > n_traits){
        temp_batch_size = ceil(n_traits/double(n_streams*gpu_id_list.size()));
    }
    
    batch_size = temp_batch_size;	
    return 0;
}
static vector<string> get_covariate_term_list(int & n_covariates){
	Covariate * c;
	n_covariates = 0;
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
        n_covariates++;
    	}
    	
    	return covariate_terms;

}
static const char * get_covariate_term_data(const char * phenotype_filename,Eigen::MatrixXd & covariate_term_matrix, vector<string> & covariate_ids, vector<string> covariate_term_names){
	const char * errmsg = 0;
	SolarFile * file = SolarFile::open("fphi", phenotype_filename, &errmsg);	
	
	if(errmsg) return errmsg;
	file->start_setup(&errmsg);
	if(errmsg) return errmsg;
	file->setup("id", &errmsg);
	if(errmsg) return errmsg;
	const int n_cols = covariate_term_names.size();
	int sex_index = -1;
	for(int i = 0 ; i < n_cols;i++){
		if(!StringCmp("sex", covariate_term_names[i].c_str(), case_ins)){
			sex_index = i + 1;
		}
		file->setup(covariate_term_names[i].c_str(), &errmsg);
		if(errmsg) return errmsg;
	}
	vector<vector<double> > covariate_data;
	char ** file_data;
	while (0 != (file_data = file->get (&errmsg))){
		string current_id = file_data[0];
		bool skip_row = false;
		vector<double> row_data(n_cols);
		for(int index = 1; index < n_cols + 1; index++){
			if(StringCmp(file_data[index], 0, case_ins)){
				if(sex_index != index){
					row_data[index - 1] = atof(file_data[index]);
				}else {
					if(!StringCmp(file_data[index],"f", case_ins)){
						row_data[index - 1] = 1;
					}else if(!StringCmp(file_data[index], "m", case_ins)){
						row_data[index - 1] = 0;
					}else{
						row_data[index - 1] = atof(file_data[index]);
					}
				}
			}else{
				skip_row = true;
				break;
			}
		}
		if(!skip_row){
			covariate_ids.push_back(current_id);
			covariate_data.push_back(row_data);
		}
	}	
				
		
		
	delete file;
	const int n_rows = covariate_ids.size();
	covariate_term_matrix.resize(n_rows, n_cols);
	for(int row = 0; row < n_rows; row++){
		vector<double> row_data = covariate_data[row];
		for(int col = 0; col < n_cols; col++){
			covariate_term_matrix(row,col) = row_data[col];
		}
	}
	
	return 0;	

}


static float * calculate_covariate_matrix(Eigen::MatrixXd eigenvectors_transposed, vector<string> ids, vector<string> covariate_ids, vector<string> covariate_term_names, Eigen::MatrixXd covariate_term_matrix, \
						const int n_covariates, const int pitch){
	
	Eigen::MatrixXd ordered_covariate_term_matrix(ids.size(), covariate_term_names.size());
	for(int i = 0; i < ids.size(); i++){
		vector<string>::iterator find_iter = find(covariate_ids.begin(), covariate_ids.end(), ids[i]);
		if(find_iter !=  covariate_ids.end()){
			const int row_index = distance(covariate_ids.begin(), find_iter);
			for(int col = 0 ; col < covariate_term_names.size(); col++){
				ordered_covariate_term_matrix(i, col) = covariate_term_matrix(row_index, col);
			}
		}else{
			float * null_matrix = 0;
			return null_matrix;
		}
	}
	
	for(int i = 0; i < covariate_term_names.size(); i++){
		if(StringCmp(covariate_term_names[i].c_str(), "sex", case_ins)){
			ordered_covariate_term_matrix.col(i) = ordered_covariate_term_matrix.col(i).array() - ordered_covariate_term_matrix.col(i).mean();
		}
	}
    	
    Eigen::MatrixXd covariate_matrix = Eigen::MatrixXd::Ones(ids.size(), n_covariates + 1);	
    Covariate * cov;
    
    for(int col = 0; (cov = Covariate::index(col)); col++){
        CovariateTerm * cov_term;
        for(cov_term = cov->terms(); cov_term; cov_term = cov_term->next){
            int index = 0;
            
            for(vector<string>::iterator cov_iter = covariate_term_names.begin(); cov_iter != covariate_term_names.end(); cov_iter++){
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
    
   covariate_matrix = 	eigenvectors_transposed*covariate_matrix;
   Eigen::MatrixXd hat_matrix = Eigen::MatrixXd::Identity(covariate_matrix.rows(), covariate_matrix.rows()) - covariate_matrix*(covariate_matrix.transpose()*covariate_matrix).inverse()*covariate_matrix.transpose();
   
   float * output_hat_matrix = new float[pitch*pitch];
   memset(output_hat_matrix, 0, sizeof(float)*pitch*pitch);
   for(int col = 0; col < ids.size(); col++){
   	for(int row = 0 ;row < ids.size(); row++){
   		output_hat_matrix[col*pitch + row] = hat_matrix(row,col);
   	}
   }

   return output_hat_matrix;

}

static const char * write_output_to_image_volume_sets(string template_filename, string output_filename, vector<string> trait_names, float * h2r, float * chi_squared, float * standard_error){
    RicVolumeSet * template_volume;
    try{
    	template_volume = new RicVolumeSet(template_filename);
    }catch(...){
    	return "Error opening template nifti volume output cannot be written";
    }
    const int nx = template_volume->nx;
    const int ny = template_volume->ny;
    const int nz = template_volume->nz;
    RicVolumeSet * h2r_volume = new RicVolumeSet(nx, ny, nz, 1);
    h2r_volume->NIFTIorientation = template_volume->NIFTIorientation;
    
    RicVolumeSet * standard_error_volume = new RicVolumeSet(nx, ny, nz, 1);
    standard_error_volume->NIFTIorientation = template_volume->NIFTIorientation;
    
    RicVolumeSet * p_value_volume = new RicVolumeSet(nx, ny, nz, 1);
    p_value_volume->NIFTIorientation = template_volume->NIFTIorientation; 
    h2r_volume->dtype = DT_FLOAT;
    standard_error_volume->dtype = DT_FLOAT;
    p_value_volume->dtype = DT_FLOAT;   
    for(int trait = 0; trait < trait_names.size() ; trait++){
    	string current_trait = trait_names[trait];
    	if(current_trait.substr(0,5) != "VOXEL"){
    		cout << current_trait << " does not contain voxel naming format VOXEL_x_y_z and will be excluded from output volumes\n";
    		continue;
    	}
    	int x,y,z;
    	int index = 6; 
    	string str;
    	if(current_trait.length() < 7){
    		cout << current_trait << " does not contain voxel naming format VOXEL_x_y_z and will be excluded from output volumes\n";
    		continue;
    	}
    	bool skip_trait = false;
    	while(current_trait[index] != '_'){
    		str+= current_trait[index++];
    		if(index == current_trait.length()){
    			skip_trait = true;
    		}
    	}
    	if(skip_trait){
    		cout << current_trait << " does not contain voxel naming format VOXEL_x_y_z and will be excluded from output volumes\n";
    		continue;
    	}    	
    	x = stoi(str);
    	index++;
    	str.clear();
    	if(current_trait.length() < index){
    		cout << current_trait << " does not contain voxel naming format VOXEL_x_y_z and will be excluded from output volumes\n";
    		continue;
    	}    	
    	while(current_trait[index] != '_'){
    		str += current_trait[index++];
    		if(index == current_trait.length()){
    			skip_trait = true;
    		}
    	}
    	if(skip_trait){
    		cout << current_trait << " does not contain voxel naming format VOXEL_x_y_z and will be excluded from output volumes\n";
    		continue;
    	}     	
    	y = stoi(str);
    	index++;
    	str.clear();
    	if(current_trait.length() < index){
    		cout << current_trait << " does not contain voxel naming format VOXEL_x_y_z and will be excluded from output volumes\n";
    		continue;
    	}    	
    	while(current_trait.length() != index){
    		str += current_trait[index++];	
    	}
    	z = stoi(str);
    	if(x >= nx || y >= ny || z >= nz){
    		cout << "Voxel with indices x: " << x << " y: " << y << " z: " << z << " is out of bounds set by template volume some must be skipped\n";
    		continue;
    	}
    	const float l_h2r = h2r[trait];
    	const float l_chi_squared = chi_squared[trait];
    	if(l_h2r == l_h2r){
    		float pvalue = 0.5f;
    		if(l_chi_squared >= 0){
    			pvalue =  gwas_chicdf(l_chi_squared, 1)/2.f;
    		}
    		 h2r_volume->VolSet[0].vox[x][y][z] = l_h2r;
    		 standard_error_volume->VolSet[0].vox[x][y][z] = standard_error[trait];
    		 p_value_volume->VolSet[0].vox[x][y][z] = pvalue;
    	}else{
     		 h2r_volume->VolSet[0].vox[x][y][z] = 0;
    		 standard_error_volume->VolSet[0].vox[x][y][z] = 0;
    		 p_value_volume->VolSet[0].vox[x][y][z] = .5f;
    	}   		
    }
    string base_filename =  output_filename + ".nii.gz";
    h2r_volume->Write("h2r-volume-"+base_filename);
    p_value_volume->Write("p-value-volume-" + base_filename);
    standard_error_volume->Write("standard-error-volume-"+base_filename);
    delete template_volume;
    delete h2r_volume;
    delete p_value_volume;
    delete standard_error_volume;
    return 0;

}
static void write_output(const char * output_filename, vector<string> trait_names, float * h2r, float * chi_squared, float * standard_error){
	ofstream output_stream(output_filename);
	output_stream << "Trait,h2r,standard_error,p-value\n";
	for(int trait = 0; trait < trait_names.size(); trait++){
		float l_h2r = h2r[trait];
		float l_chi_squared = chi_squared[trait];
		if(l_h2r == l_h2r && l_h2r >= 0.f && l_h2r <= 1.f){
			float pvalue = 0.5f;
			if(l_chi_squared >= 0){
				pvalue =  gwas_chicdf(l_chi_squared, 1)/2.f;
			}
			output_stream << trait_names[trait]  <<  "," << l_h2r << "," << standard_error[trait] << "," << pvalue << endl;
		}else{
			output_stream << trait_names[trait]  <<  "," << 0 << "," << 0 << "," << 0.5f << endl;
		}		
	}
	output_stream.close();
}
const char * run_gpu_fphi(const char * phenotype_filename, const char * output_filename, const char * evd_base_filename, const char * template_filename, vector<string> raw_trait_list, vector<int> gpu_id_list,  int n_streams,\
			  int defined_blockSize, int defined_batch_size,  bool use_covariates, const bool verbose){
	const char * error_message = 0;
	const int n_devices = gpu_id_list.size();
	vector<string> covariate_ids;
	Eigen::MatrixXd raw_covariate_term_matrix;
	vector<string> covariate_terms;
	int n_covariates;
	if(use_covariates){
		covariate_terms= get_covariate_term_list(n_covariates);
		if(n_covariates == 0){
			use_covariates = false;
			covariate_ids.clear();
		}
			
		if(use_covariates) error_message = get_covariate_term_data(phenotype_filename, raw_covariate_term_matrix,covariate_ids,  covariate_terms);
		if(error_message) return error_message;
	}
	
	Solar_Trait_Reader * trait_reader;
	if(verbose){
		cout << "Reading phenotype file for traits\n";
	}
	if(!evd_base_filename){
		try{
			trait_reader = new Solar_Trait_Reader(phenotype_filename, raw_trait_list, covariate_ids);
		}catch(Solar_Trait_Reader_Exception & e){
	
			return e.what();
		}
	}else{
		try{
			trait_reader = new Solar_Trait_Reader(phenotype_filename, evd_base_filename, raw_trait_list);
		}catch(Solar_Trait_Reader_Exception & e){
	
			return e.what();
		}
	}		
	const int total_phenotypes= trait_reader->get_n_phenotypes();
	float * h2r;// = new float[total_phenotypes];
	float * chi_squared;// = new float[total_phenotypes];
	float * standard_error;// = new float[total_phenotypes];
	try{
		cudaErrorCheck(cudaMallocHost((void**)&h2r, sizeof(float)*total_phenotypes));
	}catch(GPU_Exception & e){
		return e.what();
	}
	try{
		cudaErrorCheck(cudaMallocHost((void**)&chi_squared, sizeof(float)*total_phenotypes));
	}catch(GPU_Exception & e){
		return e.what();
	}
	try{
		cudaErrorCheck(cudaMallocHost((void**)&standard_error, sizeof(float)*total_phenotypes));
	}catch(GPU_Exception & e){
		return e.what();
	}		

	int pitch, blockSize, batch_size;
	int result_start_index = 0;
	vector<string> trait_names;
	for(int set = 0; set < trait_reader->get_n_sets(); set++){
		Eigen_Data * eigen_data = trait_reader->get_eigen_data_set(set);
		const int n_subjects = eigen_data->get_n_subjects();
		double * raw_eigenvectors_transposed = eigen_data->get_eigenvectors_transposed();
		Eigen::MatrixXd eigen_eigenvectors_transposed = Eigen::Map<Eigen::MatrixXd>(raw_eigenvectors_transposed, n_subjects, n_subjects);
		const int n_phenotypes = eigen_data->get_n_phenotypes();
		error_message = calculate_pitch_blockSize_and_batch_size(gpu_id_list, n_streams, use_covariates, pitch, blockSize, batch_size,defined_blockSize,\
							  defined_batch_size,  n_subjects,  n_phenotypes);
		if(error_message) return error_message;					
		float * hat_matrix = 0;
		if(use_covariates){
			hat_matrix = calculate_covariate_matrix(eigen_eigenvectors_transposed, eigen_data->get_ids(), covariate_ids, covariate_terms, raw_covariate_term_matrix,  n_covariates,  pitch);
			if(hat_matrix == 0){
				cout << "Unable to produce covariates proceeding without them\n";
				use_covariates = false;
			}
		}		
		
		

		//double * raw_phenotype_data = eigen_data->get_phenotype_buffer();
		
		
		double * raw_eigenvalues = eigen_data->get_eigenvalues();
		Eigen::MatrixXd auxiliary = Eigen::MatrixXd::Ones(n_subjects, 2);
		auxiliary.col(1) = Eigen::Map<Eigen::VectorXd>(raw_eigenvalues, n_subjects);
		Eigen::MatrixXd ZTZI = (auxiliary.transpose()*auxiliary).inverse();
		const float ZTZI_0 = ZTZI(0,0);
		const float ZTZI_1 = ZTZI(1,0);
		const float ZTZI_2 = ZTZI(0,1);
		const float ZTZI_3 = ZTZI(1,1);
		float * eigenvectors_transposed = new float[pitch*pitch];
		float * eigenvalues = new float[pitch];
		memset(eigenvalues, 0, sizeof(float)*pitch);
		memset(eigenvectors_transposed, 0, sizeof(float)*pitch*pitch);
		for(int col = 0; col < n_subjects; col++){
			eigenvalues[col] = raw_eigenvalues[col];
			for(int row = 0; row < n_subjects; row++){
				eigenvectors_transposed[pitch*col + row] = raw_eigenvectors_transposed[n_subjects*col + row];
			}
		}
		vector<string> set_trait_names = eigen_data->get_trait_names();
		for(int trait = 0; trait < set_trait_names.size(); trait++){
			trait_names.push_back(set_trait_names[trait]);
		} 

		
		
		GPU_Data ** gpu_data = new GPU_Data*[n_streams*n_devices];
		GPU_FPHI_Shared_Variables ** shared_variables = new GPU_FPHI_Shared_Variables*[n_devices];
		GPU_FPHI_Stream_Variables ** stream_variables = new GPU_FPHI_Stream_Variables*[n_devices*n_streams];
		GPU_FPHI_Results ** results = new GPU_FPHI_Results*[n_devices*n_streams];
		for(int gpu_index  = 0 ; gpu_index < n_devices; gpu_index++){
			try{
				cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
			}catch(GPU_Exception & e){
				return e.what();
			}
			try{
				shared_variables[gpu_index] = new GPU_FPHI_Shared_Variables(gpu_id_list[gpu_index], pitch, n_subjects, eigenvalues, eigenvectors_transposed,  hat_matrix,ZTZI_0, \
				ZTZI_1,ZTZI_2, ZTZI_3);
			}catch(GPU_Exception & e){
				return e.what();
			}
			for(int stream_index = 0; stream_index < n_streams; stream_index++){

				try{
					gpu_data[gpu_index*n_streams + stream_index] =  new GPU_Data(gpu_id_list[gpu_index], blockSize, batch_size);
				}catch(GPU_Exception & e){
					return e.what();
				}
				try{
					stream_variables[gpu_index*n_streams + stream_index] = new GPU_FPHI_Stream_Variables(gpu_data[gpu_index*n_streams + stream_index], pitch, n_subjects, batch_size);
				}catch(GPU_Exception & e){
					return e.what();
				}
				try{
					results[gpu_index*n_streams + stream_index] =  new GPU_FPHI_Results(h2r + result_start_index, chi_squared + result_start_index, standard_error + result_start_index, batch_size);
				}catch(GPU_Exception & e){
					return e.what();
				}
			}
		}
		GPU_FPHI_Estimator * estimator;
		try{
			estimator =  new GPU_FPHI_Estimator( eigen_data->get_phenotype_buffer(), n_phenotypes, n_devices, n_streams, blockSize, pitch, n_subjects, batch_size, verbose);
    		}catch(GPU_Exception & e){
    			return e.what();
    		}
    		const char ** error_messages = new const char*[n_devices*n_streams];
    		for(int i = 0; i < n_devices*n_streams; i++){
    			error_messages[i] = 0;
    		}
    		bool break_loop = false;
    		string error_message_list;
    		if(verbose){
    			cout << "Starting FPHI estimation of trait set " << set + 1 << " out of " << trait_reader->get_n_sets() << " sets which contains " << n_phenotypes << " traits out of a total of " << total_phenotypes << " traits\n";
    		}	
    		if(n_devices*n_streams == 1){
    			estimator->Thread_Launch(0, 0, gpu_data[0], stream_variables[0], shared_variables[0], results[0], &error_messages[0], break_loop);
    			if(break_loop){
    				error_message_list = string(error_messages[0]) + " error occurred during gpu kernels calls with device id: " + to_string(gpu_id_list[0]) + " and stream index: 0";
    			}
    		}else{
    		    	omp_set_dynamic(0);
        		omp_set_num_threads(n_devices*n_streams);
        	#pragma omp parallel
 		       	{	
				const int thread_index = omp_get_thread_num();	
				const int gpu_index = floor((double)thread_index/n_streams);
				const int stream_index = thread_index % n_streams;
				
				estimator->Thread_Launch(gpu_index, thread_index, gpu_data[gpu_index*n_streams + stream_index],stream_variables[gpu_index*n_streams + stream_index], shared_variables[gpu_index],\
						 results[gpu_index*n_streams + stream_index], &error_messages[gpu_index*n_devices + stream_index], break_loop);
			}
			if(break_loop){
				string error_message_list;
				for(int gpu_index = 0; gpu_index < n_devices; gpu_index++){
					for(int stream_index = 0; stream_index < n_streams; stream_index++){
						if(error_messages[gpu_index*n_devices + stream_index]){
							error_message_list += string(error_messages[gpu_index*n_devices + stream_index]) + " error occurred during gpu kernels calls with device id: " + to_string(gpu_id_list[gpu_index]) \
								 + " and stream_index: " + to_string(stream_index) + "\n";
						} 
					}
				}
			}
		}
		if(error_message_list.length() != 0){
			return error_message_list.c_str();
		}
		for(int gpu_index = 0 ; gpu_index < n_devices; gpu_index++){
			cudaErrorCheck(cudaSetDevice(gpu_id_list[gpu_index]));
			delete shared_variables[gpu_index];
			for(int stream_index = 0; stream_index < n_streams ;stream_index++){
				delete gpu_data[gpu_index*n_streams + stream_index];
				delete stream_variables[gpu_index*n_streams + stream_index];
				delete results[gpu_index*n_streams + stream_index];
				
			}
		}
		
		delete [] shared_variables;
		delete [] gpu_data;
		delete [] stream_variables;
		delete [] results;
		if(hat_matrix) delete [] hat_matrix;
		delete [] eigenvectors_transposed;
		delete [] eigenvalues;
		delete estimator;
		result_start_index += n_phenotypes;
	}
	delete trait_reader;
	if(trait_names.size() != total_phenotypes){
		return "The number of trait names is not equal to the total of phenotypes calculated.  Output cannot be written";
	}
	if(!template_filename){
		 write_output(output_filename,  trait_names, h2r, chi_squared, standard_error);
	}else{
		error_message = write_output_to_image_volume_sets(string(template_filename), string(output_filename), trait_names,  h2r, chi_squared, standard_error);
	}
	try{
		cudaErrorCheck(cudaFreeHost(h2r));
	}catch(GPU_Exception & e){
		return e.what();
	}
	try{
		cudaErrorCheck(cudaFreeHost(chi_squared));
	}catch(GPU_Exception & e){
		return e.what();
	}
	try{
		cudaErrorCheck(cudaFreeHost(standard_error));
	}catch(GPU_Exception & e){
		return e.what();
	}
	
	return error_message;	
	
}

//
//  nifti_to_csv.cpp
//  
//
//  Created by Brian Donohue on 7/26/18.
//
//

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "solar.h"
#include "RicVolumeSet.h"
using namespace std;
class exceptions
{
public:
    virtual const char* what() const throw()
    {
        return "Unkown exception";
    }
} ;
class trait_index_error
{
public:
    trait_index_error(){
    }
    virtual const char* what() const throw()
    {
        return "Trait name not found in phenotype file";
    }
} ;
class out_of_bounds_exception: public exceptions
{
private:
    size_t x;
    size_t y;
    size_t z;
public:
    out_of_bounds_exception(const size_t _x, const size_t _y, const size_t _z){
        x = _x;
        y = _y;
        z = _z;
    }
    virtual const char* what() const throw()
    {
        return string("Index "  + to_string(x)  + " x " + to_string(y) + " y " +  to_string(z) + " z out of bounds").c_str() ;
    }
} ;
 struct voxel_struct{
    size_t x;
    size_t y;
    size_t z;
    string name;
   ~voxel_struct(){
	name.clear();
    }
};

static std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
static std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}
 
static std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, chars), chars);
}
static void  get_mask_vector(vector<int> & x_vector, vector<int> & y_vector, vector<int>  & z_vector, vector<string> &  names,string mask_filename){
    RicVolumeSet mask(mask_filename);
    const int nx = mask.nx;
    const int ny = mask.ny;
    const int nz = mask.nz;
  //  vector<voxel_struct> output_vector;
    for(int x = 0; x < nx; x++)for(int y = 0; y < ny; y++)for(int z = 0; z < nz; z++){
        if(mask.VolSet[0].vox[x][y][z] != 0){
          /*  voxel_struct new_struct;
            new_struct.x = x;
            new_struct.y = y;
            new_struct.z = z;
            new_struct.name*/
	    string name  = "VOXEL_" + to_string(x) + "_" + to_string(y) + "_" + to_string(z);
            names.push_back(name);
	    x_vector.push_back(x);
	    y_vector.push_back(y);
	    z_vector.push_back(z);
        }
    }
    
    
   // return output_vector;
    
}
static float * fill_subject_data(vector<int> mask_x, vector<int> mask_y, vector<int> mask_z, string input_filename){
   // cout << input_filename << endl;
    RicVolumeSet input(input_filename);
    float * buffer = new float[mask_x.size()];
    for(size_t index = 0 ;index < mask_x.size(); index++){
   //     const voxel_struct current_voxel = mask_vector[index];
        const size_t i_x = mask_x[index];
        const size_t i_y = mask_y[index];
        const size_t i_z = mask_z[index];
        if(i_x >= input.nx || i_y >= input.ny || i_z >= input.nz){
            delete [] buffer;
            throw out_of_bounds_exception(i_x, i_y,i_z);
        }
        buffer[index] = input.VolSet[0].vox[i_x][i_y][i_z];
    }
    
    return buffer;
}

static float ** get_subject_data(vector<int> mask_x, vector<int> mask_y, vector<int> mask_z,vector<string> subject_files){
    float ** subject_data = new float*[subject_files.size()];
    
    for(size_t index = 0; index < subject_files.size(); index++){
        subject_data[index] = fill_subject_data(mask_x, mask_y, mask_z, subject_files[index]);
    }
    return subject_data;
}

static void write_data_out_to_single_file(float ** subject_data, vector<string> names, vector<string> subject_lines, string title_line, const char * output_filename, const size_t start, const size_t end){
    ofstream file_out(output_filename);
    file_out << title_line;
    for(size_t index = start; index < end; index++){
        file_out << "," << names[index];
    }
    file_out << endl;
    for(size_t subject = 0; subject < subject_lines.size(); subject++){
        file_out << subject_lines[subject];
        for(size_t index = start; index < end; index++){
            file_out << "," << subject_data[subject][index];
        }
        file_out << endl;
    }
    
    file_out.close();
    
}
static void write_header(vector<string> names, const unsigned start, const unsigned n_voxels, const char * output_filename){
	ofstream file_out(output_filename);
         file_out << names[start];
	for(unsigned i = start + 1 ; i < n_voxels; i++){
		file_out << " " << names[i];
	}
	file_out.close();

}
static void write_data_out_file(float ** subject_data,vector<string> names, vector<string> subject_lines,string title_line, const char * base_filename, const size_t max_voxels_per_file){
    if(max_voxels_per_file == 0 || max_voxels_per_file >= names.size()){
	string header_filename = string(base_filename) + ".header";
	write_header(names, 0, names.size(), header_filename.c_str());
        write_data_out_to_single_file(subject_data, names, subject_lines,title_line, base_filename, 0, names.size());
    }else{
        size_t file_index = 0;
        for(size_t voxels = 0; voxels < names.size(); voxels += max_voxels_per_file){
            string  output_filename = string(base_filename) + "." + to_string(file_index) + ".csv";

	    string header_filename = string(base_filename) + "." + to_string(file_index++) + ".header";
		 
            if(voxels + max_voxels_per_file > names.size()){
		write_header(names, voxels, names.size(), header_filename.c_str());
                write_data_out_to_single_file(subject_data, names, subject_lines,title_line, output_filename.c_str(), voxels, names.size());
            }else{
		 write_header(names, voxels, voxels + max_voxels_per_file, header_filename.c_str());
                write_data_out_to_single_file(subject_data, names, subject_lines,title_line, output_filename.c_str(), voxels, voxels + max_voxels_per_file);
            }
        }
    }
    
    
}
static vector<string> split_line(string line){
    vector<string> line_split;
    string token;
    for(size_t index = 0; index < line.length(); index++){
        const char c = line[index];
        if(c == ','){
            line_split.push_back(trim(token));
            token.clear();
        }else{
            token += c;
        }
    }
    line_split.push_back(trim(token));
    return line_split;
}

static vector<string> get_subject_lines_and_subject_files(const char * phenotype_filename, vector<string> & subject_files,string trait_name, string & title_line){
    
    ifstream file_in(phenotype_filename);
  /*  if(file_in.is_open()){
	cout << "its open\n";
	}else{
	cout << "closed\n";
	}*/
    string raw_title_line;
    getline(file_in, raw_title_line);
    vector<string> split_title_line = split_line(raw_title_line);
    int  trait_index = -1;
    for(int index = 0; index < split_title_line.size(); index++){
        if(!StringCmp(split_title_line[index].c_str(), trait_name.c_str(), case_ins)){
            trait_index = index;
        }else{
            if(title_line.length() != 0){
                title_line +=  "," + split_title_line[index];
            }else{
                title_line = split_title_line[index];
            }
            
        }
    }
    if(trait_index == -1){
       throw trait_index_error();
     }
    string line;
    vector<string> subject_lines;
    while(getline(file_in, line)){
        vector<string> line_split = split_line(line);
        subject_files.push_back(line_split[trait_index]);
	string sum_line;
        for(size_t index = 0; index < line_split.size(); index++){
	 if(index != trait_index){
                if(sum_line.length() != 0){
                    sum_line +=  "," + line_split[index];
                }else{
                    sum_line = line_split[index];
                }
            }
        }
        subject_lines.push_back(sum_line);
    }
    file_in.close(); 
    return subject_lines;
}


static void print_nifti_to_csv_help(Tcl_Interp * interp){
    RESULT_LIT("help nifti_to_csv");
}

extern "C" int nifti_to_csv_command(ClientData clientData, Tcl_Interp * interp,
                                    int argc, const char * argv[]){
    
    string mask_filename;
    const char * phenotype_filename;
    string  trait_name;
    const char * base_filename;
    size_t max_voxels_per_file = 0;
    if(argc != 5 && argc != 6){
        print_nifti_to_csv_help(interp);
        return TCL_OK;
    }
    
    mask_filename = argv[1];
    trait_name = argv[2];
    phenotype_filename = argv[3];
    base_filename = argv[4];
    
    if(argc == 6){
        max_voxels_per_file = stoi(string(argv[5]));
    }
    //cout << "here\n";    
   // vector<voxel_struct> mask_vector;
     vector<int> mask_x;
     vector<int> mask_y;
     vector<int> mask_z;
    vector<string> names;
   // cout << "here 2 \n";
  //  try{
        get_mask_vector(mask_x, mask_y, mask_z, names,mask_filename);
	
/*	for(unsigned i = 0; i < mask_vector.size(); i++) {

		voxel_struct t = mask_vector[i];
	}*/
   /* }catch(...){
        RESULT_LIT("Error loading mask");
        return TCL_ERROR;
    }*/
  //  cout << "here 3\n";
    vector<string> subject_lines;
    vector<string> subject_files;
    string title_line;
    try{
        subject_lines =  get_subject_lines_and_subject_files(phenotype_filename, subject_files, trait_name, title_line);
    }catch(...){
        RESULT_LIT("Error loading phenotype file data");
        return TCL_ERROR;
    }
  //  cout << "here 4\n";
    float ** subject_data;
  //  try{
        subject_data = get_subject_data( mask_x,  mask_y,  mask_z, subject_files);
  //  }catch(std::bad_alloc& ba){
      //  RESULT_LIT("Bad allocation of memory encountered while loading subject data");
    //    return TCL_ERROR;
  //  }catch(out_of_bounds_exception & ex){
      //  RESULT_BUF(ex.what());
    //    return TCL_ERROR;
  //  }catch(trait_index_error& err){
      //  RESULT_BUF(err.what());
    //    return TCL_OK;
  //  }catch(...){
      //  RESULT_LIT("An error occurred in loading subject data");
    //    return TCL_ERROR;
  //  }
  //  try{
        write_data_out_file(subject_data,  names,  subject_lines, title_line, base_filename, max_voxels_per_file);
  //  }catch(...){
    //    RESULT_LIT("Error occurred while writing data out to file");
  //      return TCL_ERROR;
//    }
  //  cout << "here 5\n";
    for(size_t index= 0; index < subject_lines.size(); index++){
        delete [] subject_data[index];
    }
    delete [] subject_data;
  //  cout << "here 6\n";
   
    return TCL_OK;
    
}


//
//  annotate_gwas_cmd.cpp
//  
//
//  Created by Brian Donohue on 2/22/19.
//

#include <stdio.h>
#include "solar.h"
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;

static vector<string> split_line(string line){
    char c;
    string str;
    vector<string> output;
    for(unsigned i = 0 ; i < line.length(); i++){
        c = line[i];
        if(c != ','){
            str += c;
        }else{
            output.push_back(str);
            str.clear();
        }
    }
    output.push_back(str);
    return output;
}
class gwas_data
{
private:
    double threshold;
public:
    vector<string> snp_names;
    vector<string> gwas_lines;
    string field_line;
    
    gwas_data(const char * gwas_filename, const double l_threshold = 0.05){
        threshold = l_threshold;
        ifstream data_in(gwas_filename);
        getline(data_in, field_line);
        vector<string> parsed_line = split_line(field_line);
        size_t pvalue_index;
        for(unsigned i = 0; i < parsed_line.size(); i++){
            if(!StringCmp(parsed_line[i].c_str(), "pvalue", case_ins) || !StringCmp(parsed_line[i].c_str(), "p-value", case_ins)){
                pvalue_index = i;
                break;
            }
        }
        size_t snp_name_index;
        for(unsigned i = 0; i < parsed_line.size(); i++){
            if(!StringCmp(parsed_line[i].c_str(), "snp", case_ins) || !StringCmp(parsed_line[i].c_str(), "snp", case_ins)){
                snp_name_index = i;
                break;
            }
        }
        string line;
        while(getline(data_in, line)){
            parsed_line = split_line(line);
            double temp_pvalue = atof(parsed_line[pvalue_index].c_str());
            if(temp_pvalue <= threshold){
                snp_names.push_back(parsed_line[snp_name_index]);
                gwas_lines.push_back(line);

            }
        }
        data_in.close();
    }
    void write_annotation_data(const char * annotation_filename, const char * field_filename, const char * output_filename){
        ifstream field_stream(field_filename);
        string line;
        getline(field_stream, line);
	vector<string> field_list = split_line(line);
        field_stream.close();
        size_t comma_index = 0;
        while(line[comma_index] != ',') comma_index++;
        line = line.substr(comma_index, line.length() - comma_index - 1);
        ofstream output_stream(output_filename);
        
        output_stream << field_line << line << endl;
        ifstream annotation_stream(annotation_filename);
        
        vector<string>  temp_snp_names = snp_names;
        vector<string> temp_gwas_lines = gwas_lines;
        string snp_name;
        vector<string>::iterator iter;
        while(getline(annotation_stream, line) && temp_snp_names.size() != 0){
            comma_index = 0;
            while(line[comma_index] != ',') comma_index++;
            iter = find(temp_snp_names.begin(), temp_snp_names.end(), line.substr(0, comma_index));
            if(iter != temp_snp_names.end()){
                size_t vector_index = distance(temp_snp_names.begin(), iter);
                output_stream << temp_gwas_lines[vector_index]  << line.substr(comma_index, line.length() - comma_index - 1) << endl;
                temp_snp_names.erase(iter);
                temp_gwas_lines.erase(temp_gwas_lines.begin() + vector_index);
            }
        }
	for(unsigned index = 0; index < temp_snp_names.size(); index++){
		output_stream << temp_gwas_lines[index];
		for(size_t comma = 0; comma < field_list.size(); comma++){
			output_stream << ",N/A";
		}
		output_stream << endl;
	}
        output_stream.close();
        annotation_stream.close();
        
    }
    
};
static void print_help(Tcl_Interp * interp){
    Solar_Eval(interp, "help annotate_gwas");
    
}
extern "C" int annotationCmd(ClientData clientData, Tcl_Interp *interp,
                              int argc,const char *argv[]){
    const char * gwas_filename = 0;
    const char * annotation_filename = 0;
    const char * field_filename = 0;
    const char * output_filename = 0;
    double threshold = 0.05;
    for(unsigned arg = 1 ; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--i", case_ins) || !StringCmp(argv[arg], "-i", case_ins) ) && arg + 1 < argc){
            gwas_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--o", case_ins) || !StringCmp(argv[arg], "-o", case_ins) ) && arg + 1 < argc){
            output_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--a", case_ins) || !StringCmp(argv[arg], "-a", case_ins) ) && arg + 2 < argc){
            annotation_filename = argv[++arg];
            field_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--t", case_ins) || !StringCmp(argv[arg], "-t", case_ins) ) && arg + 1 < argc){
            threshold = atof(argv[++arg]);
        }else if(!StringCmp(argv[arg], "--help", case_ins) || !StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "help", case_ins)){
            print_help(interp);
            return TCL_OK;
        }else{
            RESULT_LIT("Invalid argument wsa entered");
            return TCL_ERROR;
        }
    }
    
    gwas_data annotator(gwas_filename, threshold);
    annotator.write_annotation_data(annotation_filename,  field_filename, output_filename);
    return TCL_OK;
    
}

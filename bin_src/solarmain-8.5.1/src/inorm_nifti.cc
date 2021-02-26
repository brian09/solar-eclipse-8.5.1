
#include "solar.h"
#include <string>
#include <vector>
#include <queue>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <chrono>
using namespace std;



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

void inorm(vector< pair<double, size_t> > data_in, float * output_data){
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
			output_data[current_id] = sum;
		}
		
		sum = 0.0;
		
	}
}


/*

extern "C" int inormNiftiCmd(ClientData clientData, Tcl_Interp *interp,
                 int argc,const char *argv[]){
auto start = std::chrono::high_resolution_clock::now();					 
	string output_filename;
	string mask_filename;
	string input_filename;
	for(int arg = 1; arg < argc ; arg++){
		if((!StringCmp(argv[arg], "--mask", case_ins) || !StringCmp(argv[arg], "-mask", case_ins))  &&
				arg + 1 != argc){
			mask_filename = argv[++arg];	
		}else if((!StringCmp(argv[arg], "--out", case_ins) || !StringCmp(argv[arg], "-out", case_ins) || !StringCmp(argv[arg], "-output", case_ins) ||
		!StringCmp(argv[arg], "--output", case_ins) || !StringCmp(argv[arg], "--o", case_ins) || !StringCmp(argv[arg], "-o", case_ins))  &&
				arg + 1 != argc){
			output_filename = argv[++arg];	
		}else if((!StringCmp(argv[arg], "--in", case_ins) || !StringCmp(argv[arg], "-in", case_ins) || !StringCmp(argv[arg], "-input", case_ins) ||
		!StringCmp(argv[arg], "--input", case_ins) || !StringCmp(argv[arg], "--i", case_ins) || !StringCmp(argv[arg], "-i", case_ins))  &&
				arg + 1 != argc){
			input_filename = argv[++arg];	
		}else{
			RESULT_LIT("An invalid argument was entered.");
		}		
	}
	
	RicVolumeSet  *  Input = new RicVolumeSet(input_filename);
	//Input = new RicVolumeSet(input_filename);
	cout << "Done loading input" << endl;
	SolarVolumeSet Output(Input->nx, Input->ny, Input->nz, Input->nvol, output_filename.c_str());
		cout << "Done creating output" << endl;

	if(mask_filename.length() == 0){
		#pragma omp for ordered schedule(static,1)
		for(size_t x = 0 ;x < Input->nx; x++) for(size_t y = 0 ; y < Input->ny; y++)for(size_t z = 0 ; z < Input->nz ; z++){		
			vector< pair<double, size_t> > data(Input->nvol);
			for(size_t vol = 0; vol < Input->nvol; vol++){
				data[vol] = pair<double, size_t>(Input->VolSet[vol].vox[x][y][z], vol);
			}
			float  output_data[Input->nvol];
			inorm(data, output_data);
			#pragma omp ordered 
			{
				Output.write_voxel_set(output_data, x, y, z);	
			}	
				
		}			
	}else{
		RicVolumeSet * Mask_file = 0;
		try{
			Mask_file = new RicVolumeSet(mask_filename);
		}catch(...){
			RESULT_LIT("Error reading mask volume set");
			return TCL_ERROR;
		}
		#pragma omp parallel for collapse(3)
		for(size_t x = 0 ;x < Input->nx; x++) for(size_t y = 0 ; y < Input->ny; y++)for(size_t z = 0 ; z < Input->nz ; z++){	
			if(Mask_file->VolSet[0].vox[x][y][z] != 0) {
				vector< pair<double, size_t> > data(Input->nvol);
				for(size_t vol = 0; vol < Input->nvol; vol++){
					data[vol] = pair<double, size_t>(Input->VolSet[vol].vox[x][y][z], vol);
				}
				float  output_data[Input->nvol];
				inorm(data, output_data);
				Output.write_voxel_set(output_data, x, y, z);		
					
			}else{
				float  output_voxel_set[Input->nvol];
				memset(output_voxel_set, 0, sizeof(float)*Input->nvol);
				Output.write_voxel_set(output_voxel_set, x, y, z);
				
			}
		}			
		
		delete Mask_file;
		
	}
	
	delete Input;

	auto elapsed = std::chrono::high_resolution_clock::now() - start;
	
  auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);		
    cout << seconds.count() << endl;
	return TCL_OK;	 
}

*/

#include <string>
#include <algorithm>
#include <iostream>
#include <cstring>
#include "gpu-exception.h"
#include "gpu-selection.h"
using namespace std;

static void print_gpu_info(int gpu_id){
	cudaDeviceProp 	prop;
	cudaErrorCheck(cudaGetDeviceProperties(&prop, gpu_id));
	cout << "GPU ID: " << gpu_id << " Name: " << prop.name << endl;
	cout << "Device Architecture: " << prop.major << "." << prop.minor << endl;
	cout << "Total Global Memory: " <<  0.0009765625*prop.totalGlobalMem << " kilobytes\n";
	cout << "Clockrate: " << prop.clockRate << " kilohertz\n";
	cout << "Number of Multiprocessors: " << prop.multiProcessorCount << endl;
	
}

static vector<int> get_usable_device_ids(const int device_count){
	
	vector<int> device_list;
	cudaDeviceProp prop;
	for(int device_id = 0; device_id < device_count ; device_id++){
		cudaErrorCheck(cudaGetDeviceProperties(&prop, device_id));
		if(prop.major >= 3 && prop.computeMode == 0){
			device_list.push_back(device_id);
		}
			 
		
	}
	
	return device_list;
}
static vector<int> device_select(){
	int device_count;
	cudaErrorCheck(cudaGetDeviceCount(&device_count));
	vector<int> usable_device_list = get_usable_device_ids(device_count);
	device_count = usable_device_list.size();
	cout << "Listing NVIDIA GPUs pf architecture greather than equal to 3.0 and compute mode 0\n\n";
	for(int device_id = 0; device_id < device_count ; device_id++){
		print_gpu_info(usable_device_list[device_id]);
		cout << endl;
	}
	vector<int> device_list;
	
	bool done_selection = false;
	do{
		device_list.clear();
		for(int device_id = 0; device_id < device_count; device_id++){
			bool end_selection = false;
			print_gpu_info(usable_device_list[device_id]);
			while(!end_selection){
				cout << "Would you like to use this device? (y/n)";
				string selection;
				cin >> selection;
				cout << endl;
				std::transform(selection.begin(), selection.end(),selection.begin(), ::toupper);
				if(selection == "YES" || selection == "Y"){
					device_list.push_back(usable_device_list[device_id]);
					end_selection = true;
				}else if(selection == "NO" || selection == "NO"){
					end_selection = true;
					
				}
			}
		}
		
		cout << "The following devices were selected: \n";
		for(int device_id = 0; device_id < device_list.size() ; device_id++){
			print_gpu_info(device_list[device_id]);
			cout << endl;
		}		
		bool end_selection = false; 
		while(!end_selection){
			cout << "Would you like to proceed with these devices? (y/n)";
			string selection;
			cin >> selection;
			cout << endl;
			std::transform(selection.begin(), selection.end(),selection.begin(), ::toupper);
			if(selection == "YES" || selection == "Y"){
				done_selection = true;
				end_selection = true; 
			}else if(selection == "NO" || selection == "NO"){
				end_selection = true;
			}
		}		
		
	}while(!done_selection);
	
	
	
	return device_list;
}

static vector<int> convert_string_to_gpu_list(const char * gpu_list_cstring){
	string gpu_list_string = gpu_list_cstring;
	string str;
	vector<int> device_list;
	int id;
	for(size_t idx = 0; idx < gpu_list_string.length(); idx++){
		char c = gpu_list_string[idx];
		if(c != ','){
			str += c;
		}else{
			id = stoi(str);
			device_list.push_back(id);
			str.clear();
		}
		
	}
    id = stoi(str);
	device_list.push_back(id);
	int device_count;
	cudaErrorCheck(cudaGetDeviceCount(&device_count));
	vector<int> usable_devices = get_usable_device_ids(device_count);
	for(size_t i = 0; i < device_list.size(); i++){
		id = device_list[i];
		vector<int>::iterator find_iter = find(usable_devices.begin(), usable_devices.end(), id);
		if(find_iter == usable_devices.end()){
			cout << "Device ID " << id << " is not a usable NVIDIA GPU\n";
			device_list.clear();
			return device_list;
		}
		
	}

	return device_list;
}

void print_usable_devices(){
	
	int device_count;
	cudaErrorCheck(cudaGetDeviceCount(&device_count));
	vector<int> usable_device_list = get_usable_device_ids(device_count);
	device_count = usable_device_list.size();
	cout << "Listing NVIDIA GPUs pf architecture greather than equal to 3.0 and compute mode 0\n\n";
	for(int device_id = 0; device_id < device_count ; device_id++){
		print_gpu_info(usable_device_list[device_id]);
		cout << endl;
	}
	
}
vector<int> get_gpu_id_list(const char * gpu_list, bool select_all){
	if(gpu_list == 0 && !select_all){
		return device_select();
		
	}else if(gpu_list != 0){
		return convert_string_to_gpu_list(gpu_list);
		
	}else{
		int device_count;
		cudaErrorCheck(cudaGetDeviceCount(&device_count));
		return get_usable_device_ids(device_count);
	}
	
	
}

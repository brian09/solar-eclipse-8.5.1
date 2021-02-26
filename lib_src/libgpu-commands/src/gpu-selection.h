//
//  gpu-selection.h
//  
//
//  Created by Brian Donohue on 8/23/18.
//

#ifndef gpu_selection_h
#define gpu_selection_h
#include <vector>
void print_usable_devices();
std::vector<int> get_gpu_id_list(const char * gpu_list, bool select_all);

#endif /* gpu-selection_h */

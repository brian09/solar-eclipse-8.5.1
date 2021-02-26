
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "solar.h"
//#include <cstdlib>
#include <time.h>
#include "RicVolumeSet.h"
#include <sys/stat.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

// Please down load the necessary libraries from http://www.nitrc.org/projects/brainvisa_ext for compilation of this code


using namespace std;
static string float2string(float tmp);

//main storage class for voxel-wise data
class table_column {
public:
int N;
float *array;
string label;
int x,y,z;

table_column()
{
	N=0;
	x=y=z=0;
}

table_column(int n_col)
{
 N=n_col;
 array=new float[n_col];

}

void init(int n_col)
{
 N=n_col;
 array=new float[n_col];

}


 ~table_column(){delete [] array;}
};


// This function will partition the file name into sub-components based on "_" separator
static void split_the_name(string name, string &mr_id, string &a_id, string &age, string &sex)
{
int location[]={0,0,0,0,0,0,0};
int instance=0;

for (int i=name.find("_", 0); i!=string::npos;i=name.find("_",i))
{

     instance++;
        location[instance]=i;
        i++;
        if (instance > 5) {cerr<< " Too many variables in the name "<<name<<endl; exit(1);}

}

mr_id=name.substr(0, location[1]);
a_id=name.substr( location[1]+1, location[2]- location[1]-1);

age= name.substr( location[2]+1, location[3]-location[2]-1);

sex= name.substr( location[3]+1, location[4]-location[3]-1);

}

static long n_column(RicVolumeSet *vol, int threshold)
{
long  N_col=0;;
for (int i=0; i< vol->VolSet[0].num_cols();i ++)
        for(int j=0; j<vol->VolSet[0].num_rows();j++)
                for (int k=0;k<vol->VolSet[0].num_slices();k++)
                           if (vol->VolSet[0].vox[i][j][k]>threshold) {N_col++;   vol->VolSet[0].vox[i][j][k]=N_col;}
                                                        else   vol->VolSet[0].vox[i][j][k]=0;
return N_col;
}


static bool file_exists(const char * filename){
    
    ifstream input(filename);
    
    bool does_exist = input.good();
    
    input.close();
    
    return does_exist;
    
}

static void fill_table(table_column* table, int n_file, RicVolumeSet *index, RicVolumeSet *data)
{

long ref_n=0;

for (int i=1; i<data->VolSet[0].nx-1 ; i++){
        string x_str = to_string(i);
	for(int j=1; j<data->VolSet[0].ny-1;j++){
                string y_str = to_string(j);
		for (int k=1;k<data->VolSet[0].nz-1;k++){
               ref_n= index->VolSet[0].vox[i][j][k];

                if (ref_n > 0) {
                	// we are referencing the voxels by their sequential number in the mask
                	//cerr<<ref_n<<" "<< i<<" "<<j<<" "<<k<<" "<<data->VolSet[0].vox[i][j][k]<<endl;
                                                    table[ref_n-1].array[n_file]= data->VolSet[0].vox[i][j][k];
                                                    table[ref_n-1].x=i;
                                                    table[ref_n-1].y=j;
                                                    table[ref_n-1].z=k;
                                                    table[ref_n-1].label="VOXEL_"+ x_str +"_"+ y_str +"_"+ to_string(k);

						}
                      }
	}
}
}

static string float2string(float tmp)
{
 ostringstream oss;
oss<<tmp;
return oss.str();

}

static void usage(Tcl_Interp * interp)
{
	
	Solar_Eval(interp, "help nifti_to_csv");



}

static vector<string> get_file_list(const char * file_list_name){
    
    ifstream input(file_list_name);
    vector<string> file_list;
   
    if(input.is_open() == false){
        return file_list;
    }
    
    
    string filename;
    
    while(getline(input, filename)){
        file_list.push_back(filename);
    }
    
    
    input.close();
    
    return file_list;
    
    
}

extern "C" int niftiToSolar(ClientData clientData, Tcl_Interp *interp,
                                         int argc,const char *argv[])
{
/*
 * This is a simple program to convert voxel-wise data to csv files for solar analysis.
 * This is a bandaid until the nifti reading module works reliably in solar which might take a while
 * To save on memory, the voxels in the mask file are given sequential index
 * An array of columns is then created to store the data, the size of the array is equal to the number of voxels in the mask
 * Then I create the array of rows that is equal to the number of subjects in the data set e.g. the number of files
 * Finally, the tables with the predeterminied size are saved together with headers.
 *
 *
 *
 */


    if (argc != 4){
        usage(interp);
        return TCL_ERROR;
    }
RicVolumeSet mask(argv[1]);
vector<string> file_list = get_file_list(argv[3]);
    if(file_list.size() == 0){
        RESULT_LIT("File list is empty or does not exist.");
        return TCL_ERROR;
    }
    
int X,Y,Z;
X=mask.VolSet[0].num_cols();
Y=mask.VolSet[0].num_rows();
Z=mask.VolSet[0].num_slices();
int N_vols=file_list.size();  // total number of files
int N_col_per_table= atoi(argv[2]);
    
    if(N_col_per_table <= 0){
        RESULT_LIT("Invalid number of columns per table was entered.");
        return TCL_ERROR;
    }

int cur_file=0;
int cur_col=0;
string names[file_list.size()];
string mri_id[file_list.size()];
string a_id[file_list.size()];
string age[file_list.size()];
string sex[file_list.size()];
// The main loop that will go through the files.
cerr<<"Will load "<<file_list.size()<<" volumes"<<endl;
int N_col_total=  n_column (&mask, 0.5);
cerr<<"There are a total of N="<<N_col_total<<" in the mask"<<endl;
// initialize the storage space for the table
//table_column columns[N_col_total];
table_column * columns = new table_column[N_col_total];
for (int i=0; i< N_col_total;i++){ columns[i].init(N_vols); for (int j=0;j<N_vols;j++) columns[i].array[j]=0;}


// This loop will fill the table with values
for (int i=0;i< N_vols;i++)
{
	string tmp= file_list[i];
	names[i]=tmp.substr(0, tmp.find(".nii"));

	// lets partition the file name into the subcomponents for the table
	split_the_name( names[i], mri_id[i], a_id[i], age[i], sex[i]);
	cerr<<"Loading "<<mri_id[i]<<" "<< a_id[i]<<" "<<age[i]<<" "<<sex[i]<<endl;
    if(file_exists(file_list[i].c_str()) == false){
        string error_msg = "Could not find file " + file_list[i] + ".";
        RESULT_LIT(error_msg.c_str());
        return TCL_ERROR;
    }
	RicVolumeSet loader(file_list[i]);// Declared statically to avoid memory leaks
	//cerr<<"loaded";
	fill_table(  columns, i, &mask, &loader);
}




// This loop will print the the table into the csv/header files.
for (int i=0; i<N_col_total;i+= N_col_per_table)
{
        string f_name="out_"+float2string((i+1)/N_col_per_table)+".csv";
        string fh_name="out_"+float2string((i+1)/N_col_per_table)+".header";
        ofstream out_file(f_name.c_str());
        ofstream header_file(fh_name.c_str());
        cerr<<"Writing "<<f_name<<endl;

        for (int n=0;n< N_vols+1;n++)
        {
                if (n==0)
                {out_file<<"mrid,"<<"id,"<<"age,"<<"sex"; }

                else   out_file<<mri_id[n-1]<<","<<a_id[n-1]<<","<<age[n-1]<<","<<sex[n-1];
                for (int j=0;j< N_col_per_table;j++)
                        {
                                if (i+j< N_col_total)
                                {
                                        if (n==0) {out_file<< "," << columns[i+j].label;
                                                                        header_file<<columns[i+j].label<<"\t";}
                                        else       out_file<< "," << columns[i+j].array[n-1];
                                }
                        }
 out_file<<endl;
}
out_file.close();
header_file.close();
}
delete [] columns;


return TCL_OK;

}


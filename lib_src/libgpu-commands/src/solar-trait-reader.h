//
//  solar-trait-reader.h
//  
//
//  Created by Brian Donohue on 8/19/18.
//

#ifndef solar_trait_Reader_h
#define solar_trait_Reader_h

#include <stdio.h>
#include "solar.h"
#include <string>
#include <vector>
#include <exception>
class Solar_Trait_Reader_Exception : virtual public std::exception{
protected:
    const char * error_message;
public:
    Solar_Trait_Reader_Exception(const  char * msg): error_message(msg){

    }
    virtual ~Solar_Trait_Reader_Exception() throw() {}
    virtual const char * what() const throw()
	{
		return error_message;
	}
};
void load_phi2_matrix(Tcl_Interp * interp);
class Eigen_Data{
private:
    void calculate_eigenvectors_and_eigenvalues (double * phi2, double * eigenvectors ,int n);
    std::vector<size_t> ibdids;
    std::vector<std::string> ids;
    std::vector<std::string> trait_names;
    double * eigenvalues;
    double * eigenvectors_transposed;
    double * phenotype_buffer;
    std::vector<int> index_map;
    size_t n_phenotypes;
    size_t n_subjects;
public:
    inline double * get_phenotype_column(const size_t index) const {return (index < n_phenotypes) ? &phenotype_buffer[index*n_subjects] : 0;}
    inline std::string get_trait_name(const size_t index) { return  (index < n_phenotypes) ? trait_names[index] : std::string("");}
    inline double * get_eigenvalues() const {return eigenvalues;};
    inline double * get_eigenvectors_transposed() const {return eigenvectors_transposed;};
    inline double * get_phenotype_buffer() const {return phenotype_buffer;};
    inline std::vector<std::string> get_trait_names() {return trait_names;};
    inline std::vector<std::string> get_ids() {return ids;};
    inline  size_t get_n_subjects() {return n_subjects;};
    inline  size_t get_n_phenotypes() {return n_phenotypes;};
   
    Eigen_Data(std::vector<std::string>, const char *, const size_t);
    Eigen_Data(std::vector<std::string>, std::vector<std::string>, const size_t);
    ~Eigen_Data();
    void set_phenotype_column(const size_t column_index, std::string name,  double *  input_buffer);
};

class Solar_Trait_Reader{
private:
    void Sort_ID_Sets(const double * const , std::vector<std::string> , std::vector<std::string> , \
                      std::vector<size_t> & , std::vector<size_t> & , std::vector< std::vector <std::string> > & ,\
                      const size_t , const size_t );
    void Read_Phenotype_Data_Thread_Launch(SolarFile * const , double * const , std::vector<std::string> ,const size_t , const size_t , const size_t );
    void Read_Phenotype_Data_Thread_Launch(SolarFile * const , double * const ,std::vector<std::string> , std::vector<std::string> , std::vector<std::string> & ,\
								 const unsigned , const unsigned , const unsigned );
    size_t n_sets;
    size_t n_phenotypes;
    Eigen_Data ** eigen_data;
public:
    Solar_Trait_Reader(const char * , std::vector<std::string>, std::vector<std::string> );
    Solar_Trait_Reader(const char *, const char *, std::vector<std::string>);
    ~Solar_Trait_Reader();
    inline size_t get_n_sets() { return n_sets; };
    inline size_t get_n_phenotypes() {return n_phenotypes; };
    inline Eigen_Data * get_eigen_data_set(const size_t index) const { return  (index < n_sets) ? eigen_data[index] : 0;}
    
};
#endif /* Solar_Trait_Reader_h */

#ifndef __PUBLIC_H__
#define __PUBLIC_H__

//MJJ 2/1/09
#include <cstdlib>

#include <math.h> 
#include "cstring.h"
#include <iomanip>
#include "arrays.h"
                   

using namespace std; 

typedef double my_float; 
                        
enum Mut_Type {MICROSAT, RFLP, DNA, SNP, STANDARD};

// A routine for generating good random numbers
extern double ran3(long *idum);

extern float gammln(float xx);

// A routine to generate poisson deviates with mean xm 
extern float poidev(float xm, long *idum);

extern void quicksort(int l, int r, double* list);

// A routine to initialize the probabilities for the Gamma distributed mutation rates.
void init_gamma_weights(double* p, int n_loci, double alph);
//Discrete gamma case, using Ziheng Yang's routines
void init_gamma_weights(double* p, int n_loci, double alph, const int num_rates);


//returns gamma deviates
extern double gamma_dev(double a);
extern double igamma_dev(int ia);
typedef MY_TArrayAsVector<my_string> StringArray;

//returns gaussian deviates
double gasdev(long *idum);

//Arlequin output routines
//------------------------------------------------------------------------------
extern void write_Arlequin_header(const int& num_samples,
											  const Mut_Type& data_type, ostream& os,
                                   const int& genot_data);
extern void write_Arlequin_sample_header(const int& num_pop, const int& size,
														ostream & os, const int& genot_data);
extern void write_Arlequin_sample_footer(ostream & os);



//PHASE output routines

//------------------------------------------------------------------------------

extern void write_PHASE_Header(ostream & os, const int& numLoci, const int& numInds,

                                             const Mut_Type& data_type);

                               

//------------------------------------------------------------------------------


//PAUP* output routines

//------------------------------------------------------------------------------
//Recombination
void 
write_output_coalTimes(ostream&    os, 
                       const int&   num_stat_displayed,
                       char**       list_name_stat_time_computed,
                       double**     list_stat_time_computed,
                       const int&   num_ind_loci,
                       const int&   num_linked_loci); 

void 
write_output_RecomStats(ostream&    os, 
                       const int&   num_stat_displayed,
                       char**       list_name_stat_time_computed,
                       double**     list_stat_time_computed,
                       const int&   num_ind_loci,
                       const int&   num_linked_loci);

//------------------------------------------------------------------------------
enum tree_output_type {GENERATIONS=0, MUT_RATE=1, NUM_MUT=2};

enum Mut_Model {K80_GAMMA, K80_NOGAMMA, HKY85_GAMMA, HKY85_NOGAMMA/*, JK, TN93, HK85*/};

extern void write_PAUP_header(ostream& os, const char* logFile);
extern void write_PAUP_block(ostream & os, const Mut_Model mut_mod);
extern void write_PAUP_footer(ostream & os);

extern void write_PAUP_data_header(const int ntax, const int nchar,
												Mut_Type datatype, ostream & os);
extern void write_PAUP_end(ostream & os);

extern void write_PAUP_end_matrix(ostream & os);

extern void write_PAUP_tree_header(ostream & os, const int num_rep);

extern void write_PAUP_replicate_number(const int n, ostream & os);
extern void write_PAUP_tree_name(ostream & os, const int num_rep, 
                                               const tree_output_type& tree_type);

//Different tree between nodes
extern void write_PAUP_tree_name_perLoc(ostream & os, const int num_rep, const int& curLocus, 
                                                                         const tree_output_type& tree_type);


extern void write_PAUP_trees_section_header(ostream & os);

//Ziheng Yang routines for generating discrete gamma distributions
//------------------------------------------------------------------------------
extern int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf, double space[]);
extern double LnGamma (double alpha);
extern double DFGamma(double x, double alpha, double beta);
extern double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
extern double PointChi2 (double prob, double v);
extern double PointNormal (double prob);
extern int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);
extern double rndgamma (double s);
extern int MultiNomial (int n, int ncat, double prob[], int nobs[], double space[]);
extern double rndu (void);
extern int xtoy (double x[], double y[], int n);
extern int abyx (double a, double x[], int n);
extern int geometric(const double &p);
//------------------------------------------------------------------------------

#endif
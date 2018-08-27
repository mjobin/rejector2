#ifndef __DEME_H__
#define __DEME_H__

//#include <classlib\arrays.h>
#include <cstdlib>
#include "arrays.h"
#include "public.h"
#include "migrmat.h"
#include "genealgy.h"
#include "locus.h"
#include "cond_var.h"
#include "chain.h"

using namespace std;

const long MAX_GENERATIONS=10000000l;
const long MIN_SIZE=1;
const long MAX_SIZE=2000000000L;

typedef MY_TArrayAsVector<TLocus>	TLocusArray;
typedef TISimpleChain<TNode>	NodeList_chain;
typedef MY_TArrayAsVector<int>          TDataTypeArray; //Loro_24_03_03 : Array of data type per locus

//------------------------------------------------------------------------------
//A very simple class to store the index of an ancestral node and the pointer of its descendent node
class TDescNode {
   public:
      int ancestorIndex;
      TNode * pDescNode;
      TDescNode() {ancestorIndex=-1; pDescNode=NULL;}
      TDescNode(const int& a, TNode* pN) {ancestorIndex=a; pDescNode=pN;}
      TDescNode(const TDescNode& DN) {
         if (this!=&DN) *this=DN;
      }
      TDescNode& operator=(const TDescNode& DN) {
         ancestorIndex=DN.ancestorIndex;
         pDescNode=DN.pDescNode;
         return *this;
      }
      int operator==(const TDescNode& DN) {return ancestorIndex==DN.ancestorIndex;}
      int operator<(const TDescNode& DN) {return ancestorIndex<DN.ancestorIndex;}
};

typedef MY_TArrayAsVector<TDescNode> DescNodeVect;
//------------------------------------------------------------------------------

class TDeme {

	friend ostream& operator<<(ostream& os , const TDeme& D);

	private :
   	static TMigrationMatrix* Migrations; //A pointer to a migration matrix
      double _growth_rate;
      long  _deme_size;
      long  _min_deme_size, _max_deme_size;
      long  _sample_size;
      long  _lineages;
      double _proba_coal;
      long _lidum;    			//dummy variable for drawing random numbers



      
      //TListNode NodeList;    	//A list of the nodes currently present in the deme
      								//(as pointers to nodes in the Genealogy)
                              //Remember not to delete the elements of the NodeList
                              //as they belong to another tree
      //TListNode OriginalNodes;//A list of the nodes originally present at generation zero in the deme
      int _id;

      bool _is_migration;
      long _polym_sites;
      TLocusArray _loci;

	public :

	//MJJ
   //	int oldestLineage;
	
      /**********************************************************************
         
            For recombination : don t use  NodeList and OriginalNodes
         
       **********************************************************************/
      //TDeme() : NodeList(), OriginalNodes(), _loci(0,10)  {
      //   initialize();
      //};
      //TDeme(const TDeme& D)  : NodeList(), OriginalNodes(), _loci(0,10) {
      //	if (this!=&D) *this=D;
      //};
   	//TDeme(const long& ds, const long& ss, const float& gr) : NodeList(),
      //                                                         OriginalNodes(),
      //                                                         _loci(0,10)  {
      //   _growth_rate=gr;
      //   _sample_size=ss;
      //   _deme_size=ds;
      //   _lineages=_sample_size;
      //   _lidum=1L;
      //   _is_migration=true; 
      //   Migrations=NULL;
      // 	_polym_sites=0L;
      //};
      TDeme() : ChainNodeList(false), ChainOriginalNodes(false), _loci(0,10)  {
         initialize();
             
           
      };
      TDeme(const TDeme& D)  : ChainNodeList(false), ChainOriginalNodes(false), _loci(0,10) {
      	if (this!=&D) *this=D;
       
      };
   	TDeme(const long& ds, const long& ss, const float& gr) :  ChainNodeList(false),
      																			 ChainOriginalNodes(false),
                                                                _loci(0,10)  {
         _growth_rate=gr;
         _sample_size=ss;
         _deme_size=ds;
         _lineages=_sample_size;
         _lidum=1L;
         _is_migration=true; 
         Migrations=NULL;
       	_polym_sites=0L;
 
      };
      //virtual ~TDeme() {
         //For Debug
      //   cout << "Enter TDeme destructor\n";
      //}
      ~TDeme() {

      }
      void set_growth_rate(const float& gr) {_growth_rate=gr;}
      void set_sample_size(const long& size) {
      	_sample_size=size;
         _lineages=size;
      }
      void set_size_bounds(const long& min, const long& max) {
			_min_deme_size=min;
			_max_deme_size=max;
      }
      void set_deme_size(const long& size) {
      	_deme_size=size;
   		if (size) _proba_coal=1.0/_deme_size;
         else _proba_coal=0.0;
      }
      //Recombination
      long get_deme_size() { return _deme_size;}
      
      void set_id(const int& id) {_id=id;}
      int  get_id() const {return _id;};
      void initialize() {
      	_growth_rate=0;  //Stationary population
         _sample_size=0;
         _deme_size=1;
         _min_deme_size=MIN_SIZE;
			_max_deme_size=MAX_SIZE;
         _lineages=_sample_size;
         _lidum=1L;        
         _is_migration=true;
         Migrations=NULL; 
       	_polym_sites=0L;
      }
      long lineages() const {return _lineages;}
      void set_migration_matrix(TMigrationMatrix* const MM) {
      	Migrations=MM;
      }
      TDeme& operator=(const TDeme& D);
      int operator==(const TDeme& D) const {
      	if (_id==D._id) return 1;
         return 0;
      };
      int operator<(const TDeme& D) const {
      	return _id<D._id;
      }
      int copy_node_pointers(TTree* T, const int& curr_lin,
      							  const long& max_lin);
      //Implementation of a single migration event towards another deme (sink)
      /* //Old version
      int migrate(TDeme& sink) {
			//Select a migrating node at random
   		int randNode= (int) (ran3(&_lidum)*_lineages);
   		//Add selected node to sink and then increments sink's number of lineages
   		sink.NodeList[sink._lineages++]=NodeList[randNode];

   		Decrements last node's position, and replace selected node by last node
   		NodeList[randNode]=NodeList[--_lineages];
   		return 1;
		};
      */
      //New version with recombination
      int migrate(TDeme& sink) {
			//Select a migrating node at random
   		int randNode= (int) (ran3(&_lidum)*_lineages);
         TNode* outgoingNode=ChainNodeList.remove(randNode);
         --_lineages;
         ++sink._lineages;
         sink.ChainNodeList.Add(outgoingNode);
   		return 1;
		};

      //Check and implements migrations
      int send_migrants(TDeme& sink, const int& i);

      //Implements the coalescence process within population
      int coalesce_lineages(const long& time, TTree& GeneTree, const int& pos);

      
      float growth() {return _growth_rate;}

      //Implements an exponential growth of population for positive _growth_rate
      //and and exponential decline for negative _growth_rate values
      const long& exponential_resize() {
      	_deme_size= (long) (_deme_size*exp(_growth_rate));
         if (_deme_size<_min_deme_size) {
         	_deme_size=_min_deme_size;
				_growth_rate=0.0;
         }
         if (_deme_size>_max_deme_size) {
         	_deme_size=_max_deme_size;
				_growth_rate=0.0;
         }
   		if (_deme_size) _proba_coal=1.0/_deme_size;
         else _proba_coal=0.0;
         return _deme_size;
      }

      //Implements a resizing of the deme size according to a given factor
      void linear_resize(const float& factor) {
      	_deme_size= (long) (_deme_size*factor);
   		if (_deme_size) _proba_coal=1.0/_deme_size;
         else _proba_coal=0.0;
      }
      void reset();
      long sample_size() {return _sample_size;}
		int  check_migrations();
      bool is_migration() {return _is_migration;}
      long deme_size() {return _deme_size;}
      void print_haplotypes_for_Arlequin(ostream & os, const int& n,
      												  const Mut_Type& data_type);
      void print_haplotypes_for_PAUP(ostream & os, const Mut_Type& data_type);
      
		//Loro_19_06_01
		void print_haplotypes_to_locus(const Mut_Type& data_type);
      
      //Recombination: New function to deal with mix of different data types
      void print_haplotypes_to_locus(const TDataTypeArray* data_type_perLoc);
      
		void flushLoci() {_loci.Flush();};
      void print_loci_for_Arlequin(ostream & os, const int& n,
      										const Mut_Type& data_type,
												const int& genot_data);
      void print_loci_for_Phase(ostream & os, const int& n, const Mut_Type& data_type);
      //Recombination: Mix of different loci
      void print_loci_for_Phase(ostream & os, const int& deme, const int& num_linked_loci, 
                                                           const TDataTypeArray* data_type_perLoc);
      void print_loci_for_Haplotyper(ostream & os, const int& n);

      /*************************************************************************   

            New procedures necessary for recombinations
      
      *************************************************************************/
      
      int copy_node_pointers(TNode* N,  TTree * T);
   
      NodeList_chain ChainNodeList;   // chained list of Nodes 
      NodeList_chain ChainOriginalNodes;   // chained list of Nodes
      
      //long recLineages;                     // Num of lineages created by recombination
                                                 
      int singleCoalescentEvent(const long& time, TTree& GeneTree, TDemeCollection* Ademe );
      int multipleCoalescentEvents(const long& time, TTree& GeneTree, TDemeCollection* Ademe );
      int recombination_per_lineage(const long& time, TTree& GeneTree, TDemeCollection* myDemes);
      int One_recombination_per_generation(const long& time, TTree& GeneTree, TDemeCollection* myDemes);
      void resetAllNodeFlags(const int& numLinkedLoci);
      void resetSequenceSize();

      //void bidon(){ int test; }
};

extern ostream& operator<<(ostream& os , const TDeme& D);


//------------------------------------------------------------------------------
// A class specifying the conditions of a demographic historical change
// sometimes in the past, like the fusion of two demes, or a change in
// deme size

class THistoricalEvent {
	friend ostream& operator<<(ostream& os, const THistoricalEvent& HE);
	friend istream& operator>>(istream& is, THistoricalEvent& HE);
	public:
		long time;     			//time in the past at which the event occurs (in generations)
      int source;    			//deme which sends of migrants
      int sink;      			//deme which receives migrants
      float migrants; 			//Proportion of migrants sent to sink from source
      float new_deme_size; 	//New relative deme size for SOURCE deme
      float new_growth_rate; 	//New growth rate for SOURCE deme
      int MigMat;					//New migration matrix number
      THistoricalEvent() {};
      ~THistoricalEvent() {};
      int operator==(const THistoricalEvent& HE) {
      	if (	source==HE.source 					&&
         		sink==HE.sink							&&
               time==HE.time							&&
               migrants==HE.migrants 				&&
               new_deme_size==HE.new_deme_size 	&&
               new_growth_rate==HE.new_growth_rate &&
               MigMat==HE.MigMat
         	) return 1;
         return 0;
      }
      int operator<(const THistoricalEvent& HE) const {
      	return time<HE.time;
      }
};

extern ostream& operator<<(ostream& os, const THistoricalEvent& HE); 
extern istream& operator>>(istream& is, THistoricalEvent& HE);



//MJJ

class TSNPEvent {

	public:
	int chrom;
	int locus;
	int time;
	int deme;
	int eventnum;


};

//------------------------------------------------------------------------------
typedef MY_TArrayAsVector<TDeme> 				TDemeArray;
typedef MY_TSArrayAsVector<THistoricalEvent> TEventArray;
typedef MY_TArrayAsVector<long> 					TDemeSizeArray;
typedef MY_TArrayAsVector<float>             TGrowthArray;
typedef MY_TArrayAsVector<long>					TSampSizeArray;
typedef MY_TArrayAsVector<TMigrationMatrix>  TMigrMatArray; //Loro_16_9_99
/* TODO -olaurent -crecombinaison : Need to use the following arrays to be put into TDemeCollection */ 
typedef MY_TArrayAsVector<double>            TRecRateArray; //Loro_24_03_03 : Array of recombination rate per locus
typedef MY_TArrayAsVector<int>               TNumLociArray; //Loro_24_03_03 : Array of number of loci per blocks
typedef MY_TArrayAsVector<double>            TMutRateArray; //Loro_24_03_03 : Array of mutation rate per locus
//typedef MY_TArrayAsVector<Mut_Type>          TDataTypeArray; //Loro_24_03_03 : Array of data type per locus
typedef MY_TArrayAsVector<int>               TDataTypeArray; //Loro_24_03_03 : Array of data type per locus


//MJJ
typedef MY_TArrayAsVector<TSNPEvent> 				TSNPEventArray;



class TDemeCollection {
      friend ostream& operator<<(ostream& os,const TDemeCollection& DC);
   private :
       
      TDeme	*Demes;
      TEventArray 		Events;
      TMigrMatArray 		*MigrMatArray;
      TMigrationMatrix 	*CoalTimes; //A matrix of mean pairwise coalescence times
      TMigrationMatrix 	*PairDiff; //A matrix of mean number of pairwise differences
      TMigrationMatrix 	*MinCoalTimes; //A matrix of minimum pairwise coalescence times
      bool 					_is_migration;
      int 					_num_demes;
      long 					_polym_sites;
      long 					_tot_mut;
      int 					_cur_mig_mat;

      /*************************************************************************
                            New variables for recombination
      **************************************************************************/
      int*   recombHits;               // number of recombination hits between adjacent loci
      int*   effective_recombHits;     // recombination hits between adjacent loci 
                                       //affecting the gametic phase: before MRCA
                                    
      int num_linked_loci;          

      
      int num_linkage_blocks;
      int num_indep_loci;
      int sum_sample_size;           
      double 	rec_rate;
      //Loro_01_03_04 slight modifs to avoid name confusion...
      const TRecRateArray* _rec_rate_perLoc;
      const TNumLociArray* _num_linked_loci_perBlock;
      const TMutRateArray* _mut_rate_perLoc;
      const TMutRateArray* _geom_param_perLoc;
      const TNumLociArray* _rangeConstraintPerLoc;
      const TMutRateArray* _transRatePerLoc;


	//MJJ 
	TSNPEventArray SNPEvents;

	
                        
      const double* freq_SNP_min;
      //************************************************************************
   public :

      
   	TDemeCollection() :
      		Events(1,10),SNPEvents(1,10),
            GeneTree() {
      	CoalTimes=NULL;
         MinCoalTimes=NULL;
         PairDiff=NULL;
         Demes=NULL;
         _num_demes=0; 
       	_polym_sites=0L;
       	_tot_mut=0L;
         _cur_mig_mat=0;
		 

		 
      };
      ~TDemeCollection() {
         if (CoalTimes)    delete CoalTimes;
      	if (PairDiff)     delete PairDiff;
      	if (MinCoalTimes) delete MinCoalTimes;
         if (Demes)        delete[] Demes;

         #ifdef _COMPUTE_MOMENTS_
            //For recombination
            if (recombHits)   delete[] recombHits;  
            if (effective_recombHits)   delete[] effective_recombHits;
         #endif
      }
      TDemeCollection& operator=(const TDemeCollection& DC);
   	long count_lineages();
      int reajust_deme_size(const int i, const long& ds);
      int create_demes(const TMigrMatArray* const MA);
      int create_lineages(const TSampSizeArray& SS);
      int initialize_growth_rates(const TGrowthArray& GR);
      int initialize_deme_size(const TDemeSizeArray& DS);
      int initialize_events(const TEventArray& HE);
      int initialize_samp_size(const TSampSizeArray& SS);
      int build_tree();      
      TMigrationMatrix& migr_matrix() {return (*MigrMatArray)[_cur_mig_mat];}
      MatElemType& min_coal_time(const int i, const int j) const {
      	return (*MinCoalTimes)(i,j);
      }
      MatElemType& mean_pair_diff(const int i, const int j) const {
      	return (*PairDiff)(i,j);
      }
      void print_gene_tree(ostream& os, const tree_output_type& tree_type, const float & mu) {
      	GeneTree.print_nodes(os, tree_type, mu);
      }
      //Recombination
      void print_gene_tree(ostream& os, const tree_output_type& tree_type, const float & mu, const int curLocus) {
         GeneTree.print_nodes(os, tree_type, mu, curLocus);
      }
      
      void print_sequences(ostream& os, const Mut_Type& mut_type) {
      	GeneTree.print_sequences(os, mut_type);
      };
      int implement_event(const int& cur_event);
	
		//MJJ 1/29/09
	int implement_snp_event(const int& cur_event);
	  
	  
	  //MJJ 11/13/08
	  int implement_plant_snp_event(const int& cur_event);
	  //MJJ End 11/13/08
	  
      void reset(const TDemeSizeArray& SA, const TGrowthArray& GA);
      
      
      int compute_moments_of_coalescence_times() {
      	return GeneTree.compute_moments_of_coalescence_times();
      };
      
      //Recombination: new parameter is a list of different statistic computed 
      int compute_moments_of_coalescence_times(double** list_stat_time) {
      	return GeneTree.compute_moments_of_coalescence_times(list_stat_time);
      };
       int compute_moments_of_demes_coalescence_times();
      double mean_coal_time() {return GeneTree.mean_coalescence_time();}
      double sd_coal_time() {return GeneTree.sd_coalescence_time();}

      const TNode& tree_node(const int& pos) {return GeneTree[pos];}
      TDeme& operator[](const int& d) {return Demes[d];}
      int print_coal_time_matrix(ostream& os ) {
      	if (!CoalTimes) return 0;
         os << "\n" << *CoalTimes;
         return 1;
      }
      int print_mean_pair_diff_matrix(ostream& os ) {
      	if (!PairDiff) return 0;
         os << "\n" << *PairDiff;
         return 1;
      }
      const TMigrationMatrix& coal_time_mat() const {return *CoalTimes;}
      const TMigrationMatrix& mean_pair_diff_mat() const {return *PairDiff;}  
      bool is_migration() {return _is_migration;}
      int check_migrations();
      int check_growth();
      int num_demes() const {return _num_demes;}
      const long& sprinkle_mutations(  const double& mut_rate, const int& numLinkedLoci,
                                       const TDataTypeArray* mut_type,
                                       const double& gamma_par, const double& mut_ratio,
                                       const double& prop_sites, const double& trans_rate,
                                       const int& range_constr, ofstream &snpfs); //snpfs added by MJJ 2/1/10




                     
		void write_samples_to_Arlequin_file(ostream& os, const Mut_Type& data_type);
      void write_group_section_to_Arlequin_file(ostream& os);
      void write_samples_to_PAUP_file(ostream& os, const Mut_Type& data_type);
      int migr_mat() {return _cur_mig_mat;}
      
      //Loro_19_06_01
		void write_locus_data_to_array(const Mut_Type& data_type);
      //Recombination: New function to deal with mix of different data types
      void write_locus_data_to_array(const TDataTypeArray* data_type_perLoc);

      
      void flushLoci();
      void write_loci_to_Arlequin_file(ostream& os, const Mut_Type& data_type, const int& genot_data);
      //recombination
      void write_loci_to_Arlequin_file(ostream& os, const Mut_Type& data_type, const int& genot_data, 
                                                                               const int& num_indep_loci);
      
      
      /*************************************************************************
                            New fubctions for recombination
      **************************************************************************/
      // Create the vector of the numbers of recombination hits between partially linked loci
      void allocate_recombHits();
      void allocate_recombHits(const int num_tot_interval);

      void initialize_recombHits();
    
      void add_recombHits(const TDeme* curDeme, const int recomb_pos);
      int return_recombHits(const int linked_loci);
      int number_recombHits_perDeme(const TDeme* curDeme);
      int number_recombHits_perDeme_perIndepLoci(const TDeme* curDeme, const int & indep_loci);

      void allocate_effective_recombHits();
      void add_no_effective_recombHits(const TDeme* curDeme, const int recomb_pos);
      int return_no_effective_recombHits(const int linked_loci);
      int number_no_effective_recombHits_perDeme(const TDeme* curDeme);
      int number_no_effective_recombHits_perDeme_perIndepLoci(const TDeme* curDeme, const int & indep_loci);
                                                                             
      void set_num_linked_loci(const int& numLinkedLoci);
      int get_num_linked_loci() {return num_linked_loci;};

      void set_num_linkage_blocks(int numLinkageBlocks) {num_linkage_blocks=numLinkageBlocks;};
      int get_num_linkage_blocks() {return num_linkage_blocks;};

      void set_num_indep_loci(int numIndepLoci) {num_indep_loci=numIndepLoci;};
      int get_num_indep_loci() {return num_indep_loci;};

      void set_rec_rate(double recRate) {rec_rate=recRate;};
      double get_rec_rate() {return rec_rate;};

      void  set_rec_rate_perLoc(const TRecRateArray* recRatePerLoc) {_rec_rate_perLoc=recRatePerLoc;};
      const TRecRateArray* get_rec_rate_perLoc() {return _rec_rate_perLoc;};

      void set_rangeConstraintPerLoc(const TNumLociArray* rangeConstraintPerLoc) {
         _rangeConstraintPerLoc=rangeConstraintPerLoc;
      }
      const TNumLociArray* get_rangeConstraintPerLoc() {return _rangeConstraintPerLoc;}

      void setTransitionRatePerLoc(const TMutRateArray* transRatePerLoc) {_transRatePerLoc=transRatePerLoc;};

      void  set_mut_rate_perLoc(const TMutRateArray* mutRatePerLoc) {_mut_rate_perLoc=mutRatePerLoc;};
      const TMutRateArray* get_mut_rate_perLoc() {return _mut_rate_perLoc;};

      void  set_geom_param_perLoc(const TMutRateArray* geomParamPerLoc) {_geom_param_perLoc=geomParamPerLoc;};
      const TMutRateArray* get_geom_param_perLoc() {return _geom_param_perLoc;};
      
      void set_num_linked_loci_perBlock(const TNumLociArray* numLinkedLociPerBlock) {_num_linked_loci_perBlock=numLinkedLociPerBlock;};
      const TNumLociArray* get_num_linked_loci_perBlock() {return _num_linked_loci_perBlock;};
      
      void set_freq_SNP_min(const double* FreqSNPmin) {freq_SNP_min=FreqSNPmin;};
      const double* get_freq_SNP_min() {return freq_SNP_min;};
      
      void set_sum_sample_size(int sumSampleSize) {sum_sample_size=sumSampleSize;};
      const int get_sum_sample_size() {return sum_sample_size;};

      int get_lociBlock(int loci) {
         int downBoundary=0, upBoundary=0;
         for(int i=0;i<get_num_linkage_blocks();++i) {
            upBoundary+=(*_num_linked_loci_perBlock->elem(i)-1);
            if(loci >= downBoundary && loci <= upBoundary) {
               return i;
            }
            ++upBoundary;
            downBoundary=upBoundary; 
         }
         return 0;
      }

      const TDeme& get_deme(const int & i);
      
      //For recombination: this variables were changed from private to public
      TTree	 GeneTree; 
      int    ind_loci;

      void resetRecCounts();
      
};

extern ostream& operator<<(ostream& os,const TDemeCollection& DC);

//------------------------------------------------------------------------------


#endif


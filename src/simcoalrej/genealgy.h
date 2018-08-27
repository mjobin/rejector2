#ifndef genealogy_h              // Sentry, use file only if it's not already included.
#define genealogy_h

#include <stdlib.h>
#include "arrays.h"
#include "mutation.h"
#include "chain.h"

using namespace std;

#define SNP_TOL 500 //MJJ defines tolerance in which a SNP will be allowed to drop

//Recombination
//enum tree_output_type {GENERATIONS=0, MUT_RATE=1, NUM_MUT=2};

extern float UNIT_TIME;

typedef MY_TArrayAsVector<int> TIntVect;
//Declaration of two classes defined elsewhere but needed here

typedef MY_TArrayAsVector<int>               TNumLociArray; //Loro_24_03_03 : Array of number of loci per blocks
typedef MY_TArrayAsVector<double>            TMutRateArray; //Loro_24_03_03 : Array of mutation rate per locus
typedef MY_TArrayAsVector<int>               TDataTypeArray; //Loro_24_03_03 : Array of data type per locus

class TDemeCollection;
class TMigrationMatrix;
typedef MY_TArrayAsVector<double>            TRecRateArray; //Loro_24_03_03 : Array of recombination rate per locus
//------------------------------------------------------------------------------
class TNode {
	public:   
   	long time;
      TNode *desc1;
      TNode *desc2;
      TNode *ancestor;
      TIntVect* desc1_nodes;
      TIntVect* desc2_nodes;
      TIntVect* sequence;     //The sequence to mutate
      static TIntVect hits;   //The number of hits per site
      static GammaRates mut_rates; //A static vector of gamma distributed rates
      int seq_length;
      static double mut_rate;			//Global mutation rate along the tree
      
      long num_new_mut;       //Number of mutations as compared to the direct ancestor
      long left_mut;
      long right_mut;
      int node_number;
      int deme;         		//Remember to which deme it belongs
      static double tree_length; //Length of the tree to which the node belongs, expressed in number of mutations
      static int min_mic;
      static int max_mic; //Range constraints for microsatellite data //Loro_26_2_99   
      static int num_linked_loci;         // number of partially linked loci
	  


   	TNode() {
			time			=  0;
         desc1		   =  NULL;
         desc2       =  NULL;
         ancestor    =  NULL;
         /**********************************************************************
         
             For recombination: don t forget to initialize the asc1 and asc2

         **********************************************************************/
         asc1        =  NULL;
         asc2        =  NULL;
         flag        =  NULL;
         
         desc1_nodes	=  NULL;
         desc2_nodes =  NULL;
         num_desc1	=  0;
         num_desc2   =  0;
         node_count++;
         node_number	=	node_count;
         deme			=	-1;
         num_new_mut	=  0;
         left_mut    =  0;
         right_mut   =  0;
         seq_length	=	0;
         sequence		=	NULL;
         _lidum		=  1L;
      /*************************************************************************   
      Initializing the new properties needed for recombination 
      *************************************************************************/
         min_pos     =  0 ;                              // min_pos and max_pos are the rank
         max_pos     =  num_linked_loci - 1;             // of the loci observed in the sample
                                                         // For the first generation 
                                                         // lye between 0 and the number of linked loci - 1

         flag        =  new bool[num_linked_loci];        // create an array of the state (observed or nonobserved
                                                         // in the sample) of every mutation.
         
         for(int cpt=0; cpt<num_linked_loci; cpt++){     // for every TNode created intialize this array to true
            flag[cpt]=true;
         }
         /*
         cout
            << "Node #" << node_number << ": " << num_linked_loci << " partially linked loci\n";
         for(int cpt=0; cpt<num_linked_loci; cpt++){     
            cout << "\t" << (cpt+1) << flag[cpt] << endl;
            
         }
         */
		 

		 
      }      
      TNode(const TNode& N)   {
      	if (this!=&N) *this=N;
      }
      ~TNode() {
      	if (desc1_nodes) {
            desc1_nodes->Flush();
         	delete desc1_nodes;
            desc1_nodes=NULL;
         }
   		if (desc2_nodes) {
            desc2_nodes->Flush();
         	delete desc2_nodes;
            desc2_nodes=NULL;
         }
   		if (sequence) {
            sequence->Flush();
         	delete sequence;
            sequence=NULL;
         }
         
        /***********************************************************************
        
             For recombination : Remember to delete the array of flag 

        ************************************************************************/
        if(flag){
          delete[] flag;
          flag=NULL;
        }
      }
      int operator==(const TNode& N) const {
      	if ( 	desc1==N.desc1  	&&
         		desc2==N.desc2		&&
               ancestor==N.ancestor)  return 1;
         return 0;
      }
      TNode& operator=(const TNode& N);
      
      void set_deme(const int& d) {deme=d;}
	  


      void set_mut_rate(const double& m) {mut_rate=m;}
      
      double tree_exp_length() {return tree_length;}
      void reset_tree_length() {tree_length=0;}
      int count_desc();
      //For recombination: called for every loci                            
      int count_desc(const int& loc);
      
      int count_polym_sites();
      
      int build_lists_of_descendent_nodes();
      int compute_moments_of_pairwise_divergence(double& sum_time, double& sum_square_time);
      //For recombination: called for every loci 
      int compute_moments_of_pairwise_divergence(double& sum_time, double& sum_square_time
                                                                 , const int& loc);

                                                
   	int compute_total_coal_times_among_demes( TDemeCollection& DemeCollec,
                                                TMigrationMatrix& CoalTimes,
                                                TMigrationMatrix& PairDiff, 
                                                TMigrationMatrix& MinCoalTimes);
      /*void print_desc_nodes(TDrawingBoard& DB, int node_posx, int node_posy,
      								char node, char hor_bar, char left_corner,
                              char right_corner, char vert_bar);
      */
      int print_info(ostream& os);
   	int num_desc1, num_desc2;
      void reset_node_count() {node_count=0;};
      long add_mutations(long *lidum, const int& num_mut,
      		const int& len, const Mut_Type& mut_type, const double& gamma_par,
            const double& trans_rate, const int& range_const);
      
      //For recombination: **********************************************************************************
      //for every loci (len=1) propagate nutation according to the coalescent tree

      //never called
      long add_mutations(long * lidum, const int& curLocus, const int& loc, const int& len, 
            const int& numLinkedLoci, const Mut_Type& mut_type, const double& mutRatePerLoc,
            const double& geomParamPerLoc, const double& gamma_par, const double& trans_rate,
            const int& range_const, const TNode* const pN);
      
      //used because faster than the previous
      long add_mutations(long *lidum, const int& num_mut, const int& loc/*, const int& len*/,
      	    const int& numLinkedLoci, const Mut_Type& mut_type, const double& mutRatePerLoc,
            const double& geomParamPerLoc, const double& gamma_par, const double& trans_rate,
            const int& range_const, bool&  dejaMut, const TNode* const pN, const int& SNPtime, const int& SNPdeme, ofstream &snpfs, const int& SNPeventnum); //MJJ dejaMut and SNPtime snp deme added 11/17/08

      //Guillaume 02 04 2004
      long add_mutations_without_recombination( long *                  lidum,
                                                int*                    num_mut,
                                                int&                    curLocus,
                                                //const int&            len,
                                                const int&              numLinkedLoci,
                                                const TDataTypeArray*   data_type,
                                                const TMutRateArray*    mutRatePerLoc,
                                                const TMutRateArray*    geomParamPerLoc,
                                                const double&           gamma_par,
                                                const TMutRateArray*    trans_rate,
                                                const TNumLociArray*    range_const,
                                                      TNode**            rootSNPs,
                                                 bool*                  dejamut,
                                                const TNode*            const pN,
												const int& SNPtime, const int& SNPdeme, ofstream &snpfs, const int& SNPeventnum  //MJJ 11/17/08
												);
      //******************************************************************************************************

      //A recursive routine to print tree topology and branch length
      void print_tree_structure(ostream& os, const tree_output_type& tree_type, const float & mu);
      //For recombination: print tree for every loci
      void print_tree_structure(ostream& os, const tree_output_type& tree_type, const float & mu,
                                                                                const int curLocus, 
                                                                                const TNode* root,
                                                                                const TNode* ancNode);

      /************************************************************************* 
            New properties and procedures necessary for recombinations
      *************************************************************************/
      
      int ID_Node;                        // variable for verifications
      int event;                          // variable for verifications
                                          // 0 sample, 1 coalescence, 2 recombination, 3 migration

      long _lidum;                        // dummy variable for the random number generator

      bool *flag;                         // loci observed in the sample: flag = TRUE
                                          // Otherwise                  : flag = FALSE

      int min_pos;                        // minimum position of the mutation observed (flag = TRUE)
      int max_pos;                        // maximum position of the mutation observed (flag = FALSE)

      double rec_proba;                   // Recombination probability for the whole sequence
                                          // rec_proba = (max_pos - min_pos) * recRate

      TNode *asc1;                        // Node have now 2 possible ascendants (recombination creates 2 new ancestors)
      TNode *asc2;                        

      //For every loci, randomly put one mutation somewhere on the tree
      //Used for SNP only
      long add_SNP_mutation(    long *lidum,
                                const int& num_mut,
                                const int& loc
                                /*, const int& len*/,
      		                const int& numLinkedLoci,
                                /*const Mut_Type& mut_type,*/
                                const float& gamma_par,
                                const float& trans_rate,
                                const int& range_const,
                                const TNode* const pN,
                                const double& SNP_mut_pos,
                                double& readLength,
                                bool& dejaMut,
                                TNode*& SNProot);
       
      int count_SNP_desc(const int& loc);
      int compute_TreeSize(double& sum_time, const int& loc);
      //Guillaume _05_04_2004
      int compute_TreeSize_SizeToNode_NumDescList(      int*   i,
                                                        int&   sum_time,
                                                        int*   sizeToNode,
                                                        //int*   numDescList,
                                                        TNode** scannedNodes);
      // Uniform distribution of recombination hits: Not used      
      int recombinationPosition(const int & min,const int & max );
      // Non Uniform distribution: chromosomes blocks with different recombination rates 
      int recombinationPosition(const int& min, const int& max,                    
                                const TRecRateArray* RecRatePerLoc,
                                double recProba);
                                
      double update_recombination_proba(const int& min, 
                                  const int& max,
                                  const TRecRateArray* RecRatePerLoc);
      
      //to initialize the static member num_linked_loci, once for all
      void initialise_first_Node(const int& nLoc) {
        num_linked_loci  =   nLoc;
      }

      //to destroy TNode: Not used
      void destroy();
      
      
   private:
   	static int node_count;
};

ostream & operator<<(ostream& os, const TNode& node);   

//------------------------------------------------------------------------------
//An array of pointers to Nodes
//It is like an indirect array of Nodes where Nodes would not be deleted
//in the destructor of the indirect array

class TListNode {
	public:
		TListNode() {list=NULL; ListSize=0;}
		TListNode(const int& size, TNode * tree);
		~TListNode() {
			if (list) delete[] list;
      }
      TNode*& operator[](const int& n) const {return list[n];}
      int ListSize;
      int allocate_list(const int& size);
      void de_allocate_list() {
      	if (list) delete[] list;
         list=NULL;
      };
      const int& size() const {return _size;}
	private:
   	TNode **list;
      int _size;
};

typedef TISimpleChain<TNode>	tree_chain;

//------------------------------------------------------------------------------
class TTree {
   public:
      TTree(): chainTree(true) {
      SampleSize=0;
      /***************************************************************************

      For recombination : initiale pointer to the chain list

      ****************************************************************************/
      //tree=NULL;
      delete_nodes=1;
      _mean_coalescence_time=0;
      _sd_coalescence_time=0;

      //cout << "Enter TTRee constructor";
      // For Recombination: To stop simulation when every linked loci has coalesced
      num_linked_loci=0;
      DeathVect = NULL;    
      MRCA_list = NULL;
      };
      TTree(const int& size);

      ~TTree() {

      //cout << " enter TTree destructor" << endl;
       /***************************************************************************

      For recombination : delete the chain list chainTree

      ****************************************************************************/ 
         // Recombination: Don t forget to delete these vectors 
         if(DeathVect){
            delete[] DeathVect;
            //DeathVect=NULL;
         }
         if(MRCA_list){
            delete[] MRCA_list;
            //MRCA_list=NULL;  
         }

         //if(tree_length_perLoc){
         //   delete[] tree_length_perLoc;
         //   //MRCA_list=NULL;  
         //}
         //if(tot_mut_perLoc){
         //   delete[] tot_mut_perLoc;
         //   //MRCA_list=NULL;  
         //}

      }
      /***************************************************************************

      For recombination : the operator [] returns the TNode 
                           of the ith element of the chain list

      ****************************************************************************/
      //TNode& operator[](const int& n) {return tree[n];}
      TNode& operator[](const int& n) {
        return chainTree(n);
      }

      /***************************************************************************

      For recombination : the operator () returns a pointer to the TNode 
                           of the ith element of the chain list

      ****************************************************************************/
      //TNode* operator()(const int& i) {return tree + i;}
      TNode* operator()(const int& i) {
         TNode *pointer = (TNode*) chainTree[i];
         return pointer;
      }

      int allocate_tree(const int& sampsize);
      
      int bottleneck(double time, double factor);
      void print_nodes(ostream& os, const tree_output_type& tree_type, const float & mu);
      void print_sequences(ostream& os, const Mut_Type& mut_type);
      int sample_size() const {return SampleSize;}
      void set_delete_nodes() {delete_nodes=1;}
      void set_not_delete_nodes() {delete_nodes=0;}
      
      /***************************************************************************

      For recombination : Create a chained list of TNodes with size SampleSize

      ****************************************************************************/
      void print_nodes(ostream& os, const tree_output_type& tree_type, const float & mu, 
                                                                       const int curLocus); 

      //float total_length() {if (tree) return tree[0].tree_length; else return 0;};
      float total_length() {
         if ( chainTree.get_beg() ) return chainTree(0).tree_length; 
         else return 0;
      }

      int compute_moments_of_coalescence_times();
      //For recomnbination
      int compute_moments_of_coalescence_times(double** list_stat_time);

      const double& mean_coalescence_time() {return _mean_coalescence_time;}
      const double& sd_coalescence_time() {return _sd_coalescence_time;}

   private:
   	int SampleSize;

   /***************************************************************************

   For recombination : Dont use tree

   ****************************************************************************/
      //TNode* tree;

      
   public:
      tree_chain chainTree;
      int delete_nodes;
      double _mean_coalescence_time;
      double _sd_coalescence_time;
      double _tree_subst_length;
      double _num_mut_on_tree;
      
      int num_linked_loci;
      int *DeathVect;            // A vector storing the current numbers of TNodes
                                 // for ever linked loci (initialized by the number of
                                 // genes). These numbers decrease when a coalescence event occurs
                                 // The simulations will be stopped
                                 // when there are 1 TNode for every linked loci
                                 // This TNode is the MRCA for every linked locus

      TNode **MRCA_list;
      
      //double* tree_length_perLoc;    //kept for later
      //double* tot_mut_perLoc;
};

//------------------------------------------------------------------------------
#endif  // genealogy_h sentry.


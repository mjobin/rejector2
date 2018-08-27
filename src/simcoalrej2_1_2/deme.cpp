////////////////////////////////////////////////////////////////////////////////
//
// A class implementing a deme within a population, exchanging migrants with
// other demes
//
//
////////////////////////////////////////////////////////////////////////////////

#include "cond_var.h"

// This put in place on OS X in lieu of including values.h, which is not included.
#define	MAXLONG		LONG_MAX

#ifdef _SHOW_MUT_DISTR_
	#include <fstream> //MJJ 2/1/09
	ofstream f_site_hits("mut_hits_distr.sum");
	static count_simul=0;
#endif

#include "deme.h"
#include "genealgy.h"

#include <strstream>


//------------------------------------------------------------------------------

//Initialization of static members of TDeme
TMigrationMatrix* TDeme::Migrations=NULL;

//==============================================================================
//------------------------------------------------------------------------------
TDeme&
TDeme::operator=(const TDeme& D) {
	_growth_rate=D._growth_rate;
   _sample_size=D._sample_size;
   _deme_size=D._deme_size;
   _lineages=D._lineages;
   _lidum=D._lidum;
   _id=D._id;
   _min_deme_size=D._min_deme_size;
   _max_deme_size=D._max_deme_size;
   _proba_coal=D._proba_coal;
   Migrations=D.Migrations;
   _is_migration=D._is_migration;

   
   /**********************************************************************
         
            For recombination : don t use  NodeList and OriginalNodes
         
   **********************************************************************/
   //NodeList.de_allocate_list();
   ChainNodeList.Call_deleteChain();

   
   //int size=D.NodeList.ListSize;
   int size=D.ChainNodeList.get_size();
  
   if (size) {
   	//NodeList.allocate_list(size);
   	for (int i=0; i<size; ++i) {
      	//NodeList[i]=D.NodeList[i];
         //TNode *pointer= (TNode*) D.ChainNodeList[i]; 
         //ChainNodeList.Add( pointer );
      }
   }
   //NodeList.de_allocate_list();
   ChainNodeList.Call_deleteChain();


   
   //size=D.OriginalNodes.ListSize;
   size=D.ChainOriginalNodes.get_size();
   if (size) {
   	//OriginalNodes.allocate_list(size);
   	for (int i=0; i<size; ++i) {
      	//OriginalNodes[i]=D.OriginalNodes[i];
         //TNode *pointer= (TNode*) D.ChainOriginalNodes[i]; 
         //ChainOriginalNodes.Add( pointer ); 
      }
   }

   _polym_sites=D._polym_sites;
   return *this;
}
//------------------------------------------------------------------------------
//REM: Function should also work for TRecNode and TRecTree as is
/*
int
TDeme::copy_node_pointers(TTree* T, const int& curr_lin, const long& max_lin) {

   int base=curr_lin;

   //This node list is dynamic and will change according to migrations and coalescences
   NodeList.allocate_list(max_lin);
   

   //Keeps a copy of the original nodes
   OriginalNodes.allocate_list(_sample_size);

	for (int i=0; i<_sample_size; ++i, ++base) {
   	//Stamp the node as originally belonging to the current deme
   	(*T)(base)->set_deme(_id);
   
      //Then copy the node pointer to the local lists
   	NodeList[i]=(*T)(base);
   }

   //Filling up the node list once for all with empty pointers
   //so that no further size change is required during the generation
   //of the genealogy
   for (int i=_sample_size; i<max_lin; ++i) {
   	NodeList[i]=NULL;
   }
   
   return 1;
}
*/ 
//------------------------------------------------------------------------------
//for recombination
int
TDeme::copy_node_pointers(TNode* N, TTree * T) {

	for (int i=0; i<_sample_size; ) {
   	        //Stamp the node as originally belonging to the current deme
   	        N->set_deme(_id);
                //Then copy the node pointer to the local lists
                ChainNodeList.Add(N);
   	        ChainOriginalNodes.Add(N);
                if (++i<_sample_size) N=T->chainTree.next();
        }
  	return 1;
}
//------------------------------------------------------------------------------
//Check and implements migrations
int
TDeme::send_migrants(TDeme& sink, const int& i) {
   if (!_lineages) return 0;
   //We will allow for multiple migrations per generation
   //Loro: !!! A faster routine could be made by first getting the number of migrants
   //from a binomial distribution
   for (int j=0, num_lineages=_lineages; j<num_lineages; ++j) {
   	if (ran3(&_lidum)<=(*Migrations)(_id,i)){
         migrate(sink);
      }
   }
   return 1;
}

//------------------------------------------------------------------------------
//Implements the coalescence process within population
/**
int
TDeme::coalesce_lineages(const long& time, TTree& GeneTree, const int& pos) {

   //If there is only one lineage left, do not even try to coalesce it...
   if (_lineages<2) return 0;

	if (!_deme_size) {
      cout 	<< "TDeme::coalesce_lineages() : deme size of population "
      		<< _id << " is zero !" << endl;
   	return 0;
   }

   double proba_coal=_proba_coal*0.5*_lineages*(_lineages-1);

   if (proba_coal==0.0) return 0;

   //Coalesce two lineages with probability proba
   if (ran3(&_lidum)<=proba_coal) {
   	//Pick two lineages at random for a coalescence
      int pick1= (int) (ran3(&_lidum)*_lineages), pick2= (int) (ran3(&_lidum)*_lineages);
      long count=0;
      while (pick1==pick2) {
      	pick2=(int) (ran3(&_lidum)*_lineages);
         if (++count>1000) return 0;
      }

      TNode *desc1=NodeList[pick1],
      		*desc2=NodeList[pick2],
            *ancestor=GeneTree(pos);

      //Update node information
      ancestor->desc1=desc1;
      ancestor->desc2=desc2;
      ancestor->time=time;
      desc1->ancestor=ancestor;
      desc2->ancestor=ancestor;

		//Replace one coalescing node by the ancestor
      //and swap the second with the last node after updating the number
      //of remaining lineages
      NodeList[pick1]=ancestor;
      NodeList[pick2]=NodeList[--_lineages];
      return 1;
   }
   return 0;
}
*/

//------------------------------------------------------------------------------
//Implements the coalescence process within population
//Taking into account the recombination between loci
int
TDeme::singleCoalescentEvent(const long& time, TTree& GeneTree, TDemeCollection* myDemes ) {

   //If there is only one lineage left, do not even try to coalesce it...
   if (_lineages<2) return 0;

	if (!_deme_size) {
      cout 	<< "TDeme::coalesce_lineages() : deme size of population "
      		<< _id << " is zero !" << endl;
   	return 0;
   }

   double proba_coal=_proba_coal*0.5*_lineages*(_lineages-1);

   if (proba_coal==0.0) return 0;

   //Coalesce two lineages with probability proba
   
   if (ran3(&_lidum)<=proba_coal) {
      //Pick two lineages at random for a coalescence
      int pick1= (int) (ran3(&_lidum)*_lineages);
      int pick2= (int) (ran3(&_lidum)*_lineages);      

      long count=0;
   
      while (pick1==pick2 || pick1>=_lineages || pick2>=_lineages) {
        pick1=(int) (ran3(&_lidum)*_lineages);
        pick2=(int) (ran3(&_lidum)*_lineages);
        count++;
        if (++count>1000) return 0;
      }
   
      /****************************************************************************

      For recombination : use a chained list ChainNodeList rather than the NodeList 
   
      ****************************************************************************/

      TNode  *desc1, *desc2, *ancestor;
      
      ancestor = new TNode;                                                     // For recombination add a new
      GeneTree.chainTree.Add(ancestor);                                         // TNode in the Chained List of TTree
      ancestor->ID_Node=(GeneTree.chainTree.get_size() - 1 );
      ancestor->event = 1;

      ChainNodeList.returnTwoElemsAndReplaceFirstByNewElem(desc1, desc2, ancestor, pick1, pick2);
      --_lineages;

      //Update node information
      ancestor->desc1=desc1;
      ancestor->desc2=desc2;
      ancestor->time=time;

//MJJ
		ancestor->deme=_id;





         
      desc1->ancestor=ancestor;
      desc1->asc1=ancestor;
      desc1->asc2=NULL;

      desc2->ancestor=ancestor;
      desc2->asc1=ancestor;
      desc2->asc2=NULL;

      //Update the flag of the TNode after a coalescence event
      for (int cpt=0; cpt<ancestor->num_linked_loci; cpt++){
         if( desc1->flag[cpt] || desc2->flag[cpt] ){
            ancestor->flag[cpt]=true;
         }
         else{
            ancestor->flag[cpt]=false;
         }
         if( desc1->flag[cpt] && desc2->flag[cpt] ){
            myDemes->GeneTree.DeathVect[cpt]--;
            if(myDemes->GeneTree.DeathVect[cpt]==1){
               myDemes->GeneTree.MRCA_list[cpt]=ancestor;
            }
         }
      }
                                                               
      if( desc1->min_pos < desc2->min_pos) {
         ancestor->min_pos=desc1->min_pos;
      }
      else {
         ancestor->min_pos=desc2->min_pos;
      }
      if( desc1->max_pos > desc2->max_pos) {
         ancestor->max_pos=desc1->max_pos;
      }
      else {
         ancestor->max_pos=desc2->max_pos;
      }

      return 1;
   }
   return 0;
}
//------------------------------------------------------------------------------
//Implements the coalescence process within population
//Taking into account the recombination between loci
int
TDeme::multipleCoalescentEvents(const long& time, TTree& GeneTree, TDemeCollection* myDemes ) {

   //If there is only one lineage left, do not even try to coalesce it...
   if (_lineages<2) return 0;

	if (!_deme_size) {
      cout 	<< "TDeme::coalesce_lineages() : deme size of population "
      		<< _id << " is zero !" << endl;
   	return 0;
   }
   DescNodeVect curDescNodeVect(_lineages-1, 0);
   
   //For every Node in the current generation
   //draw 1 ancestor among the _demesize ancestors in the previous generation 
   ChainNodeList.resetIterator();
   TNode * currentNode= (TNode*) ChainNodeList[0];
   for(int cpt=0; cpt < _lineages; ++cpt, currentNode=ChainNodeList.next()) {
      int tempanc= (int) (ran3(&_lidum)*_deme_size);
      int anc= (int) tempanc;
      //int anc= ran3(&_lidum)*(_deme_size);
      curDescNodeVect.Add(TDescNode(anc,currentNode));
   }
   ChainNodeList.resetIterator();
   curDescNodeVect.quickSortArray();
   
   //Find couple of Nodes sharing the same ancestor
   bool  oneCoal=false;
   int curNumLineages=_lineages;
   for (int i=1; i<curNumLineages; ++i) {
      if(curDescNodeVect[i] == curDescNodeVect[i-1] ) {
         oneCoal=true;

         TNode  *desc1, *desc2, *ancestor;
         desc1=curDescNodeVect[i].pDescNode;
         desc2=curDescNodeVect[i-1].pDescNode;

         ancestor = new TNode;                                             // For recombination add a new
         GeneTree.chainTree.Add(ancestor);                                 // TNode in the Chained List of TTree
         ancestor->ID_Node=(GeneTree.chainTree.get_size() - 1 );
         ancestor->event = 1;

         ChainNodeList.giveTwoElemsAndReplaceFirstByNewElem(desc1, desc2, ancestor);
         --_lineages;

         curDescNodeVect[i].pDescNode=ancestor;
         
         //Update node information
         ancestor->desc1=desc1;
         ancestor->desc2=desc2;

         ancestor->time=time;
		 
		 //MJJ
		ancestor->deme=_id;

		 
         desc1->ancestor=ancestor;
         desc2->ancestor=ancestor;

         desc1->asc1=ancestor;
         desc1->asc2=NULL;
         desc2->asc1=ancestor;
         desc2->asc2=NULL;
         
         //Update the flag of the TNode after a coalescence event
         for(int cpt=0; cpt<ancestor->num_linked_loci; ++cpt){
            if( desc1->flag[cpt] || desc2->flag[cpt] ){
               ancestor->flag[cpt]=true;
            }
            else {
               ancestor->flag[cpt]=false;
            }
            if( desc1->flag[cpt] && desc2->flag[cpt] ) {
               myDemes->GeneTree.DeathVect[cpt]--;
               if(myDemes->GeneTree.DeathVect[cpt]==1) {
                  myDemes->GeneTree.MRCA_list[cpt]=ancestor;
               }
            }
         }            
         
         if( desc1->min_pos < desc2->min_pos) {
            ancestor->min_pos=desc1->min_pos;
         }
         else {
            ancestor->min_pos=desc2->min_pos;
         }

         if( desc1->max_pos > desc2->max_pos) {
            ancestor->max_pos=desc1->max_pos;
         }
         else {
            ancestor->max_pos=desc2->max_pos;
         }
      }
   }

   if (oneCoal) {
      return 1;
   }
   else {
      return 0;
   }
   
}

//------------------------------------------------------------------------------
void
TDeme::reset() {

   /****************************************************************************
                              Old version without recombination
    
   //Restoring initial number of lineages
	_lineages=_sample_size;

   int size=NodeList.size();
   //Copying original nodes into dynamic list
	for (int i=0; i<_sample_size; ++i) {
   	NodeList[i]=OriginalNodes[i];
   }
   //Filling up the node list once for all with empty pointers
   for (int i=_sample_size; i<size; ++i) {
   	NodeList[i]=NULL;
   }
   _min_deme_size=MIN_SIZE;
   _max_deme_size=MAX_SIZE;
   _proba_coal=1.0/_deme_size;
   ****************************************************************************/


  /****************************************************************************
                          New version for recombination
        Use the chained list ChainNodeList and ChainOriginalNodes
   ****************************************************************************/
   //Restoring initial number of lineages
   _lineages=_sample_size;

   

   ChainNodeList.Call_deleteChain();
   //Wilson
   if(_sample_size){
        ChainOriginalNodes.resetIterator();
        TNode * curNode=ChainOriginalNodes.current();
         for (int i=0; i<ChainOriginalNodes.get_size(); ++i,curNode=ChainOriginalNodes.next() ) {
        ChainNodeList.Add(curNode);
        }
        ChainNodeList.resetIterator();
        ChainOriginalNodes.resetIterator();
   }
   
   _min_deme_size=MIN_SIZE;
   _max_deme_size=MAX_SIZE;
   _proba_coal=1.0/_deme_size;

}
//------------------------------------------------------------------------------
//A simple procedure to check whether there are migrations from this population
int
TDeme::check_migrations(){
   if (Migrations) {
      _is_migration=false;
   	int size=Migrations->size();
      for (int i=0; i<size; ++i)
      if (_id!=i) {
      	if ((*Migrations)(_id,i)>0.0) {
         	_is_migration=true;
            return 1;
         }
      }
      return 0;
   }
   return 0;
}
//------------------------------------------------------------------------------
void
TDeme::print_haplotypes_for_Arlequin(ostream & os, const int& n, const Mut_Type& data_type) {

   char DNA_letters[4]={'A','G','C','T'};

   /****************************************************************************
   
                       For recombination Use the chained list
                       ChainNodeList and ChainOriginalNodes
                       
   ****************************************************************************/
   //TNode* anyNode=OriginalNodes[0];
   TNode* anyNode= (TNode*) ChainOriginalNodes[0];

   int len=anyNode->seq_length;

   os << "#Number of mutation hits per site\n"
   	<< "#Sites ";
   for (int j=0; j<len; ++j) {
      	os
            << setw(4)
            << setiosflags(ios::right)
         	<< j << " ";
      }      
   os << "\n";

   #ifdef _SHOW_MUT_DISTR_
   	f_site_hits << "Simulation #" << ++count_simul << "\t";
   #endif

   os << "#Hits  ";
   for (int j=0; j<len; ++j) {
      	os 
            << setw(4)
            << setiosflags(ios::right)
            << anyNode->hits[j] << " ";

         #ifdef _SHOW_MUT_DISTR_
      	f_site_hits
            << setiosflags(ios::right)
            << anyNode->hits[j] << "\t";
         #endif
      }
   os << "\n";

   #ifdef _SHOW_MUT_DISTR_
   f_site_hits  << "\n";
   #endif

   for (int i=0; i<_sample_size; ++i) {
   	//Write id and frequency
   	os << n << "_" << i << "\t1 ";
      
      //TNode* CurNode=OriginalNodes[i];
      TNode* CurNode=(TNode*) ChainOriginalNodes[i];

	   os << "Called here?" << endl;
      
      if (CurNode->sequence) {
      	int* cursite=CurNode->sequence->begin();
      	for (int j=0; j<len; ++j) {
         	if (data_type==DNA || data_type==SNP)
   					//os << DNA_letters[(*CurNode->sequence)[j]];
   					os << DNA_letters[*cursite++];

					//os << DNA_letters[*cursite++] << " ";  //MJJ 6/3/09
               else
   					//os << (*CurNode->sequence)[j];
                  //Loro_1_3_99 : In order to avoid negative numbers
               	if (data_type==MICROSAT) os << (500+*cursite++) << " ";
               	else os << *cursite++;
         }
         os << "\n";
      }
   }
}
//------------------------------------------------------------------------------
int findPolymorphicLoci(const	TLocusArray& Loci, const int& sampSize,
									const int& numLoci, bool* polymFlag     ) {
   if (!numLoci) return 0;
   int numPolymLoci=0;
   int numLinkedLoci=Loci[0].data[0].length();
	for (int i=0, pos=0; i<numLoci; ++i) {
   	for (int j=0; j<numLinkedLoci; ++j, ++pos) {
      	//Look if there is more than one allele in the sample
         bool isPolym=false;
         char ref=Loci[i].data[0][j];
         for (int k=1; k<sampSize && !isPolym; ++k){
         	if (Loci[i].data[k][j]!=ref) {
            	isPolym=true;
               ++numPolymLoci;
            }
         }
         polymFlag[pos]=isPolym;
      }
   }
   return numPolymLoci;
}
//------------------------------------------------------------------------------
//Transforme les donnees ACGT en donnees binaires 0/1.
int findPolymorphicLociAndChangeToZeroOne( TLocusArray& Loci, const int& sampSize,
									            const int& numLoci, bool* polymFlag) {
   if (!numLoci) return 0;
   int numPolymLoci=0;
   int numLinkedLoci=Loci[0].data[0].length();
	for (int i=0, pos=0; i<numLoci; ++i) {
   	for (int j=0; j<numLinkedLoci; ++j, ++pos) {
      	//Look if there is more than one allele in the sample
         bool isPolym=false;
         char ref=Loci[i].data[0][j];
         //Set reference character to Zero
         Loci[i].data[0][j]='0';
         for (int k=1; k<sampSize; ++k){
         	if (Loci[i].data[k][j]!=ref) {
            	isPolym=true;
               Loci[i].data[k][j]='1';
            }
            else Loci[i].data[k][j]='0';
         }
         if (isPolym) ++numPolymLoci;
         polymFlag[pos]=isPolym;
      }
   }
   return numPolymLoci;
}
//------------------------------------------------------------------------------
//Transforme les donnees ACGT en donnees binaires 0/1.
//Recombination: to deal with mix of different loci
int findPolymorphicLociAndChangeToZeroOne( TLocusArray& Loci, const int& sampSize,
									               const int& indLoci, 
                                          bool* polymFlag, 
                                          int num_linked_loci, 
                                          int& gap) {
   //if (!numLoci) return 0;
   int i=indLoci;
   int numPolymLoci=0;
   int numLinkedLoci=Loci[0].data[0].length();
   //int numLinkedLoci=num_linked_loci;
   if(gap) {
      ++gap;
   }
	//for (int i=0, pos=gap; i<numLoci; ++i) {
      bool Is_there_blank=false;
   	for (int j=gap, pos=gap; j<numLinkedLoci; ++j, ++pos) {
      	//Look if there is more than one allele in the sample
         bool isPolym=false;
         char ref=Loci[i].data[0][j];
         char blank=' ';
         if(ref==blank) {
            Is_there_blank=true;
            gap=j;
            polymFlag[pos]=true;
            return numPolymLoci;
         }
         //Set reference character to Zero
         Loci[i].data[0][j]='0';
         for (int k=1; k<sampSize; ++k){
         	if (Loci[i].data[k][j]!=ref) {
            	isPolym=true;
               Loci[i].data[k][j]='1';
            }
            else Loci[i].data[k][j]='0';
         }
         if (isPolym) ++numPolymLoci;
         polymFlag[pos]=isPolym;
         if(j==(numLinkedLoci-1) && !Is_there_blank) {
            gap=numLinkedLoci;
         }
      }
   //}
   return numPolymLoci;
}
//Recombination: to deal with mix of different loci
//------------------------------------------------------------------------------
int findPolymorphicLoci(const	TLocusArray& Loci, const int& sampSize,
                        const int& indLoci, 
                        bool* polymFlag, 
                        int num_linked_loci, 
                        int& gap) {
   //if (!numLoci) return 0;
   int i=indLoci;
   int numPolymLoci=0;
   int numLinkedLoci=Loci[0].data[0].length();
   //int numLinkedLoci=num_linked_loci;
   if(gap) {
      ++gap;
   }
	//for (int i=0, pos=gap; i<numLoci; ++i) {
   	for (int j=gap, pos=gap; j<numLinkedLoci; ++j, ++pos) {
      	//Look if there is more than one allele in the sample
         
         //bool isPolym=false;
         
         char ref=Loci[i].data[0][j];
         char blank=' ';
         if(ref==blank) {
            gap=j;
            polymFlag[pos]=true;
            return numPolymLoci;
         }
         //for (int k=1; k<sampSize && !isPolym; ++k){
         //	if (Loci[i].data[k][j]!=ref) {
         //   	isPolym=true;
         //      ++numPolymLoci;
         //   }
         //}
         polymFlag[pos]=true;    //flag  put only at the firts character of the Microsat name
         numPolymLoci=1;            //Recombination 25 04 03:
                                    //This is wrong !!!!
                                    //for instance we assume that every microsat loci is polymorphe
      }                             //in the sample
   //}                          
   return numPolymLoci;
}
//------------------------------------------------------------------------------
//modified for recombination
void
TDeme::print_loci_for_Arlequin(ostream & os, const int& deme, const Mut_Type& data_type,
											            const int& genot_data) {
   int num_loci=_loci.GetItemsInContainer();
   int numLinkedLoci=_loci[0].data[0].length();
   int num_polym_loci;

   bool* isPolymorphic=NULL;     //Array of polymorphism flags

   #ifdef _PRINT_SNP_ONLY_
   if (data_type==DNA || data_type==SNP) {
   	isPolymorphic= new bool[num_loci*numLinkedLoci];
		num_polym_loci=findPolymorphicLociAndChangeToZeroOne(_loci, _sample_size,
      																		num_loci, isPolymorphic);
      os << "#Number of polymorphic sites : " << num_polym_loci << "\n";
   }
   #endif

   for (int i=0, j=0; i<_sample_size; ++i) {
   	//Write id and frequency
   	if (genot_data) {
        if (!(i % 2) ) {
        	os << (deme+1) << "_" << (j+1) << "\t1\t";
         ++j;
        }
        else os << "\t\t";
     }
     else os << (deme+1) << "_" <<  (i+1) << "\t1\t";
     for (int m=0, pos=0; m<num_loci; ++m) {
         if (!isPolymorphic) os << _loci[m].data[i] << "\t";
         else {  //More elaborate printing to print SNP data only
         	for (int k=0; k<numLinkedLoci; ++k, ++pos) {
            	if (isPolymorphic[pos]) os << _loci[m].data[i][k];
            }
         }
     }
     os << "\n";
   }

   #ifdef _PRINT_SNP_ONLY_
   if ((data_type==DNA || data_type==SNP) && isPolymorphic) delete [] isPolymorphic;
   #endif

}
//------------------------------------------------------------------------------
//modified for recombination
void
TDeme::print_loci_for_Phase(ostream & os, const int& deme, const Mut_Type& data_type) {

   int num_loci=_loci.GetItemsInContainer();
   int numLinkedLoci=_loci[0].data[0].length();

   bool* isPolymorphic=new bool[num_loci*numLinkedLoci];     //Array of polymorphism flags

   int numPolymLoci;

   if (data_type==DNA || data_type==SNP) {
		numPolymLoci=findPolymorphicLociAndChangeToZeroOne(_loci, _sample_size,
      																	num_loci, isPolymorphic);
   }
   else numPolymLoci=findPolymorphicLoci(_loci, _sample_size, num_loci, isPolymorphic);

   //Print header
   os << (_sample_size/2) << "\n" << numPolymLoci << "\n";

	for (int i=0; i<numPolymLoci; ++i)  {

   	switch (data_type) {

      	case DNA       : os << "S"; break;

      	case SNP       : os << "S"; break;

         case MICROSAT  : os << "M"; break;

         default 			: os << "S"; break;

      }

   }

   os << "\n";

   for (int i=0, j=0; i<_sample_size; ++i) {
   	//Write id
      if (!(i % 2) ) {
      	os << (deme+1) << "_" << (j+1) << "\n";
         ++j;
      }
     	for (int m=0, pos=0; m<num_loci; ++m) {
         if (!isPolymorphic)
        	os << _loci[m].data[i] << "\t";
         else {  //More elaborate printing to print SNP data only
         	for (int k=0; k<numLinkedLoci; ++k, ++pos) {
            	if (isPolymorphic[pos]) os << _loci[m].data[i][k];
            }
         }
     }
     os << "\n";
   }

   if (isPolymorphic) delete [] isPolymorphic;
}
//------------------------------------------------------------------------------
//Overload for recombination: Mix of different loci
void
TDeme::print_loci_for_Phase(ostream & os, const int& deme, const int& num_linked_loci, 
                                                           const TDataTypeArray* data_type_perLoc) {

   int num_loci=_loci.GetItemsInContainer();
   int numLinkedLoci=_loci[0].data[0].length();
   //int numLinkedLoci=num_linked_loci;

   
   bool* isPolymorphic=new bool[num_loci*numLinkedLoci];     //Array of polymorphism flags

   int numPolymLoci;

   
   
   numPolymLoci=0;
   int first_loci_in_curBlock=0;
   for(int i=0; i<num_loci; ++i) {
      for(int j=0; j<num_linked_loci;) {
         if (*data_type_perLoc->elem(j)==2 || *data_type_perLoc->elem(j)==3) {
		      numPolymLoci+=findPolymorphicLociAndChangeToZeroOne(_loci, _sample_size,
      																	      i,
                                                               isPolymorphic, num_linked_loci,   
                                                               first_loci_in_curBlock);
            j=first_loci_in_curBlock;
         }
         else {
            numPolymLoci+=findPolymorphicLoci(_loci, _sample_size, 
                                                      i, 
                                                      isPolymorphic, num_linked_loci,   
                                                      first_loci_in_curBlock);
            ++j;
         }
      }
   }
   //Print header
   os << (_sample_size/2) << "\n" << numPolymLoci << "\n";

	for (int i=0; i<num_linked_loci; ++i)  {     // all markers are polymorphes

   	   switch (*data_type_perLoc->elem(i)) { 

            case 0         : os << "M"; break; 

      	   case 1         : os << "S"; break;

      	   case 2         : os << "S"; break;

      	   case 3         : os << "S"; break;

            default        : os << "S"; break;

         }

   }

   os << "\n";

   for (int i=0, j=0; i<_sample_size; ++i) {
   	//Write id
      if (!(i % 2) ) {
      	os << (deme+1) << "_" << (j+1) << "\n";
         ++j;
      }
     	for (int m=0, pos=0; m<num_loci; ++m) {
         if (!isPolymorphic)
        	os << _loci[m].data[i] << "\t";
         else {  //More elaborate printing to print SNP data only
         	for (int k=0; k<numLinkedLoci; ++k, ++pos) {
            	if (isPolymorphic[pos]) os << _loci[m].data[i][k];
            }
         }
     }
     os << "\n";
   }

   if (isPolymorphic) delete [] isPolymorphic;
}
//------------------------------------------------------------------------------
//Modified for recombination
#pragma argsused
void
TDeme::print_loci_for_Haplotyper(ostream & os, const int& n) {

   int num_loci=_loci.GetItemsInContainer();
   int numLinkedLoci=_loci[0].data[0].length();

   bool* isPolymorphic=new bool[num_loci*numLinkedLoci];     //Array of polymorphism flags

   findPolymorphicLociAndChangeToZeroOne(_loci, _sample_size,
   																		 num_loci, isPolymorphic);
   
   char ref, all[2];
	for (int i=0; i<_sample_size; ++i) {
     	for (int m=0, pos=0; m<num_loci; ++m) {
      	for (int k=0; k<numLinkedLoci; ++k, ++pos) {
      		if (isPolymorphic[pos]) {
               //Consider the allele of the first gamete as the reference
               //like for PHASE
         		if (i==0) ref=_loci[m].data[0][k];
               all[0]=_loci[m].data[i][k];
               all[1]=_loci[m].data[i+1][k];
            	if (all[0]==all[1]) {
               	if (all[0]==ref) os << '1';
                  else os << '2';
               }
               else os << '0';
            }
         }
      }
      ++i; //Also increment i here
   	os << "\n";
   }

   if (isPolymorphic) delete [] isPolymorphic;
}
//------------------------------------------------------------------------------ 
//Modified for recombination
void
TDeme::print_haplotypes_to_locus(const Mut_Type& data_type) {

   TLocus curLocus(_sample_size);
 
   char DNA_letters[4]={'A','G','C','T'};
   
   //int len=OriginalNodes[0]->seq_length;
   int len=ChainOriginalNodes[0]->seq_length;

   

   for (int i=0; i<_sample_size; ++i) {
      //TNode* CurNode=OriginalNodes[i];
      TNode* CurNode= (TNode*) ChainOriginalNodes[i];

      
      strstream curHaplStream;
      my_string curHapl;
	   
	   

      if (CurNode->sequence) {
      	int* cursite=CurNode->sequence->begin();
      	for (int j=0; j<len; ++j) {
         	if (data_type==DNA  || data_type==SNP)
   					curHaplStream << DNA_letters[*cursite++];
               else
               	if (data_type==MICROSAT) curHaplStream << (500+*cursite++) << " ";
               	else curHaplStream << *cursite++;
         }
      }
      curHapl.read_line(curHaplStream);
      curLocus.data.Add(curHapl);
   }
   //Add the current locus to the array of loci
   _loci.Add(curLocus); 
}
//------------------------------------------------------------------------------
//Modified for recombination
//mix of different data types (ex SNP linked to MICROSAT)
void
TDeme::print_haplotypes_to_locus(const TDataTypeArray* data_type_perLoc) {

   TLocus curLocus(_sample_size);
 
   char DNA_letters[4]={'A','G','C','T'};
   
   //int len=OriginalNodes[0]->seq_length;
   int len=ChainOriginalNodes[0]->seq_length;

   

   for (int i=0; i<_sample_size; ++i) {
      //TNode* CurNode=OriginalNodes[i];
      TNode* CurNode= (TNode*) ChainOriginalNodes[i];

      
      strstream curHaplStream;
      my_string curHapl;

      if (CurNode->sequence) {
      	int* cursite=CurNode->sequence->begin();
      	for (int j=0; j<len; ++j) {
            /*
            bool endBlock=false;
            if(j==(len-1)) {
            }
            else {
               if(*data_type_perLoc->elem(j)!=*data_type_perLoc->elem(j+1)) {
                  endBlock=true;
               }
            }
            */
         	if (*data_type_perLoc->elem(j)==2 /*DNA*/ ) {
   					curHaplStream << DNA_letters[*cursite++];
				//this might be it, though this is just adding/reading in
				
            } 
         	else if (*data_type_perLoc->elem(j)==3 /*SNP*/) {
   					curHaplStream << DNA_letters[*cursite++];
                  curHaplStream << " ";
            }         	
            else {
               	if (*data_type_perLoc->elem(j)==0 /*MICROSAT*/) {
                     curHaplStream << " " << (500+*cursite++) << " ";
                     }
               	else if(*data_type_perLoc->elem(j)==1 /*RFLP*/){
                  curHaplStream << " " << *cursite++ << " " ;
                  }
            }
         }
      }

      curHapl.read_line(curHaplStream);
      curLocus.data.Add(curHapl);
   }
   //Add the current locus to the array of loci
   _loci.Add(curLocus); 
}
//------------------------------------------------------------------------------
//Modified for recombination
void
TDeme::print_haplotypes_for_PAUP(ostream & os, const Mut_Type& data_type) {

   char DNA_letters[4]={'A','G','C','T'};
   
   //int len=OriginalNodes[0]->seq_length;
   int len=ChainOriginalNodes[0]->seq_length;

   for (int i=0; i<_sample_size; ++i) {
   
   	//TNode* CurNode=OriginalNodes[i];
      TNode* CurNode= (TNode*) ChainOriginalNodes[i];
      
      if (CurNode->sequence) {    
      	int* cursite=CurNode->sequence->begin();
			//Write id
      	os << CurNode->node_number << "." << (CurNode->deme+1) << "\n";
      	for (int j=0; j<len; ++j) {
         	if (data_type==DNA  || data_type==SNP)
   					//os << DNA_letters[(*CurNode->sequence)[j]];
   					os << DNA_letters[*cursite++];
               else
   					//os << (*CurNode->sequence)[j];
                  //Loro_1_3_99 : In order to avoid negative numbers
               	if (data_type==MICROSAT) os << (500+*cursite++) << " ";
               	else os << *cursite++;
         }
         os << "\n";
      }
   }
}
//------------------------------------------------------------------------------
ostream&
operator<<(ostream& os , const TDeme& D) {
	os  	<< "\n#Deme size   : " << D._deme_size
   		<< "\n#Sample size : " << D._sample_size
         << "\n#Growth rate : " << D._growth_rate
         << "\n";
   return os;
}

//------------------------------------------------------------------------------
ostream&
operator<<(ostream& os, const THistoricalEvent& HE) {
	os	<< "\n#Time             : " << HE.time
   	<< "\n#Source           : " << HE.source
      << "\n#Sink             : " << HE.sink
      << "\n#Migrants         : " << HE.migrants
      << "\n#New size         : " << HE.new_deme_size
      << "\n#New growth rate  : " << HE.new_growth_rate
      << "\n#New migr. matrix : " << HE.MigMat
      << endl;
   return os;
}
//------------------------------------------------------------------------------
istream&
operator>>(istream& is, THistoricalEvent& HE) {
	is	>> HE.time
   	>> HE.source
      >> HE.sink
      >> HE.migrants
      >> HE.new_deme_size
      >> HE.new_growth_rate
      >> HE.MigMat;
   return is;
}
//------------------------------------------------------------------------------
//==============================================================================
//------------------------------------------------------------------------------
long
TDemeCollection::count_lineages() {
   long count=0;
   for (int i=0 ;i<_num_demes; ++i) {
   	count+=Demes[i].lineages();
   }
   return count;
}
//------------------------------------------------------------------------------
 TDemeCollection&
 TDemeCollection::operator=(const TDemeCollection& DC) {
   _is_migration=DC._is_migration;
   _num_demes=DC._num_demes;
   _cur_mig_mat=DC._cur_mig_mat;
   if (Demes) delete[] Demes;
   if (DC.Demes) {
      Demes= new TDeme[DC._num_demes];
   	for (int i=0; i<DC._num_demes; ++i) {
      	Demes[i]=DC.Demes[i];;
      }
   }
   Events.Flush();
   if (DC.Events.GetItemsInContainer()) {
   	for (int i=0; i<DC.Events.GetItemsInContainer(); ++i) {
      	Events.Add(DC.Events[i]);
      }
   }
   MigrMatArray=DC.MigrMatArray;
   GeneTree=DC.GeneTree;
   if (CoalTimes) {
   	delete CoalTimes;
      CoalTimes=NULL;
   }
   if (DC.CoalTimes) {
      	try {
         	CoalTimes= new TMigrationMatrix();
         }
         catch (...) {
         	cout 	<< "TDeme::operator=(const TDeme& D): unable to allocate memory"
                  << endl;
            if (CoalTimes) {
      			delete CoalTimes;
      			CoalTimes=NULL;
      		}
            return *this;
         }
   	*CoalTimes=*DC.CoalTimes;
   }
   if (PairDiff) {
   	delete PairDiff;
      PairDiff=NULL;
   }
   if (DC.PairDiff) {
      	try {
         	PairDiff= new TMigrationMatrix();
         }
         catch (...) {
         	cout 	<< "TDeme::operator=(const TDeme& D): unable to allocate memory"
                  << endl;
            if (PairDiff) {
      			delete PairDiff;
      			PairDiff=NULL;
      		}
            return *this;
         }
   	*PairDiff=*DC.PairDiff;
   }
   if (MinCoalTimes) {
   	delete MinCoalTimes;
      MinCoalTimes=NULL;
   }
   if (DC.MinCoalTimes) {
      	try {
         	MinCoalTimes= new TMigrationMatrix();
         }
         catch (...) {
         	cout 	<< "TDeme::operator=(const TDeme& D): unable to allocate memory"
                  << endl;
   			if (MinCoalTimes) {
      			delete MinCoalTimes;
      			MinCoalTimes=NULL;
      		}
            return *this;
         }
   	*MinCoalTimes=*DC.MinCoalTimes;
   }
   _polym_sites=DC._polym_sites;
   _tot_mut=DC._tot_mut;
   return *this;
}
//------------------------------------------------------------------------------
//Check for migrations and update coalescence probabilities within each population
int
TDemeCollection::check_migrations(){
   _is_migration=false;
   //Examine all populations in turns to know which
   //one experience migrations. !! Do not stop the checking after first migration found
   for (int i=0; i<_num_demes; ++i)  {
   	if (Demes[i].check_migrations()) {
      	_is_migration=true;
      }
   }
   if (_is_migration) return 1;
   return 0;
}
//------------------------------------------------------------------------------
int
TDemeCollection::check_growth(){
   for (int i=0; i<_num_demes; ++i)  {
   	if (Demes[i].growth()!=0.0) {
         return 1;
      }
   }
   return 0;
}
//------------------------------------------------------------------------------
int
TDemeCollection::create_demes(const TMigrMatArray* const MA) {


	MigrMatArray=(TMigrMatArray*)MA;
   
	int size=(*MA)[0].size(); //Get number of demes from migration matrix

   TDeme aDeme;

   //Copies the pointer to the migration matrix in each deme
   aDeme.set_migration_matrix(&(*MigrMatArray)[0]);

   if (_num_demes!=size) {
   	if (Demes) delete[] Demes;
      _num_demes=size;
      Demes= new TDeme[_num_demes];
   }

   for (int i=0; i<size; ++i) {
   	//Set different ids for each deme
   	aDeme.set_id(i);
   	Demes[i]=aDeme;
   }
   return _num_demes;
}
/*
//------------------------------------------------------------------------------
int
TDemeCollection::create_lineages(const TSampSizeArray& SS) {

	initialize_samp_size(SS);

   int size=SS.GetItemsInContainer();
   if (size!=_num_demes) return 0;

   //Count total initial number of lineages (sum of sample sizes)
   long num_lineages=count_lineages();

   //Allocate room for the tree structure
   if (num_lineages) GeneTree.allocate_tree(num_lineages);
   else {
   	cout << "TDemeCollection::create_lineages() : Tree is empty " << endl;
      return 0;
   }

   //Copy the pointers to the nodes of the genealogy in each deme
	for (int i=0, current_lineage=0; i<size; ++i) {
   	Demes[i].copy_node_pointers(&GeneTree,current_lineage,num_lineages);
      current_lineage+=Demes[i].lineages();
   }
 	return 1;
}
*/
//------------------------------------------------------------------------------
//For recombination 
int
TDemeCollection::create_lineages(const TSampSizeArray& SS) {

   initialize_samp_size(SS);

   int size=SS.GetItemsInContainer();
   if (size!=_num_demes) return 0;

   //Count total initial number of lineages (sum of sample sizes)
   long num_lineages=count_lineages();

   //Allocate room for the tree structure
   if (num_lineages) GeneTree.allocate_tree(num_lineages);
   else {
   	cout << "TDemeCollection::create_lineages() : Tree is empty " << endl;
      return 0;
   }

   //Copy the pointers to the nodes of the genealogy in each deme
   TNode * curNode=GeneTree(0);
   GeneTree.chainTree.resetIterator();
   for (int i=0; i<size; ++i) {
        //Wilson
        if(Demes[i].sample_size()){

                Demes[i].copy_node_pointers(curNode, &GeneTree);
                GeneTree.chainTree.next();
                curNode=GeneTree.chainTree.get_curItem()->elem;

        }
        else{
        
        }

   }
   return 1;
}
//------------------------------------------------------------------------------
int
TDemeCollection::initialize_growth_rates(const TGrowthArray& GR) {
	int size=GR.GetItemsInContainer(); //Get number of demes from migration matrix
   if (size!=_num_demes) return 0;
  for (int i=0; i<size; ++i) {
   	Demes[i].set_growth_rate(GR[i]);
//      cout << i << "\t" << Demes[i].growth() << "\n";
   }
   return 1;
}
//------------------------------------------------------------------------------
int
TDemeCollection::initialize_deme_size(const TDemeSizeArray& DS) {
	int size=DS.GetItemsInContainer();
   if (size!=_num_demes) return 0;
   for (int i=0; i<size; ++i) {
   	Demes[i].set_deme_size(DS[i]);
   }
   return 1;
}
//------------------------------------------------------------------------------
int
TDemeCollection::initialize_samp_size(const TSampSizeArray& SS) {
	int size=SS.GetItemsInContainer();
   if (size!=_num_demes) return 0;
   for (int i=0; i<size; ++i) {
   	Demes[i].set_sample_size(SS[i]);
   }
   return 1;
}
//------------------------------------------------------------------------------
extern ofstream ArlSimulConditions;
int
TDemeCollection::initialize_events(const TEventArray& HE) {
	int size=HE.GetItemsInContainer();
   for (int i=0; i<size; ++i) {
   	Events.Add(HE[i]); 
      #ifdef _VERBOSE_
      cout << Events[i] << endl;
      #endif
   }
   return 1;
}

//------------------------------------------------------------------------------
int
TDemeCollection::implement_snp_event(const int& cur_event) {
	THistoricalEvent& CurEvent=Events[cur_event];
	if (Demes[CurEvent.source].deme_size() && CurEvent.migrants == ind_loci) { //if there are any in this deme now and this is the right chromosome
		
		
		//cout << "SNP specified time " <<  CurEvent.time  << " Source deme: " << CurEvent.source << " locus: " << CurEvent.new_deme_size <<  endl;
		
		
		
		
		
		
		
		
		TSNPEvent NewSNPEvent;
		NewSNPEvent.chrom = CurEvent.migrants;
		NewSNPEvent.locus = CurEvent.new_deme_size;
		NewSNPEvent.time = CurEvent.time;
		NewSNPEvent.deme = CurEvent.source;
		NewSNPEvent.eventnum = cur_event;
		
		//cout << "and this is event " << NewSNPEvent.eventnum << endl;
		
		
		SNPEvents.Add(NewSNPEvent);	
	}
	return 1;
}

//------------------------------------------------------------------------------
int
TDemeCollection::implement_event(const int& cur_event) {
   THistoricalEvent& CurEvent=Events[cur_event];
   
   
   
   
   //Historical Event to fix a mutation in a deme MJJ 11/08/08
   
	if(CurEvent.sink == -1) {
	
		//cout << "ind_loci inside implement event: " << ind_loci << endl;
		
		/*
		if (Demes[CurEvent.source].deme_size() && CurEvent.migrants == ind_loci) { //if there are any in this deme now and this is the right chromosome
		
		
			//cout << "SNP specified time " <<  CurEvent.time  << " Source deme: " << CurEvent.source << " locus: " << CurEvent.new_deme_size <<  endl;
			

			
		

		
			
				
				TSNPEvent NewSNPEvent;
				NewSNPEvent.chrom = CurEvent.migrants;
				NewSNPEvent.locus = CurEvent.new_deme_size;
				NewSNPEvent.time = CurEvent.time;
				NewSNPEvent.deme = CurEvent.source;
				
				
				SNPEvents.Add(NewSNPEvent);		
		
		//END if 
		}
		*/
		
		}
   
   
   
   // END MJJ 11/08/08
   
   else {
   
   //MJJ 11/25/08
   //cout << "Time: " << CurEvent.time << " lineages in " << CurEvent.source << " " << Demes[CurEvent.source].lineages() << " lineages in " << CurEvent.sink << " " << Demes[CurEvent.sink].lineages() << endl;

   
   //Step 1: Send migrants from source to sink
   // Check also that the sink population has a deme size >0
   if (CurEvent.migrants>0.0 && Demes[CurEvent.sink].deme_size()) {
   //int num_migrants= (int) (Demes[CurEvent.source].lineages()*CurEvent.migrants);
   //Loro_07_06_04 BBBB UUUUU GGGGG BUG BUG BUG BUG
        int numLineages=Demes[CurEvent.source].lineages();
        int num_migrants=0;
        long lidum=1;
        for (int i=0; i<numLineages; ++i) {
                if (ran3(&lidum)<CurEvent.migrants) ++num_migrants;
        }
        for (int i=0; i<num_migrants; ++i) {
                Demes[CurEvent.source].migrate(Demes[CurEvent.sink]);
        }
   }
   //Step 2: Resize sink deme size
   Demes[CurEvent.sink].linear_resize(CurEvent.new_deme_size);
   //Step 3: Readjust sink growth rate
   Demes[CurEvent.sink].set_growth_rate(CurEvent.new_growth_rate);
   //Step 4: Set new migration matrix
   if (_cur_mig_mat!=CurEvent.MigMat) {
   	TDeme aDeme;
      //Copies the static pointer to the migration matrix in each deme
   	aDeme.set_migration_matrix(&(*MigrMatArray)[CurEvent.MigMat]);
      _cur_mig_mat=CurEvent.MigMat;
   }
   
   }
   return 1;
}







//------------------------------------------------------------------------------
int
TDemeCollection::build_tree() {

	//MJJ
	//bool notten = true;

   long num_lineages=count_lineages();                       // sum over demes
   int next_event=0, num_event=Events.GetItemsInContainer();

   //Loro_21_9_99 :
   //By default take the first migration matrix of the array: zero migration matrix
   // if no migration matrix is provided
   TDeme aDeme;
   //Copies the static pointer to the migration matrix in each deme
   aDeme.set_migration_matrix(&(*MigrMatArray)[0]);
   _cur_mig_mat=0;

   //Check that there are migrations in the genealogical process
   bool 	implement_migration=check_migrations(),
   		implement_growth=check_growth(),
         implement_recombination;

   if (num_linked_loci<2) implement_recombination=false;
   else {
      if (rec_rate>1.0e-10) implement_recombination=true;
      else implement_recombination=false;
   }

   long Sum_DeathVect=0;
   GeneTree.num_linked_loci=num_linked_loci;

   if (GeneTree.DeathVect) delete [] GeneTree.DeathVect;
   if (GeneTree.MRCA_list) delete [] GeneTree.MRCA_list;
   GeneTree.DeathVect = new int[num_linked_loci];
   GeneTree.MRCA_list = new TNode* [num_linked_loci];
      //GeneTree.tree_length_perLoc = new double[num_linked_loci];
      //GeneTree.tot_mut_perLoc = new double[num_linked_loci];

   for(int i=0; i<num_linked_loci; i++){
      GeneTree.DeathVect[i]=count_lineages();
      GeneTree.MRCA_list[i]=NULL; 
      //GeneTree.tree_length_perLoc[i]=0;
      //GeneTree.tot_mut_perLoc[i]=0;
   }

   Sum_DeathVect=num_linked_loci*num_lineages;
	
	
	//MJJ 1/29/09
	//Check for SNP Historical events
	for(int isnp = 0; isnp <num_event; isnp++){

		if(Events[isnp].sink == -1){
			
			
         	implement_snp_event(isnp);
		}
		
	}

	

   //We'll simulate a maximum of MAX_GENERATIONS generations
   long time;
   for (time=1 ; time<MAX_GENERATIONS && num_lineages>1 &&
                           Sum_DeathVect>num_linked_loci; ++time) {


      
      //if(time==100){
      //  int tempstop=1;
      //}

	//MJJ
	//if(count_lineages() == 25 && notten) {
	//if (time == 5250 && notten) {
//	notten = false;
	//cout << "Time: " << time << " Total Lin: " << count_lineages() << endl;
	//}

      //Step 1: Check for possible historical events
      if (next_event<num_event) {
      	while (Events[next_event].time==time) {
			//if(Events[next_event].sink == -1) numdonesnpevents++;

			//cout << "event time: " << time << endl;
			
         	implement_event(next_event);
         	//Check that the migration matrix has not changed
         	implement_migration=check_migrations();
         	//Check that the growth rates have not changed
         	implement_growth=check_growth();
            //Loro_21_9_99: To prevent the comparison of an unexisting event !
         	if (++next_event==num_event) break;
         }
      }

		//Step 2: Migration round
      if (implement_migration) {
      	for (int j=0; j<_num_demes; ++j) {
            if (Demes[j].is_migration()) {
               for (int k=0; k<_num_demes; ++k) {
         			if (j!=k) Demes[j].send_migrants(Demes[k],k);
               }
            }
      	}
      }

   	//Step 3: Coalescence round
      for (int j=0; j<_num_demes; ++j) {
         #ifdef _MCE_
            Demes[j].multipleCoalescentEvents(time, GeneTree, this);
         #else
            #ifdef _MIXT_
               if(Demes[j].lineages()*(Demes[j].lineages()-1)*0.5*(1./Demes[j].get_deme_size()) < 1.) {
                  Demes[j].singleCoalescentEvent(time, GeneTree, this);
               }
               else {
                  Demes[j].multipleCoalescentEvents(time, GeneTree, this);
               }
            #else
               Demes[j].singleCoalescentEvent(time, GeneTree, this);
            #endif
         #endif
      }

      /*************************************************************************

            Call a new function to compute recombination for each TNode

      *************************************************************************/
      //Step 4: Recombination round
      if (implement_recombination) {
         for (int j=0; j<_num_demes; ++j) {
            #ifdef _MRE_
               Demes[j].recombination_per_lineage(/*coalNodeBlocked[j],*/ time, GeneTree, this);
            #else
               Demes[j].One_recombination_per_generation(/*coalNodeBlocked[j],*/ time, GeneTree, this);
            #endif
         }
      }

      //Loro_03_03_04  Not very efficient... could be optimized...
      //Implement a static number at the TDemeCollection level
      num_lineages=0;
      for(int j=0; j<_num_demes; ++j) {
         num_lineages+=Demes[j].lineages();
      }
      //Debug
      //cout << "num_lineages: " << num_lineages;
      /*************************************************************************/

   	//Step 5: Resize deme size after growth for next generation


      
      //long bidon2=Demes[0].deme_size();
      //if( bidon2 <= 50000 ){
      // time;
      // int stopped=0;
      //}



      if(time == 8060){
        long bidon=Demes[0].deme_size();
        bidon;
        int stopped=0;
      }


      if (implement_growth)
      for (int j=0; j<_num_demes; ++j) {
      	if (Demes[j].growth()) {
            if (!Demes[j].exponential_resize()) {
            	cout 	<< "TDemeCollection::build_tree(): size of deme "
               		<< j << "is zero !" << endl;
               return 0;
            }
         }
      }


      
      //Loro_03_03_04  Not very efficient... could be optimized...
      //Also implement a static member at the TDemeCollection level
      Sum_DeathVect=0;
      for (int loc=0; loc<num_linked_loci; ++loc) {
           Sum_DeathVect+=GeneTree.DeathVect[loc];
      }

/*
	   cout << "Time: " << time << " ";
	   for(int j=0; j<_num_demes; ++j) {
		   cout  << "deme" << j << ": " << Demes[j].lineages() << " ";
		   if (Demes[j].lineages() > 0) {
			   Demes[j].oldestLineage = time;
		   }
	   }
	   cout << endl;
	   */
	

      //if(Demes[0].lineages()){
      //  int tempstop=1;
      //}

   }

		   //cout << "end gen: " << time << "   #lineages: " << num_lineages <<endl;
   if (!rec_rate){
      if (num_lineages>1) {
   	   cout 	<< "TDemeCollection:: build_tree(): Coalescence process did not converge" << endl;
         return 0;
      }
   }
   else {
      if (num_lineages>1) {
         return 1;
      }
   }
   
   return 1;
}
//------------------------------------------------------------------------------
void
TDemeCollection::reset(const TDemeSizeArray& SA, const TGrowthArray& GA) {

   initialize_deme_size(SA);
   initialize_growth_rates(GA);

   int SampleSize=GeneTree.sample_size();
   GeneTree.chainTree.Call_deleteChain(SampleSize);
   
   for (int i=0; i<_num_demes; ++i) {
        Demes[i].reset();
   }

}
//------------------------------------------------------------------------------
int
TDemeCollection::compute_moments_of_demes_coalescence_times() {

   //Get the root node
	//TNode& root=GeneTree[2*GeneTree.sample_size()-2];     // old version without recom
   TNode& root=*GeneTree.chainTree.last();      
   

   //Build the list of descendent nodes
   root.build_lists_of_descendent_nodes();

   if (!CoalTimes) {
   	try {
   		CoalTimes=new TMigrationMatrix(_num_demes,ISLAND,0.0);
   	}
      catch (...) {
      	if (CoalTimes) delete CoalTimes; CoalTimes=NULL;
         cout  << "TDemeCollection::compute_moments_of_demes_coalescence_times()"
         		<< ": unable to allocate memory"
               << endl;
         return 0;
      }
   }
   else CoalTimes->set_elems(0.0);

  	if (!PairDiff) {
   	try {
   		PairDiff=new TMigrationMatrix(_num_demes,ISLAND,0.0);
   	}
      catch (...) {
      	if (PairDiff) delete PairDiff; PairDiff=NULL;
         cout  << "TDemeCollection::compute_moments_of_demes_coalescence_times()"
         		<< ": unable to allocate memory"
               << endl;
         return 0;
      }
   }
   else PairDiff->set_elems(0.0);

   if (!MinCoalTimes) {
   	try {
   		MinCoalTimes=new TMigrationMatrix(_num_demes,ISLAND,0.0);
   	}
      catch (...) {
      	if (MinCoalTimes) delete MinCoalTimes; MinCoalTimes=NULL;
         cout  << "TDemeCollection::compute_moments_of_demes_coalescence_times()"
         		<< ": unable to allocate memory"
               << endl;
         return 0;
      }
   }
   else MinCoalTimes->set_elems(MAXLONG);

   //Compute total coalescent times within and among demes
   root.compute_total_coal_times_among_demes(*this,*CoalTimes, *PairDiff, *MinCoalTimes);

   //Compute mean coalescent times
   for (int i=0; i<_num_demes; ++i) {
      long sampsize1=Demes[i].sample_size();
   	//Loro_20_9_99
      if (sampsize1)
      //Among pop treatment
   	for (int j=0; j<i; ++j) {     
      	long sampsize2=Demes[j].sample_size();
      	if (sampsize2) {
      		(*CoalTimes)(i,j)/=sampsize1*Demes[j].sample_size();
      		(*PairDiff)(i,j)/=sampsize1*Demes[j].sample_size();
         }
      }
      //Within pop treatment
      if (sampsize1>1L) {
      	(*CoalTimes)(i,i)/=sampsize1*(sampsize1-1)*0.5;
      	(*PairDiff)(i,i)/=sampsize1*(sampsize1-1)*0.5;
      }
   }
	return 1;
}
//------------------------------------------------------------------------------
//A procedure to sprinkle mutations once the tree branch lengths are known
//To simulate a mix between different data types
//Impossible to create a table of enum variable ????
const long&
TDemeCollection::sprinkle_mutations(const double&           mut_rate,
                                    const int&              numLinkedLoci,
                                    const TDataTypeArray*   data_type,
                                    const double&           gamma_par,
                                    const double&           mut_ratio,
                                    const double&           prop_sites,
                                    const double&           trans_rate,
                                    const int&              range_constr,
									ofstream &snpfs) {






   //Loro24_03_03: The next few lines may be unnecessary????!!!! (Remove...)
   TNode* root0=GeneTree.MRCA_list[0];
   int size=root0->hits.GetItemsInContainer(), to_add=numLinkedLoci-size;
   //Loro_29_8_98 //Resize sequence length if needed
   for (int i=0; i<to_add; ++i){
      root0->hits.Add(0);
   }
   for (int i=0; i<numLinkedLoci; ++i){
      root0->hits[i]=0;
   }
double WholeRecRate=0;
//double *RecRatePerLoc=_rec_rate_perLoc->begin();

//Recombination: different mutation rate between linkage blocks
//Reset the current mutation rate.
int Homogenous=1;
for(int i=0;i<numLinkedLoci;i++) {
        WholeRecRate+=*_rec_rate_perLoc->elem(i);
        if(i > 0){
                if(*_mut_rate_perLoc->elem(i-1) != *_mut_rate_perLoc->elem(i) ){
                        Homogenous = 0;
                }
                if(*_geom_param_perLoc->elem(i-1) != *_geom_param_perLoc->elem(i) ){
                        Homogenous = 0;
                }
                if(*_transRatePerLoc->elem(i-1) != *_transRatePerLoc->elem(i) ){
                        Homogenous = 0;
                } 
                if(*_rangeConstraintPerLoc->elem(i-1) != *_rangeConstraintPerLoc->elem(i) ){
                        Homogenous = 0;
                }
        }
}

//Debug
//WholeRecRate=1;
	/*
	cout << "Make oldest lineages" << endl;
	//MJJ 1/29/10 create array
	int oldestLineages [_num_demes];
	for (int i=0 ;i<_num_demes; ++i) {
		oldestLineages[i]=Demes[i].oldestLineage;
		cout << i << " " << oldestLineages[i] << endl;
	}
	*/

	//MJJ 1/29/10
	

//Guillaume 05 04 2004 Test if global recombination of the chromosome segemnt is equal to 0
//AND block homogeneous block
if( (WholeRecRate  <= 0.000000000001) && (Homogenous) ){ //Sprinkle mutation for the whole chromsome segment


//MJJ 11/13/08
	int SNPtime = -1;
	int SNPdeme = -1;
	int SNPeventnum = -1;
	



        //MICROSAT in the data
        bool MICROSATinDATA=false;
        int number_MICROSAT=0;
        for(int i=0;i<numLinkedLoci;i++){
             if(*data_type->elem(i)==0){
                MICROSATinDATA=true;
                number_MICROSAT++;
             }
        }
        if (MICROSATinDATA) {
                int    rangeConstrPerLocus=*_rangeConstraintPerLoc->elem(0);
                if (rangeConstrPerLocus) {//Fix minimum and maximum size for microsats
                        double y=2.0, x=rangeConstrPerLocus;
                        if (fmod(x,y)==0.0) { //Even number
            	                root0->min_mic=-(rangeConstrPerLocus/2)+1;
                                root0->max_mic=rangeConstrPerLocus/2;
                        }
                        else {
            	                root0->min_mic=(-rangeConstrPerLocus+1)/2;
                                root0->max_mic=(rangeConstrPerLocus-1)/2;
                        }
                }
                else {
                        root0->min_mic=0;
                        root0->max_mic=0;
                }
        }


        //SNPs in the data
        bool SNPinDATA=false;
        //Initialize hits to zero for this site
        root0->reset_tree_length();
        //Any positive number will not reset the random number generator
        long my_lidum=1L;
        //Create a vector of the number of mutation per locus
        int *num_mut_perLoc= new int[numLinkedLoci];


        int number_SNPs=0;
        for(int i=0;i<numLinkedLoci;i++){
             num_mut_perLoc[i]=0;
             if(*data_type->elem(i)==3){
                SNPinDATA=true;
                number_SNPs++;
             }
        }
        double *Freq_min_par_SNP= new double[numLinkedLoci];
        double *Freq_SNP= new double[numLinkedLoci];
        for(int i=0;i<numLinkedLoci;i++){
             if(*data_type->elem(i)==3){
                int k=get_lociBlock(i);
                Freq_min_par_SNP[i]=freq_SNP_min[k];
                Freq_SNP[i]=0;
             }
             else{
                Freq_min_par_SNP[i]=0;
             }
        }

        bool *dejamut=new bool[numLinkedLoci];
        if(SNPinDATA){
                int sumTreebranch=0.0;
                int *sizeToNode= new int[(2*sum_sample_size)];
                int *SNP_mut_pos= new int[numLinkedLoci];
                TNode **scannedNodes= new TNode*[(2*sum_sample_size)];
                TNode **rootSNPs= new TNode*[numLinkedLoci];

                int IdNode=0;
                int* pIdNode=&IdNode;
                //*pIdNode=5;
                sizeToNode[0]=0;
                scannedNodes[0]=root0;
                root0->compute_TreeSize_SizeToNode_NumDescList( pIdNode,
                                                                sumTreebranch,
                                                                sizeToNode,
                                                                scannedNodes);

                int tempverif=0;
                for(int i=0;i<(2*sum_sample_size-1);i++){
                        tempverif+=sizeToNode[i];
                }


                bool redo=false;
                int stopRedo=0;

                for(int i=0;i<numLinkedLoci;i++){           //

                dejamut[i]=false;
                if(*data_type->elem(i)==3){
                        //dejamut[i]=false;
                        if(redo){
                                i--;
                        }
                        redo=false;
					
						//MJJ 11/08/08 shut-off for SNP mutation if freq min snp = -1	
					int k=get_lociBlock(i);
					//cout << "freq_SNP_min for locus " << k << " is " << freq_SNP_min[k] << endl;
						 if(freq_SNP_min[k] == -1 ){ //MJJ 11/09/08 Shuts off SNP mutation if freq set to -1 exactly
							SNP_mut_pos[i] = -1;
							
						}
						else {
							SNP_mut_pos[i]= ran3(&my_lidum)*(sumTreebranch);
						}
						
						//cout << "And thus SNP_mut_pos set to " << SNP_mut_pos[i] << endl;
						
                        int branch=0;
                        for(int j=0;j< ( 2*sum_sample_size );j++){
                                double Verif_Freq_SNP;
                                branch=sizeToNode[j];
                                if(branch>SNP_mut_pos[i]){
                                        Freq_SNP[i]=scannedNodes[j]->count_SNP_desc(0)/(double)sum_sample_size;
                                        if (Freq_min_par_SNP[i] > 0 ) {

                                                if (Freq_SNP[i] < Freq_min_par_SNP[i]) {
                                                        redo=true;
                                                        stopRedo++;
                                                        if(stopRedo>1000){
                                                                redo=false;
                                                        }
                                                        rootSNPs[i]=scannedNodes[j];
                                                        break;
                                                }
                                        }
										

										
                                        else{
                                                double temp_freq=Freq_min_par_SNP[i];
                                                temp_freq=(-1 - temp_freq);
                                                temp_freq=(1 + temp_freq);
                                                if (Freq_SNP[i] != temp_freq ) {
                                                        redo=true;
                                                        stopRedo++;
                                                        if(stopRedo>1000){
                                                                redo=false;
                                                        }
                                                        rootSNPs[i]=scannedNodes[j];
                                                        break;
                                                 }

                                        }
                                        rootSNPs[i]=scannedNodes[j];
                                        break;
                                }
                        }
                }
                else{
                        SNP_mut_pos[i]=0;
                        rootSNPs[i]=NULL;
                }
                }          
		int zero=0;
                root0->add_mutations_without_recombination(&my_lidum,
                                                           num_mut_perLoc,
                                                           zero,
                                                           //1,
                                                           numLinkedLoci,
                                                           data_type,
                                                           _mut_rate_perLoc,
                                                           _geom_param_perLoc,
                                                           gamma_par,
                                                           _transRatePerLoc,
                                                           _rangeConstraintPerLoc,
                                                           rootSNPs,
                                                           dejamut,
                                                           NULL, SNPtime, SNPdeme, snpfs, SNPeventnum);

                if(sizeToNode){
                        delete[] sizeToNode;
                }
                if(SNP_mut_pos){
                        delete[] SNP_mut_pos;
                }
                if(scannedNodes){
                        delete[] scannedNodes;
                }
                if(rootSNPs){
                        delete[] rootSNPs;
                }
         }
         else{
		int zero=0;
                root0->add_mutations_without_recombination(&my_lidum,
                                                           num_mut_perLoc,
                                                           zero,
                                                           //1,
                                                           numLinkedLoci,
                                                           data_type,
                                                           _mut_rate_perLoc,
                                                           _geom_param_perLoc,
                                                           gamma_par,
                                                           _transRatePerLoc,
                                                           _rangeConstraintPerLoc,
                                                           NULL,
                                                           dejamut,
                                                           NULL,  SNPtime, SNPdeme, snpfs, SNPeventnum);


        }
        if(num_mut_perLoc){
                delete[] num_mut_perLoc;
        }
        if(Freq_min_par_SNP){
                delete[] Freq_min_par_SNP;
        }
        if(Freq_SNP){
                delete[] Freq_SNP;
        }
        if(dejamut){
                delete[] dejamut;
        }

}
else { //Sprinkle mutation marker per marker

   bool redo=false;
   //Loro_01_04_04 Bug
   //for (int curLocus=0; curLocus<num_linked_loci; ++curLocus){
   int stopRedo=0;
   for (int curLocus=0; curLocus<numLinkedLoci; ++curLocus){
      if (redo){
         curLocus--;
      }
      _tot_mut=0;

      Mut_Type mut_type;
      if(*data_type->elem(curLocus)==0) {
         mut_type=MICROSAT;
      }
      else if(*data_type->elem(curLocus)==1) {
         mut_type=RFLP;
      }
      else if(*data_type->elem(curLocus)==2) {
         mut_type=DNA;
      }
      else if(*data_type->elem(curLocus)==3) {
         mut_type=SNP;
      }
   
      TNode* curRoot=GeneTree.MRCA_list[curLocus];
      //Initialize hits to zero for this site
      curRoot->reset_tree_length();
      //Any positive number will not reset the random number generator
      long my_lidum=1L;

      //Recombination: different mutation rate between linkage blocks
      //Reset the current mutation rate.
      double mutRatePerLocus=*_mut_rate_perLoc->elem(curLocus);
      double geomParPerLocus=*_geom_param_perLoc->elem(curLocus);
      double transRatePerLocus=*_transRatePerLoc->elem(curLocus);

      //Loro_31_03_04
      int    rangeConstrPerLocus=*_rangeConstraintPerLoc->elem(curLocus);
      
      //Loro_02_03_04 potential bug because each microsat loci should be checked for potential range constraint...
      if (mut_type==MICROSAT) {
         if (rangeConstrPerLocus) {//Fix minimum and maximum size for microsats
					double y=2.0, x=rangeConstrPerLocus;
            if (fmod(x,y)==0.0) { //Even number
            	curRoot->min_mic=-(rangeConstrPerLocus/2)+1;
               curRoot->max_mic=rangeConstrPerLocus/2;
            }
            else {
            	curRoot->min_mic=(-rangeConstrPerLocus+1)/2;
               curRoot->max_mic=(rangeConstrPerLocus-1)/2;
            }
         }
         else {
            curRoot->min_mic=0;
            curRoot->max_mic=0;
         }
      }

      /* TODO -olaurent -crecombinaison : reset current data type for that locus */
      //curRoot->setDatType(curDataType);

      //Call a recursive procedure to add mutation
      if(mut_type==SNP) {                                                       //SNP case
         double sumTreebranch=0.0;                                              //independant from mutation rate
         double SNP_mut_pos=0;
         TNode* SNP_root=NULL;
         int num_SNP_desc=0;
         double Freq_SNP=0;
         curRoot->compute_TreeSize(sumTreebranch,curLocus);

         while(!SNP_mut_pos) {
              SNP_mut_pos= ran3(&my_lidum)*(sumTreebranch);  // integer random value between 0 and sumtime-1
         }


         



         const int& zero=0;
	 const int& Loc=curLocus;
	 const TNode* pNULL=NULL;
	 const double& pSNP_mut_pos=SNP_mut_pos;
	 double dblZero=0.;
	 bool dejamut=false;
	 
//MJJ 11/13/08

	int SNPtime = -1;
	int SNPdeme = -1;
		  	int SNPeventnum = -1;



	 int k=get_lociBlock(curLocus);
	 
	 if(freq_SNP_min[k] == -1 ){ //MJJ 11/09/08 Shuts off SNP mutation if freq set to -1 exactly
		SNP_mut_pos = -1; //MJJ 
		
		
		//MJJ Test SNPEvents
			int num_snpevent=SNPEvents.GetItemsInContainer();
      for (int next_snpevent=0; next_snpevent<num_snpevent; next_snpevent++){ 
      	if (SNPEvents[next_snpevent].locus==curLocus) { //MJJ THERE CAN BE ONLY ONE (PER LOCUS)
			
			SNPtime = SNPEvents[next_snpevent].time;
			SNPdeme = SNPEvents[next_snpevent].deme;
			SNPeventnum = SNPEvents[next_snpevent].eventnum;
			

         }
      }




//MJJ End
		
		
		

		
		
		
		         _tot_mut+=curRoot->add_mutations(&my_lidum, 0, curLocus,
                                          //1,
                                          numLinkedLoci,
                                          mut_type,
                                          mutRatePerLocus,
                                          geomParPerLocus, 
                                          gamma_par, transRatePerLocus,
                                          rangeConstrPerLocus, dejamut, NULL, SNPtime, SNPdeme, snpfs, SNPeventnum); //MJJ Added snptime, snpdme and dejamut and oldestlineages
		
		

		
		
		}
		
		
		else{
			

			
		         _tot_mut+=curRoot->add_SNP_mutation(&my_lidum, zero, Loc,
                                             	numLinkedLoci,
                                                gamma_par, transRatePerLocus,
                                             	range_constr, pNULL, 
                                             	pSNP_mut_pos, dblZero, dejamut, SNP_root
                                                );
		
		
		
	 
/*
         _tot_mut+=curRoot->add_SNP_mutation(&my_lidum, zero, Loc,
                                             	numLinkedLoci,
                                                gamma_par, transRatePerLocus,
                                             	range_constr, pNULL, 
                                             	pSNP_mut_pos, dblZero, dejamut, SNP_root
                                                );
*/
//MJJ end

         num_SNP_desc=SNP_root->count_SNP_desc(curLocus);
         Freq_SNP=num_SNP_desc/(double)sum_sample_size;
         

         if(freq_SNP_min[k] > 0 ){

                if (Freq_SNP < freq_SNP_min[k]) {
                        //curLocus--;
                        redo=true;
                        stopRedo++;
                        if(stopRedo>1000){
                                redo=false;
                                stopRedo=0;
                        }
                }
                else {
                        redo=false;
                        stopRedo=0;
                }
         }
         else{
                double temp_freq=freq_SNP_min[k];
                temp_freq=(-1 - temp_freq);
                temp_freq=(1 + temp_freq);
                if (Freq_SNP != temp_freq) {
                        //curLocus--;
                        redo=true;
                        stopRedo++;
                        if(stopRedo>1000){
                                redo=false;
                                stopRedo=0;
                        }
                }
                else {
                        redo=false;
                        stopRedo=0;
                }

         }
		 
		 //MJJ End else put around original (non planted) SNP mutation calls
		 }
		 
      }
      else {
	  bool dejamut=false;
         _tot_mut+=curRoot->add_mutations(&my_lidum, 0, curLocus,
                                          //1,
                                          numLinkedLoci,
                                          mut_type,
                                          mutRatePerLocus,
                                          geomParPerLocus, 
                                          gamma_par, transRatePerLocus,
                                          rangeConstrPerLocus, dejamut, NULL, -1 , -1, snpfs, -1);
      }
   }

}

   return _tot_mut;


}
//------------------------------------------------------------------------------
extern ostream& operator<<(ostream& os,const TDemeCollection& DC) {

	for (int i=0, numdemes=DC.num_demes(); i<numdemes; ++i) {
   	os << "\nDeme #" << (i+1)
      	<< "\n----------\n\n"
         << DC.Demes[i]
         << "\n---------------------------------------\n";
   }

   int size=DC.MigrMatArray->GetItemsInContainer();
   for (int i=0; i<size; ++i) {
   	os << "\nMIGRATION MATRIX " << i
   		<< "\n=================\n"
      	<< (*DC.MigrMatArray)[i]
      	<< "\n";
   }

   os << "\nHISTORICAL EVENTS"
   	<< "\n=================\n";

   for (int i=0, numevents=DC.Events.GetItemsInContainer(); i<numevents; ++i) {
   	os << "\nEvent #" << (i+1)
      	<< "\n----------\n\n"
         << DC.Events[i]
         << "\n---------------------------------------\n";
	}
   return os;
}

//------------------------------------------------------------------------------
void
TDemeCollection::write_samples_to_Arlequin_file(ostream& os, const Mut_Type& data_type) {
   for (int i=0; i<_num_demes; ++i)
   if (Demes[i].sample_size())  //Loro_20_9_99
   {
   	write_Arlequin_sample_header(i+1, Demes[i].sample_size(), os, 0 /*haplotypic data*/);
      Demes[i].print_haplotypes_for_Arlequin(os, i, data_type);
      write_Arlequin_sample_footer(os);
   }
}
//------------------------------------------------------------------------------
void
TDemeCollection::write_locus_data_to_array(const Mut_Type& data_type) {
   for (int i=0; i<_num_demes; ++i)
   if (Demes[i].sample_size())  //Loro_20_9_99
   {
      Demes[i].print_haplotypes_to_locus(data_type);
   }
}
//New function to deal with a mix of different data types (ex SNP linked to MICROSAT)
//------------------------------------------------------------------------------
void
TDemeCollection::write_locus_data_to_array(const TDataTypeArray* data_type_perLoc) {
   for (int i=0; i<_num_demes; ++i)
   if (Demes[i].sample_size())  //Loro_20_9_99
   {
      Demes[i].print_haplotypes_to_locus(data_type_perLoc);
   }
}                                                                                 
//------------------------------------------------------------------------------
void
TDemeCollection::write_loci_to_Arlequin_file(ostream& os, const Mut_Type& data_type, const int& genot_data) {
   for (int i=0; i<_num_demes; ++i)
   if (Demes[i].sample_size())  //Loro_20_9_99
   {
   	write_Arlequin_sample_header(i+1, Demes[i].sample_size(), os, genot_data);
      Demes[i].print_loci_for_Arlequin(os, i, data_type, genot_data);
      write_Arlequin_sample_footer(os);
   }
}
                                                                               
//------------------------------------------------------------------------------
//add for recombination 
void
TDemeCollection::write_loci_to_Arlequin_file(ostream&          os, 
                                             const Mut_Type&   data_type, 
                                             const int&        genot_data,
                                             const int&        num_indep_loci) {
   for (int i=0; i<_num_demes; ++i)
   if (Demes[i].sample_size())  //Loro_20_9_99
   {  
      TDeme* curDeme=&Demes[i];
      os
      << "#\t\tNumber of recombinations and number of effective recombinations: " 
      << number_recombHits_perDeme(curDeme)
      << ", " 
      << (number_recombHits_perDeme(curDeme)-number_no_effective_recombHits_perDeme(curDeme))
      << "\n";
      for(int indep_loci=0;indep_loci<num_indep_loci;indep_loci++) {
      os
      << "#\t\tRecombinations hits between adjacent loci: \n"
      << "#\t\t";
         //int temp= Demes[i].get_id();
         for(int cpt1=0;cpt1<(num_linked_loci-1);cpt1++) {
            os << return_recombHits(curDeme->get_id()*num_indep_loci*(num_linked_loci-1)
                                    + (indep_loci*(num_linked_loci-1) + cpt1)) << " ";
         }
      os << ": Total " << number_recombHits_perDeme_perIndepLoci(curDeme,indep_loci) << "\n";
      os
      << "#\t\tEffective recombinations hits between adjacent loci: \n"
      << "#\t\t";
         //temp= Demes[i].get_id();
         for(int cpt1=0;cpt1<(num_linked_loci-1);cpt1++) {
            os << ( return_recombHits(curDeme->get_id()*num_indep_loci*(num_linked_loci-1)
                                    + (indep_loci*(num_linked_loci-1) + cpt1))
                  - return_no_effective_recombHits(curDeme->get_id()*num_indep_loci*(num_linked_loci-1)
                                    + (indep_loci*(num_linked_loci-1) + cpt1)) ) << " ";
         }  
      os << ": Total " << ( number_recombHits_perDeme_perIndepLoci(curDeme,indep_loci) 
                            - number_no_effective_recombHits_perDeme_perIndepLoci(curDeme,indep_loci) )
                       << "\n";
      }
       
   	write_Arlequin_sample_header(i+1, Demes[i].sample_size(), os, genot_data);
      Demes[i].print_loci_for_Arlequin(os, i, data_type, genot_data);
      write_Arlequin_sample_footer(os);
   }
}
//------------------------------------------------------------------------------
void
TDemeCollection::write_group_section_to_Arlequin_file(ostream& os) {
	os
   	<< "\n[[Structure]]\n"
      << "\n\tStructureName=\"Simulated data\""
      << "\n\tNbGroups=1"
      << "\n\tGroup={";
   for (int i=0; i<_num_demes; ++i) {
   	if (Demes[i].sample_size())   //Loro_20_9_99
   	os
      	<< "\n\t   \"Sample " << (i+1) << "\"";
   }
   os
   	<< "\n\t}" << endl;
}
//------------------------------------------------------------------------------
void
TDemeCollection::write_samples_to_PAUP_file(ostream& os, const Mut_Type& data_type) {
   for (int i=0; i<_num_demes; ++i) {
		Demes[i].print_haplotypes_for_PAUP(os, data_type);
   }
}
//------------------------------------------------------------------------------
void
TDemeCollection::flushLoci() {
   for (int i=0; i<_num_demes; ++i) {
		Demes[i].flushLoci();
   }
}
//------------------------------------------------------------------------------
void
TDemeCollection::resetRecCounts() {
   for(int cpt=0;cpt<_num_demes*num_indep_loci*(num_linked_loci-1);cpt++) {
      recombHits[cpt]=0;
   }
}      
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------





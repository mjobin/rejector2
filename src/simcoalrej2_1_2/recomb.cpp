//---------------------------------------------------------------------------

#pragma hdrstop

//#include "recomb.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)

#include "genealgy.h"
#include "deme.h"

using namespace std;

/* TODO 5 -olaurent -crecombinaison : Implement a mapping function, in order to avoid additivity of recombination rates */

//------------------------------------------------------------------------------
// Implements the recombination process process within population
// For every Node implement the recombination with a Probability
// rec_proba = 2(min_pos - max_pos)*rec_rate
int
TDeme::recombination_per_lineage(const long& time, TTree& GeneTree, TDemeCollection* myDemes) {

   //If there is only one lineage left, do not even try to recombinate it
   //(it is posible but this shows the stop of the simulation run)
   if (_lineages<2) {
      return 0;
   }
   if (!_deme_size) {
      cout 	<< "TDeme::recombination_per_lineage(CoalNode,time,GeneTree,pos) : deme size of population "
      		<< _id << " is zero !" << endl;
   	return 0;
   }

   #ifdef _COMPUTE_MOMENTS_
   TDeme* curDeme=this;
   #endif

   TNode *Casc1, *Casc2;
   ChainNodeList.resetIterator();
   TNode * currentNode= (TNode*) ChainNodeList[0];

   //Debug
   //int NumRecombPerGeneration=0;

   for(int cpt=0; cpt < _lineages; ++cpt, currentNode=ChainNodeList.next()){
      int min_pos = currentNode->min_pos;                                       // Remember that the position begin to 0
      int max_pos = currentNode->max_pos;                                       // and that max_pos is the Number






      //Loro_01_03_04
      if (min_pos==max_pos) continue;
                                                                                // of linked loci minus 1
      //double recRate = myDemes->get_rec_rate();
      double rec_proba=currentNode->update_recombination_proba(min_pos,         // Decide if the recombination occurs
                                                               max_pos,
                                                               myDemes->get_rec_rate_perLoc());
      /*
      //Loro_01_03_04  !!!! wrong, to be replaCED BY A MAPPING FUNCTION...
      if (rec_proba>0.5) rec_proba=0.5;
      */
	  
	  

	  
	  

      double temp =ran3(&_lidum);
      if (temp < rec_proba && min_pos < max_pos){                                                    // Decide where the recombination occurs
      //Debug
      //NumRecombPerGeneration++;
         //int recomb_pos = currentNode->recombinationPosition(min_pos,max_pos);// Uniform recombination rate

         int recomb_pos = currentNode->recombinationPosition(min_pos,           // different recombination rate between loci
                                                             max_pos,
                                                             myDemes->get_rec_rate_perLoc(),
                                                             rec_proba);
			




         #ifdef _COMPUTE_MOMENTS_
            myDemes->add_recombHits(curDeme,recomb_pos);                        //Compute between loci distribution of the recombination Hits
            if(( myDemes->GeneTree.MRCA_list[recomb_pos] && myDemes->GeneTree.MRCA_list[recomb_pos+1])
            || ( !currentNode->flag[recomb_pos] || !currentNode->flag[recomb_pos + 1]) ) {
               myDemes->add_no_effective_recombHits(curDeme,recomb_pos);
            }
         #endif

         Casc1 = new TNode();
         Casc1->min_pos=currentNode->min_pos;                                   // create the first recombinating TNode
         for (int f=recomb_pos; f>=Casc1->min_pos; --f) {                       // Update the max_pos of asc1
            if (currentNode->flag[f]) {
               Casc1->max_pos=f;
               break;
            }
         }
         Casc1->time=time;
		 
		 //MJJ
		Casc1->deme=_id;
		 

		 		 
         Casc1->desc1=currentNode;
         Casc1->desc2=NULL;

         for (int cpt=0; cpt<=recomb_pos; ++cpt){                               // Update its flag array
            Casc1->flag[cpt]=currentNode->flag[cpt];
         }
         for (int cpt=recomb_pos+1; cpt<Casc1->num_linked_loci; ++cpt){
            Casc1->flag[cpt]=false;
         }

         Casc2 = new TNode();                                                   // create the second recombining TNode
         Casc2->max_pos=currentNode->max_pos;                                   // Update the min_pos of asc2
         for (int f=recomb_pos+1; f<=Casc2->max_pos; ++f) {
            if (currentNode->flag[f]) {
               Casc2->min_pos=f;
               break;
            }
         }
         Casc2->time=time;
		 
		 //MJJ
		Casc2->deme=_id;

		 
         Casc2->desc1=currentNode;
         Casc2->desc2=NULL;

         // Update its flag array
         for (int cpt=0; cpt<=recomb_pos; cpt++){
            Casc2->flag[cpt]=false;
         }
         for (int cpt=recomb_pos+1; cpt<Casc2->num_linked_loci; ++cpt) {
            Casc2->flag[cpt]=currentNode->flag[cpt];
         }

         currentNode->ancestor=Casc1;
         currentNode->asc1=Casc1;
         currentNode->asc2=Casc2;

         // add asc1 and asc2 in the Chained List chainTree
         GeneTree.chainTree.Add(Casc1);
         Casc1->ID_Node=(GeneTree.chainTree.get_size() - 1 );
         Casc1->event = 2;

         GeneTree.chainTree.Add(Casc2);
         Casc2->ID_Node=(GeneTree.chainTree.get_size() - 1 );
         Casc2->event = 2;

         ChainNodeList.replace_elem(Casc1);
         ChainNodeList.Add(Casc2);
		 

         ++_lineages;
       }
   }
   //Debug
   //cout << "Time: " << time << " Num Recomb:" << NumRecombPerGeneration << endl;
  // cout << "Time: " << time << endl;
   return 1;
}
//------------------------------------------------------------------------------
// HUDSON assumption: 1 recombination per generation
int
TDeme::One_recombination_per_generation(const long& time, TTree& GeneTree, TDemeCollection* myDemes) {

   //If there is only one lineage left, do not even try to recombinate it
   //(it is posible but this shows the stop of the simulation run)
   if (_lineages<2) {
      return 0;
   }
   if (!_deme_size) {
      cout 	<< "TDeme::recombination_per_lineage(CoalNode,time,GeneTree,pos) : deme size of population "
      		<< _id << " is zero !" << endl;
   	return 0;
   }

   #ifdef _COMPUTE_MOMENTS_
   TDeme* curDeme=this;
   #endif
   TNode *Casc1, *Casc2;
   ChainNodeList.resetIterator();
   TNode * currentNode= (TNode*) ChainNodeList[0];

   double Tot_rec_proba=0;

   bool *canRecomb= new bool[_lineages];

   for   (int cpt=0; cpt < _lineages; ++cpt, currentNode=ChainNodeList.next())   {
      int min_pos = currentNode->min_pos;                                       // Remember that the position begin to 0
      int max_pos = currentNode->max_pos;                                       // and that max_pos is the Number
                                                                                // of linked loci minus 1
      //double recRate = myDemes->get_rec_rate();
      double rec_proba=currentNode->update_recombination_proba(min_pos,         // Decide if the recombination occurs
                                                               max_pos,
                                                               myDemes->get_rec_rate_perLoc());
      if(min_pos != max_pos){
         canRecomb[cpt]=true;
         Tot_rec_proba+=rec_proba;
      }
      else{
         canRecomb[cpt]=false;
      }
   }

   bool DoRecomb=true;
   long recombNode=0;
   if (Tot_rec_proba > 1) {
      recombNode= (int) ran3(&_lidum)*_lineages;
   }
   else {
      double isRecomb=ran3(&_lidum);
      if (isRecomb < Tot_rec_proba)  {
         recombNode=  (int) ran3(&_lidum)*_lineages;
         while(!canRecomb[recombNode]){
            recombNode=  (int) ran3(&_lidum)*_lineages;
         }
      }
      else{
         DoRecomb=false;
      }
   }

   if(DoRecomb)   {
      ChainNodeList.resetIterator();
      currentNode= (TNode*) ChainNodeList[0];

      for(int cpt=0; cpt < _lineages; ++cpt, currentNode=ChainNodeList.next())   {

         if(cpt == recombNode)   {
            int min_pos = currentNode->min_pos;                                       // Remember that the position begin to 0
            int max_pos = currentNode->max_pos;                                       // and that max_pos is the Number
                                                                                   // of linked loci minus 1
            //double recRate = myDemes->get_rec_rate();
            double rec_proba=currentNode->update_recombination_proba(   min_pos,         // Decide if the recombination occurs
                                                                  max_pos,
                                                                  myDemes->get_rec_rate_perLoc());
            //double temp =ran3(&_lidum);
            //if (temp < rec_proba){                                                    // Decide where the recombination occurs
            //int recomb_pos = currentNode->recombinationPosition(min_pos,max_pos);// Uniform recombination rate

            int recomb_pos = currentNode->recombinationPosition(  min_pos,           // different recombination rate between loci
                                                            max_pos,
                                                            myDemes->get_rec_rate_perLoc(),
                                                            rec_proba);

               #ifdef _COMPUTE_MOMENTS_
               myDemes->add_recombHits(curDeme,recomb_pos);                           //Compute between loci distribution of the recombination Hits
               if(( myDemes->GeneTree.MRCA_list[recomb_pos] && myDemes->GeneTree.MRCA_list[recomb_pos+1])
               || ( !currentNode->flag[recomb_pos] || !currentNode->flag[recomb_pos + 1]) ) {
                  myDemes->add_no_effective_recombHits(curDeme,recomb_pos);
               }
               #endif

            Casc1 = new TNode();
            Casc1->min_pos=currentNode->min_pos;                                   // create the first recombinating TNode
            for (int f=recomb_pos; f>=Casc1->min_pos; --f) {                       // Update the max_pos of asc1
               if (currentNode->flag[f]) {
                  Casc1->max_pos=f;
                  break;
               }
            }
            Casc1->time=time;
			
			//MJJ
		Casc1->deme=_id;

			
            Casc1->desc1=currentNode;
            Casc1->desc2=NULL;

            for (int cpt=0; cpt<=recomb_pos; ++cpt){                               // Update its flag array
               Casc1->flag[cpt]=currentNode->flag[cpt];
            }
            for (int cpt=recomb_pos+1; cpt<Casc1->num_linked_loci; ++cpt) {
               Casc1->flag[cpt]=false;
            }

            Casc2 = new TNode();                                                   // create the second recombinaing TNode
            Casc2->max_pos=currentNode->max_pos;                                   // Update the min_pos of asc2
            for (int f=recomb_pos+1; f<=Casc2->max_pos; ++f) {
               if (currentNode->flag[f]) {
                  Casc2->min_pos=f;
                  break;
               }
            }
            Casc2->time=time;

//MJJ
		Casc2->deme=_id;			

			
            Casc2->desc1=currentNode;
            Casc2->desc2=NULL;

            for (int cpt=0; cpt<=recomb_pos; cpt++){                               // Update its flag array
               Casc2->flag[cpt]=false;
            }
            for (int cpt=recomb_pos+1; cpt<Casc2->num_linked_loci; ++cpt) {
               Casc2->flag[cpt]=currentNode->flag[cpt];
            }

            currentNode->ancestor=Casc1;
            currentNode->asc1=Casc1;
            currentNode->asc2=Casc2;
            // add asc1 and asc2
            //in the Chained List chainTree

            GeneTree.chainTree.Add(Casc1);
            Casc1->ID_Node=(GeneTree.chainTree.get_size() - 1 );
            Casc1->event = 2;

            GeneTree.chainTree.Add(Casc2);
            Casc2->ID_Node=(GeneTree.chainTree.get_size() - 1 );
            Casc2->event = 2;

            ChainNodeList.replace_elem(Casc1);
            ChainNodeList.Add(Casc2);

            ++_lineages;
            break;
         }
      }
   }
   if(canRecomb)   delete[]    canRecomb;
   return 1;
}
//----------------------------------------------------------------------------
//Update the recombination rate for a given TNode
double
TNode::update_recombination_proba(const int& min,
                                  const int& max,
                                  const TRecRateArray* RecRatePerLoc) {
   double rec_proba =0;
   for (int i=min; i<max; ++i) {
      rec_proba += *RecRatePerLoc->elem(i);
   }
   return rec_proba;
}
//----------------------------------------------------------------------------
//Put a recombination hit in a TNode: Only for uniform recombination rate
inline
int
TNode::recombinationPosition(const int & min, const int & max){
    int Pos_break= (int) ran3(&_lidum)*(max-min); // integer random draw between 0 and (max-min-1)
    return (min + Pos_break) ;

}
//----------------------------------------------------------------------------
//Put a recombination hit in a TNode: different recombination rates between linkage blocks
int
TNode::recombinationPosition(const int& min,
                             const int& max,
                             const TRecRateArray* RecRatePerLoc,
                                   double recProba){

    int Pos_break=-1;
    double random= ran3(&_lidum)*(1); // float random draw between 0 and 1
    double  Proba_rec=0;
    for(int i=min; i<max;++i) {
       Proba_rec+=*RecRatePerLoc->elem(i)/recProba;
       if( Proba_rec >= random ) {
         Pos_break=i;
         break;
       }
    }
    return (Pos_break);
}
//----------------------------------------------------------------------------
//A destructor of trees with recombination
void
TNode::destroy(){
    int temp=-1;
    if(asc1) temp=this->asc1->ID_Node;
    //cout << (this->ID_Node +1) << " and " << (temp +1) << endl;

    if (desc1) {
    	desc1->destroy();
    	if (desc1->asc2) desc1->asc2->desc1=NULL;
    	delete desc1;
    	desc1=NULL;
    }

    if (desc2) {
    	desc2->destroy();
    	if (desc2->asc2) desc2->asc2->desc1=NULL;
    	delete desc2;
    	desc2=NULL;
    }
}
//----------------------------------------------------------------------------
//Allocate the variable recombHits
void
TDemeCollection::allocate_recombHits(){
    recombHits = new int[_num_demes*num_indep_loci*(num_linked_loci-1)];
    effective_recombHits = new int[_num_demes*num_indep_loci*(num_linked_loci-1)];
}

//----------------------------------------------------------------------------
//Allocate the variable recombHits
void
TDemeCollection::allocate_recombHits(const int num_tot_interval){
    recombHits = new int[_num_demes*num_tot_interval];
    effective_recombHits = new int[_num_demes*num_tot_interval];
}
//----------------------------------------------------------------------------
//Allocate the variable recombHits
void
TDemeCollection::initialize_recombHits(){
    for(int cpt=0;cpt<_num_demes*num_indep_loci*(num_linked_loci-1);cpt++) {
      recombHits[cpt]=0;
      effective_recombHits[cpt]=0;
    }
}
//----------------------------------------------------------------------------
//increment when a recombination occurs
inline
void
TDemeCollection::add_recombHits(const TDeme* curDeme, const int recomb_pos){
   //int winBlock;
   //cout << endl << ind_loci <<endl;
   //cout << recomb_pos <<endl;
   int id_pop=curDeme->get_id();
   recombHits[id_pop*(num_indep_loci)*(num_linked_loci-1)+(ind_loci*(num_linked_loci-1)+recomb_pos)]+=1;
   //cin >> winBlock;
}
//----------------------------------------------------------------------------
//increment when a recombination occurs
inline
void
TDemeCollection::add_no_effective_recombHits(const TDeme* curDeme, const int recomb_pos){
   //int winBlock;
   //cout << endl << ind_loci <<endl;
   //cout << recomb_pos <<endl;
   int id_pop=curDeme->get_id();
   effective_recombHits[id_pop*(num_indep_loci)*(num_linked_loci-1)+(ind_loci*(num_linked_loci-1)+recomb_pos)]+=1;
   //cin >> winBlock;
}
//----------------------------------------------------------------------------
//return number of hits due to recombination between two ajacent loci
int
TDemeCollection::return_recombHits(const int linked_loci){
   return recombHits[linked_loci];
}
//----------------------------------------------------------------------------
//return number of hits due to recombination between two ajacent loci
int
TDemeCollection::return_no_effective_recombHits(const int linked_loci){
   return effective_recombHits[linked_loci];
}
//----------------------------------------------------------------------------
//return number of hits due to recombination between two ajacent loci
int
TDemeCollection::number_recombHits_perDeme(const TDeme* curDeme){
   int sum=0;
   int id_pop=curDeme->get_id();
   for(int cpt=0;cpt<num_indep_loci*(num_linked_loci-1);cpt++) {
      sum+=recombHits[id_pop*num_indep_loci*(num_linked_loci-1) + cpt];
   }
   return sum;
}
//----------------------------------------------------------------------------
//return number of hits due to recombination between two ajacent loci
int
TDemeCollection::number_no_effective_recombHits_perDeme(const TDeme* curDeme){
   int sum=0;
   int id_pop=curDeme->get_id();
   for(int cpt=0;cpt<num_indep_loci*(num_linked_loci-1);cpt++) {
      sum+=effective_recombHits[id_pop*num_indep_loci*(num_linked_loci-1) + cpt];
   }
   return sum;
}
//----------------------------------------------------------------------------
//return number of hits due to recombination between two ajacent loci
int
TDemeCollection::number_recombHits_perDeme_perIndepLoci(const TDeme* curDeme, const int & indep_loci){
   int sum=0;
   int id_pop=curDeme->get_id();
   for(int cpt=0;cpt<(num_linked_loci-1);cpt++) {
      sum+=recombHits[id_pop*num_indep_loci*(num_linked_loci-1) + (indep_loci*(num_linked_loci-1)+cpt)];
   }
   return sum;
}
//----------------------------------------------------------------------------
//return number of hits due to recombination between two ajacent loci
int
TDemeCollection::number_no_effective_recombHits_perDeme_perIndepLoci(const TDeme* curDeme, const int & indep_loci){
   int sum=0;
   int id_pop=curDeme->get_id();
   for(int cpt=0;cpt<(num_linked_loci-1);cpt++) {
      sum+=effective_recombHits[id_pop*num_indep_loci*(num_linked_loci-1) + (indep_loci*(num_linked_loci-1)+cpt)];
   }
   return sum;
}
//----------------------------------------------------------------------------
//return a reference to the Deme i
const TDeme&
TDemeCollection::get_deme(const int & i){
   return Demes[i];
}
//----------------------------------------------------------------------------
void TDeme::resetAllNodeFlags(const int& numLinkedLoci) {
   ChainNodeList.resetIterator();
   TNode * curNode = (TNode*) ChainNodeList[0];
   for (int j=0; j<_lineages; ++j, curNode=ChainNodeList.next()) {
      curNode->num_linked_loci=numLinkedLoci;
      curNode->min_pos=0;
      curNode->max_pos=numLinkedLoci-1;
      if (curNode->flag) delete [] curNode->flag;
      curNode->flag = new bool[curNode->num_linked_loci];
      for (int k=0; k<curNode->num_linked_loci; ++k) {
         curNode->flag[k]=true;
      }
   }
}
//----------------------------------------------------------------------------
// Need to resize all vectores for a potentially new sie if different chromosome structure
void TDemeCollection::set_num_linked_loci(const int& numLinkedLoci) {
   num_linked_loci=numLinkedLoci;
   for (int i=0; i<num_demes(); ++i) {
      TDeme & curDeme = Demes[i];
      //Readjust size of flag vector
      curDeme.resetAllNodeFlags(numLinkedLoci);
      //Readjust size of sequence if needed
      curDeme.resetSequenceSize();
   }
}
//----------------------------------------------------------------------------
void TDeme::resetSequenceSize() {
   ChainNodeList.resetIterator();
   TNode * curNode = (TNode*) ChainNodeList[0];
   for (int j=0; j<_lineages; ++j, curNode=ChainNodeList.next()) {
      if (curNode->sequence) delete curNode->sequence;
      curNode->sequence = new TIntVect(curNode->num_linked_loci, 0);
   }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

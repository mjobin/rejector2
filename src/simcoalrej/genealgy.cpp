#include "public.h"
#include "genealgy.h"
#include "deme.h"

#include <stdexcept>
#include <exception>

class Exception;

enum UNITS {WIDTH=6};



//Initialization of static members of TNode
TIntVect 	TNode::hits(0,10);
GammaRates 	TNode::mut_rates(0,0.0,1.0,0.0, 0);
double 		TNode::mut_rate=0.0;
double 		TNode::tree_length=0.0;
int 			TNode::node_count=0;
int 			TNode::min_mic=0;
int 			TNode::max_mic=0;
/*****************************************************************************

                       For recombination
            Initialization of the new static members of TNode

******************************************************************************/
//double 		TNode::rec_rate=0.0;
int 			TNode::num_linked_loci=0;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
TNode&
TNode::operator=(const TNode& N) {
	time=N.time;
   desc1=N.desc1;
   desc2=N.desc2;
   ancestor=N.ancestor;
   num_desc1=N.num_desc1;
   num_desc2=N.num_desc2;
   node_number=N.node_number;
   deme=N.deme;
   num_new_mut=N.num_new_mut;
   left_mut=N.left_mut;
   right_mut=N.right_mut;
   seq_length=N.seq_length;
   
   
   // Affectation of the new properties needed for the recombination
   
   asc1=N.asc1;
   asc2=N.asc2;
   ID_Node=N.ID_Node;
   min_pos=N.min_pos;
   max_pos=N.max_pos;
   event=N.event;
   flag=N.flag;
   _lidum=N._lidum;
   
   
   
   if (desc1_nodes) 	delete desc1_nodes; 	desc1_nodes=NULL;
   if (desc2_nodes) 	delete desc2_nodes; 	desc2_nodes=NULL;
   if (sequence) 		delete sequence; 		sequence=NULL;
   if (N.desc1_nodes) {
   	try {
      	int size=N.desc1_nodes->GetItemsInContainer();
      	desc1_nodes= new TIntVect(size,10);
         for (int i=0; i<size; ++i) {
         	desc1_nodes[i]=N.desc1_nodes[i];
         }
      }
      catch (...) {
   		if (desc1_nodes) 	delete desc1_nodes; 	desc1_nodes=NULL;
         cout << "\nTNode::operator=(): unable to allocate memory\n";
         return *this;
      }
   }
   
   if (N.desc2_nodes) {
   	try {
      	int size=N.desc2_nodes->GetItemsInContainer();
      	desc2_nodes= new TIntVect(size,10);
         for (int i=0; i<size; ++i) {
         	desc2_nodes[i]=N.desc2_nodes[i];
         }
      }
      catch (...) {
   		if (desc2_nodes) 	delete desc2_nodes; 	desc2_nodes=NULL;
         cout << "\nTNode::operator=(): unable to allocate memory\n";
         return *this;
      }
   }

   if (N.sequence) {
   	try {
      	int size=N.sequence->GetItemsInContainer();
      	sequence= new TIntVect(size,10);
         for (int i=0; i<size; ++i) {
         	sequence[i]=N.sequence[i];
         }
      }
      catch (...) {
   		if (sequence) 	delete sequence; 	sequence=NULL;
         cout << "\nTNode::operator=(): unable to allocate memory\n";
         return *this;
      }
   }
   return *this;
}
//----------------------------------------------------------------------------
// A recursive method to count the number of nodes below a given node,
// including himself
int
TNode::count_desc() {
	num_desc1=num_desc2=0;
	if (desc1) {
      num_desc1=desc1->count_desc();
   }
	if (desc2) {
      num_desc2=desc2->count_desc();
   }
   //Counts only the basal tips
   if (!desc1 && !desc2) { 
      return num_desc1+num_desc2+1;
   }
   return num_desc1+num_desc2;
}
//----------------------------------------------------------------------------
// A recursive method to count the number of nodes below a given node,
// including himself
// modification for recombination: compute for each linked loci and exclude recombinating Nodes
int
TNode::count_desc(const int& loc) {
	num_desc1=num_desc2=0;
   bool wrongNode=false;
   //TNode* CurrNode;
   //CurrNode=this;
   TNode* tempNode;
	if (desc1) {
      tempNode=desc1;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[loc] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc] && tempNode->desc2->flag[loc] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }   
      }
      num_desc1=tempNode->count_desc(loc);  
   }
	if (desc2) {
      tempNode=desc2;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[loc] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc] && tempNode->desc2->flag[loc] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }    
      }
      num_desc2=tempNode->count_desc(loc);  
   }
   //Counts only the basal tips
   if (!desc1 && !desc2) {
      return num_desc1+num_desc2+1;
   }
   return num_desc1+num_desc2;
}
//----------------------------------------------------------------------------
// A recursive method to count the number of nodes below a SNP mutation
int
TNode::count_SNP_desc(const int& loc) {
	//num_desc1=num_desc2=0;
   int num_SNP_desc=0;
   bool wrongNode=false;
   TNode* tempNode;
	if (desc1) {
      tempNode=desc1;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[loc] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc] && tempNode->desc2->flag[loc] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }   
      }
      //num_desc1=tempNode->count_desc(loc);
      num_SNP_desc+=tempNode->count_desc(loc);  
   }
	if (desc2) {
      tempNode=desc2;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[loc] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc] && tempNode->desc2->flag[loc] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[loc]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }    
      }
      //num_desc2=tempNode->count_desc(loc);
      num_SNP_desc+=tempNode->count_desc(loc);  
   }
   //Counts only the basal tips
   if (!desc1 && !desc2) {
      //return num_desc1+num_desc2+1;
      return num_SNP_desc+1;
   } 
   //return num_desc1+num_desc2;
   return num_SNP_desc;
}

//----------------------------------------------------------------------------
int
TNode::count_polym_sites() {
	int polym_loci=0, size=hits.GetItemsInContainer();
	for (int i=0; i<size; ++i) {
   	if (hits[i]) ++polym_loci;
   }
   return polym_loci;
}
//----------------------------------------------------------------------------
//A recursion procedure to establish the list of descendent nodes to the left
//and to the right of the current node
int
TNode::build_lists_of_descendent_nodes() {

	//Creates empty list of descendent nodes

		try {
   		if (!desc1_nodes) desc1_nodes= new TIntVect(1,10);
         else desc1_nodes->Flush();
      	if (!desc2_nodes) desc2_nodes= new TIntVect(1,10);
         else desc2_nodes->Flush();
   	}
   	catch (...) {
   		if (desc1_nodes) delete desc1_nodes; desc1_nodes=NULL;
      	if (desc2_nodes) delete desc2_nodes; desc2_nodes=NULL;
      	cout 	<< "TNode:: build_list_of_descendent_nodes() : Unable to allocate memory"
      			<< endl;
      	return 0;
   	}

   //Begin recursion
   if (desc1) {
   	//Adds the descendent node list of descent nodes to the current one
   	if (!desc1->build_lists_of_descendent_nodes()) return 0;
      int num_nodes=desc1->desc1_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc1_nodes->Add((*desc1->desc1_nodes)[i]);
      }
      num_nodes=desc1->desc2_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc1_nodes->Add((*desc1->desc2_nodes)[i]);
      }
   }
   else {
   	//Adds its id in its own list of descendent nodes !
   	desc1_nodes->Add(node_number);
   }
   if (desc2) {
   	//Adds the descendent node list of descent nodes to the current one
   	if (!desc2->build_lists_of_descendent_nodes()) return 0;
      int num_nodes=desc2->desc1_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc2_nodes->Add((*desc2->desc1_nodes)[i]);
      }
      num_nodes=desc2->desc2_nodes->GetItemsInContainer();
      for (int i=0; i<num_nodes; ++i) {
      	desc2_nodes->Add((*desc2->desc2_nodes)[i]);
      }
   }
   return 1;
}
//----------------------------------------------------------------------------
//Also a recursion method, assuming that method count_desc()
//has been performed before
int
TNode::compute_moments_of_pairwise_divergence(double& sum_time, double& sum_square_time) {
  if (desc1 && desc2) {
      //Step 1: The sum of divergence times is simply obtained from the knowledge
      //of the number of descendents
      sum_time+=2.0*time*num_desc1*num_desc2;
      sum_square_time+=4.0*time*time*num_desc1*num_desc2;
      //!!!Beware: This is certainly wrong. No analogy with time
      //Step 2: Now simply recurse through the tree
      desc1->compute_moments_of_pairwise_divergence(sum_time,sum_square_time);
      desc2->compute_moments_of_pairwise_divergence(sum_time,sum_square_time);
   }
   return 1;
}

//----------------------------------------------------------------------------
//Also a recursion method, assuming that method count_desc()
//has been performed before
//modified for recombination: computation for every linked loci
int                                                              
TNode::compute_moments_of_pairwise_divergence(double& sum_time, 
                                              double& sum_square_time, const int& loc) {

  //if (desc1 && desc2) {
  //Step 1: The sum of divergence times is simply obtained from the knowledge
  //of the number of descendents
  sum_time+=2.0*time*num_desc1*num_desc2;
  sum_square_time+=4.0*time*time*num_desc1*num_desc2;
  //!!!Beware: This is certainly wrong. No analogy with time
  //Step 2: Now simply recurse through the tree 
  //after testing if Nodes are true coalescing Nodes
  bool wrongNode=false;
  TNode* CurrNode;
  //CurrNode=this;
  TNode* tempNode;
  if (desc1) {
     tempNode=desc1;
     wrongNode=true;
     while(wrongNode) {
        if(tempNode->flag[loc] ) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[loc] && tempNode->desc2->flag[loc] ) { 
                 wrongNode=false;
              }
           }
           else {
              if (!tempNode->desc1 && !tempNode->desc2) {
                 wrongNode=false;
              }
           }         
        }
        if(wrongNode) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[loc]) {
                 tempNode=tempNode->desc1;
              }
              else {
                 tempNode=tempNode->desc2;
              }
           }
           else {
              if(!tempNode->desc1) {
                    tempNode=tempNode->desc2;
              }
              else {
                 tempNode=tempNode->desc1;
              }
           }
        }   
     }
     tempNode->compute_moments_of_pairwise_divergence(sum_time,sum_square_time,loc);
  }
  if (desc2) {
     tempNode=desc2;
     wrongNode=true;
     while(wrongNode) {
        if(tempNode->flag[loc] ) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[loc] && tempNode->desc2->flag[loc] ) { 
                 wrongNode=false;
              }
           }
           else {
              if (!tempNode->desc1 && !tempNode->desc2) {
                 wrongNode=false;
              }
           }         
        }
        if(wrongNode) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[loc]) {
                 tempNode=tempNode->desc1;
              }
              else {
                 tempNode=tempNode->desc2;
              }
           }
           else {
              if(!tempNode->desc1) {
                 tempNode=tempNode->desc2;
              }
              else {
                 tempNode=tempNode->desc1;
              }
           }
        }    
     }
     tempNode->compute_moments_of_pairwise_divergence(sum_time,sum_square_time,loc);
  }
  //}
   return 1;
}
//----------------------------------------------------------------------------
//WITHOUT RECOMBINATION
//Guillaume 05 _04_2004
//Compute total size of Tree, size of Tree for every TNode and List of number of
//gene in the sample for every TNode
int
TNode::compute_TreeSize_SizeToNode_NumDescList(int*   i,
                                               int&   sum_time,
                                               int*   sizeToNode,
                                               //int*   numDescList,
                                               TNode** scannedNodes) {
  if (desc1 && desc2) {
      int j=*i+1;
      *i=j;
      sum_time+=time-desc1->time;
      //sizeToNode[*i]=time-desc1->time;
      sizeToNode[*i]=sum_time;
      //numDescList[*i]=desc1->count_SNP_desc(0);
      scannedNodes[*i]=desc1;
      desc1->compute_TreeSize_SizeToNode_NumDescList(   i,
                                                        sum_time,
                                                        sizeToNode,
                                                        //numDescList,
                                                        scannedNodes);

      j=*i+1;
      *i=j;
      sum_time+=time-desc2->time;
      //sizeToNode[*i]=time-desc2->time;
      sizeToNode[*i]=sum_time;
      //numDescList[*i]=desc2->count_SNP_desc(0);
      scannedNodes[*i]=desc2;
      desc2->compute_TreeSize_SizeToNode_NumDescList(   i,
                                                        sum_time,
                                                        sizeToNode,
                                                        //numDescList,
                                                        scannedNodes);
   }
   return 1;
}
//----------------------------------------------------------------------------
//Function to compute the size of the coalescent tree (sum of of all the branches)
//for the current locus
//modified for recombination: computation for every linked loci
int
TNode::compute_TreeSize(double& sum_time, const int& curLocus) {

  bool wrongNode=false;
  TNode* tempNode;
  if (desc1) {
     tempNode=desc1;
     wrongNode=true;
     while(wrongNode) {
        if(tempNode->flag[curLocus] ) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) {
                 wrongNode=false;
              }
           }
           else {
              if (!tempNode->desc1 && !tempNode->desc2) {
                 wrongNode=false;
              }
           }
        }
        if(wrongNode) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[curLocus]) {
                 tempNode=tempNode->desc1;
              }
              else {
                 tempNode=tempNode->desc2;
              }
           }
           else {
              if(!tempNode->desc1) {
                    tempNode=tempNode->desc2;
              }
              else {
                 tempNode=tempNode->desc1;
              }
           }
        }   
     }
     sum_time+=time-tempNode->time;
     tempNode->compute_TreeSize(sum_time,curLocus);
  }
  if (desc2) {
     tempNode=desc2;
     wrongNode=true;
     while(wrongNode) {
        if(tempNode->flag[curLocus] ) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) {
                 wrongNode=false;
              }
           }
           else {
              if (!tempNode->desc1 && !tempNode->desc2) {
                 wrongNode=false;
              }
           }
        }
        if(wrongNode) {
           if(tempNode->desc1 && tempNode->desc2) {
              if(tempNode->desc1->flag[curLocus]) {
                 tempNode=tempNode->desc1;
              }
              else {
                 tempNode=tempNode->desc2;
              }
           }
           else {
              if(!tempNode->desc1) {
                 tempNode=tempNode->desc2;
              }
              else {
                 tempNode=tempNode->desc1;
              }
           }
        }    
     }
     sum_time+=time-tempNode->time;
     tempNode->compute_TreeSize(sum_time,curLocus);
  }
  return 1;
}

//----------------------------------------------------------------------------
//A recursion method to compute the total coalescent times,
//both within and among demes
//Could be easily modified to get the distribution of coalescent times
int
TNode::compute_total_coal_times_among_demes( TDemeCollection& DemeCollec,
   													  TMigrationMatrix& CoalTimes,
                                            TMigrationMatrix& PairDiff,
   													  TMigrationMatrix& MinCoalTimes) {
   if (desc1 && desc2) {
   	//Step 1: Count the number of descendent nodes
   	int 	left_desc=desc1_nodes->GetItemsInContainer(),
      		right_desc=desc2_nodes->GetItemsInContainer();

      int deme1, deme2;
      //Step 2: Explore the list of descendent nodes to update coalescence times
      for (int i=0; i<left_desc; ++i) {
      	deme1=DemeCollec.tree_node((*desc1_nodes)[i]-1).deme;
         //Count the descendent nodes belonging to the same deme as deme1
         for (int j=0; j<right_desc; ++j) {
      		deme2=DemeCollec.tree_node((*desc2_nodes)[j]-1).deme;
            if (deme2<deme1)	{  //Between pop treatment
            	CoalTimes(deme1,deme2)+=time;
            	if (time<MinCoalTimes(deme1,deme2)) MinCoalTimes(deme1,deme2)=time;
            }
            else {
            	CoalTimes(deme2,deme1)+=time;
               if (time<MinCoalTimes(deme2,deme1)) MinCoalTimes(deme2,deme1)=time;
            }
         }
      }

      //Step 3: Now count the number of pairwise differences
      int 	size=DemeCollec.num_demes();
      int 	*deme_count_left,*deme_count_right;                                                      
      
      try {
      	deme_count_left=new int[size],
      	deme_count_right=new int[size];
         for (int i=0; i<size; ++i) {
         	deme_count_left[i]=0;
         	deme_count_right[i]=0;
         }
      }
      catch (...) {
      	cout << "TNode::compute_total_coal_times_among_demes: unable to allocate memory\n";
         if (deme_count_left) delete[] deme_count_left;
         if (deme_count_right) delete[] deme_count_right;
         return 0;
      }
      //Count the number of nodes of each deme below this node
      for (int i=0; i<left_desc; ++i)
      	++deme_count_left[DemeCollec.tree_node((*desc1_nodes)[i]-1).deme];
      
      for (int i=0; i<right_desc; ++i)
      	++deme_count_right[DemeCollec.tree_node((*desc2_nodes)[i]-1).deme];

      int descnodesi, descnodesj, samp_sizei;
      for (int i=0; i<size; ++i) { 
      	descnodesi=deme_count_left[i]+deme_count_right[i];
         samp_sizei=DemeCollec[i].sample_size();
      	//Between pop treatment
      	for (int j=0; j<i; ++j) {
         	descnodesj=deme_count_left[j]+deme_count_right[j];
            PairDiff(i,j)+=num_new_mut*
            					( descnodesi*(DemeCollec[j].sample_size()-descnodesj)+
                            (descnodesj)*(samp_sizei-descnodesi) );
         }
         //Within pop treatment
         PairDiff(i,i)+=num_new_mut*descnodesi*(samp_sizei-descnodesi); //OK
      }
      //Recurse through the tree
      desc1->compute_total_coal_times_among_demes(DemeCollec,CoalTimes, PairDiff, MinCoalTimes);
      desc2->compute_total_coal_times_among_demes(DemeCollec,CoalTimes, PairDiff, MinCoalTimes);

      if (deme_count_left) delete[] deme_count_left;
      if (deme_count_right) delete[] deme_count_right;

   }
   else {
   	PairDiff(deme,deme)+=num_new_mut*(DemeCollec[deme].sample_size()-1);
      for (int i=0; i<DemeCollec.num_demes(); ++i) {
      	if (deme>i)
         	PairDiff(deme,i)+=num_new_mut*DemeCollec[i].sample_size();
         else if (i>deme)
         	PairDiff(i,deme)+=num_new_mut*DemeCollec[i].sample_size();
      }
   }
	return 1;
}
//----------------------------------------------------------------------------
ostream & operator<<(ostream& os, const TNode& node) {
	if (node.time>0) {
		os << "\nInternal node No " << node.node_number
      	<< "\n----------------"
      	<< "\nTime since present       : " << node.time
   		<< "\nNo. of left descendents  : " << node.num_desc1
   		<< "\nNo. of right descendents : " << node.num_desc2
/*      	<< "\nLeft descendent address  : " << node.desc1
      	<< "\nRight descendent address : " << node.desc2
      	<< "\nAncestor address         : " << node.ancestor */<< endl;
	}
   else os << "\nTip node No " << node.node_number
           << "\n-----------" << endl;
   return os;
}
//----------------------------------------------------------------------------
// A recursive method to print node information recursively
int
TNode::print_info(ostream& os) {
	os << *this;
	if (desc1) desc1->print_info(os);
	if (desc2) desc2->print_info(os);
   return 1;
}
/*
void
TNode::print_desc_nodes(TDrawingBoard& DB, int node_posx, int node_posy,
      								char node, char hor_bar, char left_corner,
                              char right_corner, char vert_bar) {
   //Draw Node number
   char num[20];
   itoa(this->node_number, num, 10);
   DB.draw_text(node_posx, node_posy+1,num);
   if (desc1 && desc2) {
   	int l11, l12, r11, r12;
      l11=node_posx-(int)((float)(num_desc1+num_desc2-1)*WIDTH)/2;
      l12=l11+(num_desc1-1)*WIDTH;
      r12=l11+(num_desc1+num_desc2-1)*WIDTH;
      r11=r12-(num_desc2-1)*WIDTH;
   	int pos_vert_line1=(int)((float)(l11+l12)/2),
          pos_vert_line2=(int)((float)(r11+r12)/2),
          vert_line1_height=(int)((this->time-desc1->time)*UNIT_TIME)-2,
          vert_line2_height=(int)((this->time-desc2->time)*UNIT_TIME)-2;
   	//Draw left corner
   	DB.draw_char(pos_vert_line1, node_posy, left_corner);
		//Draw first line segment                         mod
      DB.draw_hor_line(pos_vert_line1+1,node_posy, node_posx-pos_vert_line1-1, hor_bar);
   	//Draw node
   	DB.draw_char(node_posx, node_posy, node);
		//Draw second line
   	DB.draw_hor_line(node_posx+1,node_posy, pos_vert_line2-node_posx-1, hor_bar);
   	//Draw right corner
      DB.draw_char(pos_vert_line2, node_posy, right_corner);
      //Draw vertical lines until next nodes
      DB.draw_vert_line(pos_vert_line1,node_posy+1, vert_line1_height, vert_bar);
      DB.draw_vert_line(pos_vert_line2,node_posy+1, vert_line2_height, vert_bar);
		//Call descendents methods
      desc1->print_desc_nodes(DB,pos_vert_line1,node_posy+vert_line1_height+2,
      								node,hor_bar,left_corner, right_corner, vert_bar);
      desc2->print_desc_nodes(DB,pos_vert_line2,node_posy+vert_line2_height+2,
      								node,hor_bar,left_corner, right_corner, vert_bar);
   }
}
*/

//----------------------------------------------------------------------------
//Version without recombination
//A recursive procedure to add mutations to the tree starting from the root
long
TNode::add_mutations(long              *lidum, 
                     const int&        num_mut,
      		         const int&        len, 
                     const Mut_Type&   mut_type, 
                     const double&     gamma_par,
                     const double&     trans_rate,
                     const int&        range_const) {

   long tot_mut=0L, desc_mut;

   float subst_length;

   seq_length=len;

   //Create new sequence and propagate the ancestor's mutations
   if (!sequence)  {
   	try {
    		sequence= new TIntVect(len, 0);
   	}
   	catch (...) {
   		if (sequence) delete sequence; sequence=NULL;
      	seq_length=0;
      	return 0;
   	}
   }
   if (ancestor) {
   	int *cursite=sequence->begin(), *cursite_anc=ancestor->sequence->begin();
   	for (int i=0; i<len; ++i) {
      	*cursite++=*cursite_anc++; //Inheritance of ancestor's mutations
      	//(*sequence)[i]=(*ancestor->sequence)[i]; //Inheritance of ancestor's mutations
      }
   }
   else  {    //Loro_15_9_98 Create the ancestor's genetic data
   	int* cursite=sequence->begin();
   	if (mut_type==DNA) {
      	 //for (int i=0; i<len; ++i) (*sequence)[i]=ran3(lidum)*4; //A random number between 0 and 3
			for (int i=0; i<len; ++i) *cursite++= (int) (ran3(lidum)*4); //A random number between 0 and 3
        }
        //else for (int i=0; i<len; ++i) (*sequence)[i]=0;
        else for (int i=0; i<len; ++i) *cursite++=0;

   }

   //This is the number of new mutations to generate as compared to parent node
   num_new_mut=num_mut;

   double drand_num;
   //Generate those mutations
   for (int i=0, pos; i<num_mut; ++i) {
      //Find the position of the site to be hit by a mutation
   	if (fabs(gamma_par)<1e-7) { //Close to zero, even mutation rates are assumed among sites
      	pos= (int) (ran3(lidum)*len);
      }
      else {
      	drand_num=ran3(lidum);
         double tot_prob=0.0;
         /*
         for (pos=-1; pos<(len-1); ) {
            ++pos;
         	tot_prob+=mut_rates[pos];
            if (drand_num<tot_prob) break;
         }
         */
         //Explore from the end of the sequence because it is where mutatins are
         //more likely to occur
         for (pos=len-1; pos>-1; --pos) {
         	tot_prob+=mut_rates[pos];
            if (drand_num<tot_prob) break;
         }
      }
      ++hits[pos];
   	int& site=(*sequence)[pos];
   	switch (mut_type) {
      	case MICROSAT  :	//Loro_04_03_04
                        if (range_const) {
                           	//There are bouncing walls at min_mic and max_mic
                                if (site==min_mic) ++site;
                                else if (site==max_mic) --site;
                                //otherwise they are free to move randomly
                                else if (ran3(lidum)<0.5) --site;
                                else ++site;
                        }
                        else {
                                if (ran3(lidum)<0.5) --site;
                                else ++site;
                        }
         				 		break;
         case DNA			:  drand_num=ran3(lidum);
         						//Here we implement a 95% transition bias
                           //0 and 1 : A and G
                           //2 and 3 : C and T 
                           switch (site) {
                           case 0: 	if (drand_num<trans_rate) site=1;
                           			else
                                    	if (ran3(lidum)<0.5) site=2;
                                       else site=3;
                                    break;
                           case 1: if (drand_num<trans_rate) site=0;
                           			else
                                    	if (ran3(lidum)<0.5) site=2;
                                       else site=3;
                                    break;
                           case 2: if (drand_num<trans_rate) site=3;
                           			else
                                    	if (ran3(lidum)<0.5) site=0;
                                       else site=1;
                                    break;
                           case 3: if (drand_num<trans_rate) site=2;
                           			else
                                    	if (ran3(lidum)<0.5) site=0;
                                       else site=1;
                                    break;
                           }
         						break;
         case RFLP		: 	site=!site;
         						break;
         default /*DNA*/: break;
      }
   }

   //Continue recursion
   if (desc1) {
   	subst_length=mut_rate*(time-desc1->time);
      desc_mut= (int) poidev(subst_length, lidum);
      tree_length+=subst_length;
      tot_mut+=desc_mut+desc1->add_mutations(lidum, desc_mut, len, mut_type,
      													gamma_par,trans_rate, range_const);
   }
   if (desc2) {
   	subst_length=mut_rate*(time-desc2->time);
      desc_mut= (int) poidev(subst_length, lidum);
      tree_length+=subst_length;
      tot_mut+=desc_mut+desc2->add_mutations(lidum, desc_mut, len, mut_type,
      													gamma_par,trans_rate, range_const);
   }
   return tot_mut;
}
/*
//----------------------------------------------------------------------------
//Function to add mutation locus per locus
//for every loci, this function transforme the whole tree into a coalescent tree 
//Inplemented for recombination (DON T USED)
long
TNode::add_mutations(   long *            lidum,
                        const int&        num_mut,
                        const int&        curLocus,
                        const int&        len,
      		            const int&        numLinkedLoci,
                        const Mut_Type&   mut_type,
                        const double&     mutRatePerLoc,
                        const double&     geomParamPerLoc,
                        const double&     gamma_par,
                        const double&     trans_rate,
                        const int&        range_const,
                        const TNode*      const pN) //Pointer to the calling node (ancestor without recombination)
                        {


   long tot_mut=0L, desc_mut, cur_mut;
   double subst_length;

   seq_length=numLinkedLoci;
   
   //bool firstTime=false;

   //Create new sequence and propagate the ancestor's mutations
   if (!sequence)  {
   	try {
    		sequence= new TIntVect(numLinkedLoci, 0);
         int *cursite=sequence->begin();
         for (int i=0; i<numLinkedLoci; ++i, ++cursite){
            *cursite= 999; //Initialise sequence to 999
         }
   	}
   	catch (...) {
   		if (sequence) delete sequence; sequence=NULL;
      	seq_length=0;
      	return 0;
   	}
   }
   
   if (pN) {
      int *cursite=sequence->elem(curLocus), *cursite_anc=pN->sequence->elem(curLocus);
      
      for (int i=curLocus; i<curLocus+len; ++i, ++cursite, ++cursite_anc) {
         if (pN->flag[i]) { //To be sure that we only copy the site of interest
            *cursite=*cursite_anc; //Inheritance of ancestor's mutations
            //cout << ID_Node  << ", " << loc << ": " <<*cursite;
         }
         //cout << endl;
      }
      
   }
   else  {    //Loro_15_9_98 Create the ancestor's genetic data
      int* cursite=sequence->elem(curLocus);
      if (mut_type==DNA) {
         for (int i=curLocus; i<curLocus+len; ++i, ++cursite){
            *cursite= (int) (ran3(lidum)*4); //A random number between 0 and 3
            //cout << ID_Node  << ", " << loc << ": " <<*cursite;
         }
         //cout << endl;
      }
      else for (int i=curLocus; i<curLocus+len; ++i) *cursite++=0;
   }

   //This is the number of new mutations to generate as compared to parent node
   num_new_mut=num_mut;

   cur_mut=0L;
   double drand_num;
   //Generate those mutations
   for (int i=0, pos; i<num_mut; ++i) {
      //Find the position of the site to be hit by a mutation
   	if (fabs(gamma_par)<1e-7) { //Close to zero, even mutation rates are assumed among sites
      	pos= (int) (ran3(lidum)*len);
      }
      else {
      	drand_num=ran3(lidum);
         double tot_prob=0.0;
         //Explore from the end of the sequence because it is where mutations are
         //more likely to occur
         for (pos=len-1; pos>-1; --pos) {
         	tot_prob+=mut_rates[pos];
            if (drand_num<tot_prob) break;
         }
      }
      
      if (pN->flag[pos]) { //We put the mutation only if it is transmitted to the current nodes
      	++cur_mut;
      	++hits[pos];
   		int& site=(*sequence)[pos];
   		switch (mut_type) {
      		case MICROSAT  :	//Loro_04_03_04
                              if (geomParamPerLoc>0) {
                                 int step=1+geometric(geomParamPerLoc);
                                 if (range_const)  {
                              	   //There are bouncing walls at min_mic and max_mic
                              	   if (site==min_mic) site+=step;
                                    else if (site==max_mic) site-=step;
                                    //otherwise they are free to move randomly
                                    else if (ran3(lidum)<0.5) site-=step;
      									   else site+=step;
                                    if (site<min_mic) site=min_mic;
                                    else if (site>max_mic) site=max_mic;
                                 }
                                 else {
                                    if (ran3(lidum)<0.5) site+=step;
                                    else site-=step;
                                 }
                              }
         					  	  	if (range_const) {
                              	//There are bouncing walls at min_mic and max_mic
                              	if (site==min_mic) ++site;
                                 else if (site==max_mic) --site;
                                 //otherwise they are free to move randomly
                                 else if (ran3(lidum)<0.5) --site;
         								else ++site;
         							}
                              else {
         								if (ran3(lidum)<0.5) --site;
         								else ++site;
                              }
         				 			break;
         	case DNA			:  drand_num=ran3(lidum);
                              //Here we implement a 95% transition bias
                              //0 and 1 : A and G
                              //2 and 3 : C and T
                              switch (site) {
                                 case 0: 	if (drand_num<trans_rate) site=1;
                                 			else
                                          	if (ran3(lidum)<0.5) site=2;
                                             else site=3;
                                          break;
                                 case 1:  if (drand_num<trans_rate) site=0;
                                 			else
                                          	if (ran3(lidum)<0.5) site=2;
                                             else site=3;
                                          break;
                                 case 2:  if (drand_num<trans_rate) site=3;
                                 			else
                                          	if (ran3(lidum)<0.5) site=0;
                                             else site=1;
                                          break;
                                 case 3:  if (drand_num<trans_rate) site=2;
                                 			else
                                          	if (ran3(lidum)<0.5) site=0;
                                             else site=1;
                                          break;
                              }
                              break;
         	case RFLP		: 	site=!site;
         							break;
         	default : break; ///DNA
         }
      }
   }

   //Continue recursion
   //testing if Nodes are true coalescing Nodes
   bool wrongNode=false;
   
   TNode* tempNode;
   if (desc1) {
      tempNode=desc1;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[curLocus] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                     tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }   
      }
      
      subst_length=mutRatePerLoc*(time-desc1->time);
      desc_mut= (int) poidev(subst_length, lidum);
      tree_length+=subst_length;
      //cout << "Curr Node " << this->ID_Node << " ID of desc1 " << desc1->ID_Node<< endl;
      tot_mut+=tempNode->add_mutations(lidum, desc_mut, curLocus,
                                                //len,
                                                numLinkedLoci, mut_type, mutRatePerLoc, geomParamPerLoc,
    													      gamma_par,trans_rate, range_const, this);
   }
   if (desc2) {
      tempNode=desc2;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[curLocus] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }    
      }
      
   	subst_length=mutRatePerLoc*(time-desc2->time);
      desc_mut= (int) poidev(subst_length, lidum);
      tree_length+=subst_length;
      //cout << "Curr Node " << this->ID_Node << " ID of desc2 " << desc2->ID_Node<< endl;
      tot_mut+=tempNode->add_mutations(lidum, desc_mut, curLocus,
                                                // len,
                                                numLinkedLoci, mut_type, mutRatePerLoc, geomParamPerLoc,
    													      gamma_par,trans_rate, range_const, this);
   }

   
   tot_mut+=cur_mut;
   return tot_mut;
}
*/
//----------------------------------------------------------------------------
//Function to add mutation for the current locus
//Inplemented for recombination
long
TNode::add_mutations(   long *            lidum,
                        const int&        num_mut,
                        const int&        curLocus,
                        //const int&        len,
      		            const int&        numLinkedLoci,
                        const Mut_Type&   mut_type,
                        const double&     mutRatePerLoc,
                        const double&     geomParamPerLoc,
                        const double&     gamma_par,
                        const double&     trans_rate,
                        const int&        range_const,
						bool&             dejaMut, //MJJ dejaMut added 11/17/08

                        const TNode*      const pN,
						
						const int& SNPtime, const int& SNPdeme, ofstream &snpfs, const int& SNPeventnum   //MJJ 11/17/08
						
						) //Pointer to the calling node (ancestor without recombination)
                        {

   //cout << "Enter TNode::add_mutations";
	//cout <<  "curlocu: " << curLocus  << " deme: " << deme << " time: " << time << endl;


   long tot_mut=0L, desc_mut, cur_mut;
   //long time=this->time;
   double subst_length;

   seq_length=numLinkedLoci;

   //bool firstTime=false;

   if (!sequence)  { //Create new sequence
   	try {
    		sequence= new TIntVect(numLinkedLoci, 0);
         int *cursite=sequence->begin();
         for (int i=0; i<numLinkedLoci; ++i, ++cursite){
            *cursite= 0; //Initialise sequence
         }
   	}
   	catch (...) {
   		if (sequence) delete sequence; sequence=NULL;
      	seq_length=0;
      	return 0;
   	}
   }
                                        
   int& cursite=*sequence->elem(curLocus);

   if (!pN) {  //Implies we are in the MRCA
      switch (mut_type) {
         case DNA:   cursite= (int) (ran3(lidum)*4); break;
         case RFLP:  cursite= (int) (ran3(lidum)*2); break;
      }
   }
   else { //Propagate the ancestral state at the current locus
      cursite=*pN->sequence->elem(curLocus);
   }

   //This is the number of new mutations to generate as compared to parent node
   num_new_mut=num_mut;

   cur_mut=0L;
   double drand_num;

   int& valueAtCurLocus=(*sequence)[curLocus];

   //Generate those mutations
   for (int i=0; i<num_mut; ++i) { //Assumes num_mut should be zero for the MRCA
      //For the moment we do not implement the gammadistributed rates, as we are assigning mutations
      //at each locus separately. We shall assume that the computation of the mutation rate
      //for this locus has been done before
      
      ++cur_mut;
      ++hits[curLocus];
   	switch (mut_type) {
      	case MICROSAT  :	//Loro_04_03_04
                           if (geomParamPerLoc>0) {
                              int step=1+geometric(geomParamPerLoc);
                              if (range_const)  {
                           	   //There are bouncing walls at min_mic and max_mic
                           	   if (valueAtCurLocus==min_mic) valueAtCurLocus+=step;
                                 else if (valueAtCurLocus==max_mic) valueAtCurLocus-=step;
                                 //otherwise they are free to move randomly
                                 else if (ran3(lidum)<0.5) valueAtCurLocus-=step;
      								   else valueAtCurLocus+=step;
                                 if (valueAtCurLocus<min_mic)valueAtCurLocus=min_mic;
                                 else if (valueAtCurLocus>max_mic)valueAtCurLocus=max_mic;
                              }
                              else {
                                 if (ran3(lidum)<0.5) valueAtCurLocus+=step;
                                 else valueAtCurLocus-=step;
                              }
                           }
                           else
      					  	  	if (range_const) {
                           	//There are bouncing walls at min_mic and max_mic
                           	if (valueAtCurLocus==min_mic) ++valueAtCurLocus;
                              else if (valueAtCurLocus==max_mic) --valueAtCurLocus;
                              //otherwise they are free to move randomly
                              else if (ran3(lidum)<0.5) --valueAtCurLocus;
      								else ++valueAtCurLocus;
      							}
                           else {
      								if (ran3(lidum)<0.5) --valueAtCurLocus;
      								else ++valueAtCurLocus;
                           }
      				 			break;
      	case DNA			:  drand_num=ran3(lidum);
                           //0 and 1 : A and G
                           //2 and 3 : C and T
                           switch (valueAtCurLocus) {
                               case 0: if (drand_num<=trans_rate) {
                                          valueAtCurLocus=1;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                          else                 valueAtCurLocus=3;
                                       }
                                       break;
                               case 1: if (drand_num<=trans_rate) {
                                          valueAtCurLocus=0;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                          else                 valueAtCurLocus=3;
                                       }
                                       break;
                               case 2: if (drand_num<=trans_rate) {
                                          valueAtCurLocus=3;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                          else                 valueAtCurLocus=1;
                                       }
                                       break;
                               case 3: if (drand_num<=trans_rate) {
                                          valueAtCurLocus=2;  
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                          else                 valueAtCurLocus=1;
                                       }
                                       break;
                           }
                           break;
				
						   //MJJ 11/13/08		   
		case SNP		:
		
		//cout << "curlocus: " << curLocus  << " deme: " << deme << " time: " << time << endl;
		//cout << "In add mutations and curlocus is: " << curLocus << endl;
		/*REMOVED TEMPORARILY MJJ 2/22/11
		                switch (valueAtCurLocus) {
                        case 0: if (drand_num<=trans_rate) {
                                        valueAtCurLocus=1;
                                }
                    		else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                        else                 valueAtCurLocus=3;
                                }
                                break;
                        case 1: if (drand_num<=trans_rate) {
                                        valueAtCurLocus=0;
                                }
                    		else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                        else                 valueAtCurLocus=3;
                                }
                                break;
                        case 2: if (drand_num<=trans_rate) {
                                        valueAtCurLocus=3;
                                }
                                else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                        else                 valueAtCurLocus=1;
                                }
                                break;
                        case 3: if (drand_num<=trans_rate) {
                                        valueAtCurLocus=2;
                                }
                    		else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                        else                 valueAtCurLocus=1;
                                }
                                break;
                }
	*/
			
			switch (valueAtCurLocus) {
				case 0:
					valueAtCurLocus=1;
					
					break;
					
				default:
					cerr << "ERROR: All ancestral SNPs should be A!" << endl;
					abort();
					break;
			
				
			}

			for (int cL=0; cL<numLinkedLoci; ++cL){
				int& valueAtCurLocus=(*sequence)[cL];
				//0 and 1 : A and G
				//2 and 3 : C and T
				if (valueAtCurLocus == 0) snpfs << "A ";
				if (valueAtCurLocus == 1) snpfs << "G ";
				if (valueAtCurLocus == 2) snpfs << "C ";
				if (valueAtCurLocus == 3) snpfs << "T ";
			}
			snpfs << endl;
		
		
		
		
						break;
						
						//MJJ 11/13/08 end
						   
						   
      	case RFLP		: 	valueAtCurLocus=!valueAtCurLocus;
      							break;
      	default : break; ///DNA
      }
   }

   //Continue recursion
   if (desc1 && desc1->flag[curLocus] ) {
      //subst_length=mut_rate*(time-desc1->time);
      subst_length=mutRatePerLoc*(time-desc1->time);
	  

// MJJ will have to add dejamut here... and something like 	  if(deme == 0 && time < 100  && !dejaMut) { //dummies for test

   //MJJ 11/14/08
   	switch (mut_type) {
	
							   		   
		case SNP		:
			
		//cout <<	"deme: " << deme << " time: " << time << endl;
		
		//cout << "Node: " << this->ID_Node << " value: " <<  valueAtCurLocus << " Desc1: " << desc1->ID_Node << endl;
		
		if(deme == SNPdeme && time <= SNPtime && !dejaMut){ //MJJ test
		snpfs << "Event "  << SNPeventnum << " curlocu: " << curLocus  << " deme: " << deme << " time: " << time << " ";
		
			
			//cout << "Event "  << SNPeventnum << " curlocu: " << curLocus  << " deme: " << deme << " time: " << time << endl;
		
		dejaMut = true;
		desc_mut= 1;
		}
		else {
			desc_mut = 0;
			}
						break;
						
						//MJJ 11/13/08 end
	
	default : 
	desc_mut= (int) poidev(subst_length, lidum);
		break; ///DNA
	
	}

//MJJ end
	  
      //desc_mut= (int) poidev(subst_length, lidum);
      tree_length+=subst_length;
      //cout << "Curr Node " << this->ID_Node << " ID of desc1 " << desc1->ID_Node<< endl;
      tot_mut+=desc1->add_mutations(lidum, desc_mut, curLocus,
                                                   //1,
                                                   numLinkedLoci, mut_type, mutRatePerLoc,
                                                   geomParamPerLoc,
       													      gamma_par,trans_rate, range_const, dejaMut, this, SNPtime, SNPdeme, snpfs, SNPeventnum);
   }

   if (desc2  && desc2->flag[curLocus]) {
      //subst_length=mut_rate*(time-desc2->time);
      subst_length=mutRatePerLoc*(time-desc2->time);
	  
	  
   //MJJ 11/14/08
   	switch (mut_type) {
	
							   		   
		case SNP		:
		
		//cout << "Node: " << this->ID_Node << " value: " <<  valueAtCurLocus << " Desc2: " << desc2->ID_Node << endl;
		
		//if(deme == SNPdeme && time <= (SNPtime + SNP_TOL) && time >= (SNPtime - SNP_TOL)  && !dejaMut){ //MJJ test
		if(deme == SNPdeme && time <= SNPtime && !dejaMut){ //MJJ test
		snpfs << "Event "  << SNPeventnum << " curlocu: " << curLocus  << " deme: " << deme << " time: " << time << " ";
			
		//cout << "Event "  << SNPeventnum << " curlocu: " << curLocus  << " deme: " << deme << " time: " << time << endl;

		dejaMut = true;
		desc_mut= 1;
		}
		else {
			desc_mut = 0;
			}
						break;
						
						//MJJ 11/13/08 end
	
	default : 
	desc_mut= (int) poidev(subst_length, lidum);
		break; ///DNA
	
	}

//MJJ end
	  
      //desc_mut= (int) poidev(subst_length, lidum);
      tree_length+=subst_length;
      //cout << "Curr Node " << this->ID_Node << " ID of desc2 " << desc2->ID_Node<< endl;
      tot_mut+=desc2->add_mutations(lidum, desc_mut, curLocus,
                                                   //1,
                                                   numLinkedLoci, mut_type, mutRatePerLoc,
                                                   geomParamPerLoc,
       													      gamma_par,trans_rate, range_const, dejaMut, this, SNPtime, SNPdeme, snpfs, SNPeventnum);
   }

   tot_mut+=cur_mut;
   return tot_mut;
}

//----------------------------------------------------------------------------
//Guillaume 02 04 2004
//Function to add mutation without recombination
//Inplemented with different mutation rates between markers
long
TNode::add_mutations_without_recombination(long *                       lidum,
                                           int*                         num_mut,
                                           int&                         curLocus,
                                           //const int&                 len,
                                           const int&                   numLinkedLoci,
                                           const TDataTypeArray*        data_type,
                                           const TMutRateArray*         mutRatePerLoc,
                                           const TMutRateArray*         geomParamPerLoc,
                                           const double&                gamma_par,
                                           const TMutRateArray*         trans_rate,
                                           const TNumLociArray*         range_const,
                                                 TNode**                rootSNPs,
                                                 bool*                  dejamut,
                                           const TNode*      const pN, 
										   
										   const int& SNPtime, const int& SNPdeme, ofstream &snpfs, const int& SNPeventnum   //MJJ 11/17/08
										   
										   ){ //Pointer to the calling node (ancestor without recombination)
                                                                                

   //cout << "Enter TNode::add_mutations";
   

   

   long tot_mut=0L, desc_mut, cur_mut;
   //long time=this->time;
   double subst_length;

   seq_length=numLinkedLoci;

   //bool firstTime=false;

   if (!sequence)  { //Create new sequence
   	try {
    		sequence= new TIntVect(numLinkedLoci, 0);
         int *cursite=sequence->begin();
         for (int i=0; i<numLinkedLoci; ++i, ++cursite){
            *cursite= 0; //Initialise sequence
         }
   	}
   	catch (...) {
   		if (sequence) delete sequence; sequence=NULL;
      	seq_length=0;
      	return 0;
   	}
   }

for(curLocus=0;curLocus<numLinkedLoci;curLocus++){



   int& cursite=*sequence->elem(curLocus);

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

   if (!pN) {  //Implies we are in the MRCA
      switch (mut_type) {
         case DNA:   cursite= (int) (ran3(lidum)*4); break;
         case RFLP:  cursite= (int) (ran3(lidum)*2); break;  
         case SNP:  
			// cursite= (int) (ran3(lidum)*4); break;
				   //rejector mod. all ancestrals are A. MJJ 11/08/08

		cursite= 0;
		break;
      }
   }
   else { //Propagate the ancestral state at the current locus
      cursite=*pN->sequence->elem(curLocus);
   }

   //This is the number of new mutations to generate as compared to parent node
   //num_new_mut=num_mut;
   num_new_mut+=num_mut[curLocus];

   cur_mut=0L;
   double drand_num;

   int& valueAtCurLocus=(*sequence)[curLocus];

   //Generate those mutations
   for (int i=0; i<num_mut[curLocus]; ++i) { //Assumes num_mut should be zero for the MRCA
      //For the moment we do not implement the gammadistributed rates, as we are assigning mutations
      //at each locus separately. We shall assume that the computation of the mutation rate
      //for this locus has been done before

      ++cur_mut;
      ++hits[curLocus];
   	switch (mut_type) {
      	case MICROSAT  :	//Loro_04_03_04
                           if (geomParamPerLoc>0) {
                              int step=1+geometric(*geomParamPerLoc->elem(curLocus));
                              if (*range_const->elem(curLocus))  {
                           	   //There are bouncing walls at min_mic and max_mic
                           	   if (valueAtCurLocus==min_mic) valueAtCurLocus+=step;
                                 else if (valueAtCurLocus==max_mic) valueAtCurLocus-=step;
                                 //otherwise they are free to move randomly
                                 else if (ran3(lidum)<0.5) valueAtCurLocus-=step;
      								   else valueAtCurLocus+=step;
                                 if (valueAtCurLocus<min_mic)valueAtCurLocus=min_mic;
                                 else if (valueAtCurLocus>max_mic)valueAtCurLocus=max_mic;
                              }
                              else {
                                 if (ran3(lidum)<0.5) valueAtCurLocus+=step;
                                 else valueAtCurLocus-=step;
                              }
                           }
                           else
                                if (*range_const->elem(curLocus)) {
                           	        //There are bouncing walls at min_mic and max_mic
                           	        if (valueAtCurLocus==min_mic) ++valueAtCurLocus;
                                        else if (valueAtCurLocus==max_mic) --valueAtCurLocus;
                                        //otherwise they are free to move randomly
                                        else if (ran3(lidum)<0.5) --valueAtCurLocus;
      								else ++valueAtCurLocus;
                                }
                                else {
      								if (ran3(lidum)<0.5) --valueAtCurLocus;
      								else ++valueAtCurLocus;
                           }
      				 			break;
      	case DNA			:  drand_num=ran3(lidum);
                           //0 and 1 : A and G
                           //2 and 3 : C and T
                           switch (valueAtCurLocus) {
                               case 0: if (drand_num<=*trans_rate->elem(curLocus)) {
                                          valueAtCurLocus=1;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                          else                 valueAtCurLocus=3;
                                       }
                                       break;
                               case 1: if (drand_num<=*trans_rate->elem(curLocus)) {
                                          valueAtCurLocus=0;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                          else                 valueAtCurLocus=3;
                                       }
                                       break;
                               case 2: if (drand_num<=*trans_rate->elem(curLocus)) {
                                          valueAtCurLocus=3;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                          else                 valueAtCurLocus=1;
                                       }
                                       break;
                               case 3: if (drand_num<=*trans_rate->elem(curLocus)) {
                                          valueAtCurLocus=2;
                                       }
                               			else {
                                        	if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                          else                 valueAtCurLocus=1;
                                       }
                                       break;
                           }
                           break;
      	case RFLP		: 	valueAtCurLocus=!valueAtCurLocus;
      							break;

        case SNP :
        if(this->ID_Node==rootSNPs[curLocus]->ID_Node && !dejamut[curLocus]){

		


		
		
                dejamut[curLocus]=true;
                //0 and 1 : A and G
                //2 and 3 : C and T
                switch (valueAtCurLocus) {
                        case 0: if (drand_num<=*trans_rate->elem(curLocus)) {
                                        valueAtCurLocus=1;
                                }
                    		else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                        else                 valueAtCurLocus=3;
                                }
                                break;
                        case 1: if (drand_num<=*trans_rate->elem(curLocus)) {
                                        valueAtCurLocus=0;
                                }
                    		else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=2;
                                        else                 valueAtCurLocus=3;
                                }
                                break;
                        case 2: if (drand_num<=*trans_rate->elem(curLocus)) {
                                        valueAtCurLocus=3;
                                }
                                else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                        else                 valueAtCurLocus=1;
                                }
                                break;
                        case 3: if (drand_num<=*trans_rate->elem(curLocus)) {
                                        valueAtCurLocus=2;
                                }
                    		else {
                             	        if (ran3(lidum)<0.5) valueAtCurLocus=0;
                                        else                 valueAtCurLocus=1;
                                }
                                break;
                }
                break;
        }
        break;
      	default : break; ///DNA
      }
   }

}


   //Continue recursion
   // if (desc1 && desc1->flag[curLocus] ) {
   if (desc1) {
      for(curLocus=0;curLocus<numLinkedLoci;curLocus++){
        //subst_length=mut_rate*(time-desc1->time);
        if(*data_type->elem(curLocus)==3) {
                num_mut[curLocus]=1;
        }
        else{
                subst_length=*mutRatePerLoc->elem(curLocus)*(time-desc1->time);
                num_mut[curLocus]= (int) poidev(subst_length, lidum);
                tree_length+=subst_length;
        }
      }
      //cout << "Curr Node " << this->ID_Node << " ID of desc1 " << desc1->ID_Node<< endl;
      tot_mut+=desc1->add_mutations_without_recombination(lidum,
                                                          num_mut,
                                                          curLocus,
                                                          //1,
                                                          numLinkedLoci,
                                                          data_type,
                                                          mutRatePerLoc,
                                                          geomParamPerLoc,
                                                          gamma_par,
                                                          trans_rate,
                                                          range_const,
                                                          rootSNPs,
                                                          dejamut,
                                                          this,  SNPtime, SNPdeme, snpfs, SNPeventnum);
   }

   //if (desc2  && desc2->flag[curLocus]) {
   if (desc2) {
      for(curLocus=0;curLocus<numLinkedLoci;curLocus++){
        if(*data_type->elem(curLocus)==3) {
                num_mut[curLocus]=1;
        }
        else{//subst_length=mut_rate*(time-desc2->time);
                subst_length=*mutRatePerLoc->elem(curLocus)*(time-desc2->time);
                num_mut[curLocus]= (int) poidev(subst_length, lidum);
                tree_length+=subst_length;
        }
      }
      //cout << "Curr Node " << this->ID_Node << " ID of desc1 " << desc1->ID_Node<< endl;
      tot_mut+=desc2->add_mutations_without_recombination(lidum,
                                                          num_mut, curLocus,
                                                          //1,
                                                          numLinkedLoci,
                                                          data_type,
                                                          mutRatePerLoc,
                                                          geomParamPerLoc,
                                                          gamma_par,
                                                          trans_rate,
                                                          range_const,
                                                          rootSNPs,
                                                          dejamut,
                                                          this,  SNPtime, SNPdeme, snpfs, SNPeventnum);
   }

   tot_mut+=cur_mut;
   return tot_mut;
}

//----------------------------------------------------------------------------
//Function to add mutation for the current locus
//Only for SNP marker (1 mutation ramdomly put somewhere on the tree)
//Inplemented for recombination
long
TNode::add_SNP_mutation(long  *            lidum,
                        const int&        num_mut,
                        const int&        curLocus,
                        //const int&        len,
                        const int&        numLinkedLoci,
                        //const Mut_Type&   mut_type,
                        const float&      gamma_par,
                        const float&      trans_rate,
                        const int&        range_const,
                        const TNode*      const pN,
                        const double&     SNP_mut_pos,
                        double&           readLength,
                        bool&             dejaMut,
                        TNode*&      SNProot) //Pointer to the calling node (ancestor without recombination)
                        {

  // cout << "Enter TNode::add_SNPmutation" << endl;
   

   

   long tot_mut=0L, desc_mut, cur_mut;
   //long time=this->time;

   double subst_length;
   seq_length=numLinkedLoci;
   //bool firstTime=false;

   if (!sequence)  { //Create new sequence
   	try {
    		sequence= new TIntVect(numLinkedLoci, 0);
         int *cursite=sequence->begin();
         for (int i=0; i<numLinkedLoci; ++i, ++cursite){
            *cursite= 9999; //Initialise sequence to 9999
         }
   	}
   	catch (...) {
   		if (sequence) delete sequence; sequence=NULL;
      	seq_length=0;
      	return 0;
   	}
   }

   int& cursite=*sequence->elem(curLocus);
   

   

   if (!pN) {  //Implies we are in the MRCA
	   //rejector mod. all ancestrals are A. MJJ 12/15/05
	   //cursite= (int) (ran3(lidum)*4);
	   cursite= 0;
   }
   else { //Propagate the ancestral state at the current locus
      cursite=*pN->sequence->elem(curLocus);
   }

   //This is the number of new mutations to generate as compared to parent node
   num_new_mut=num_mut;
   


   cur_mut=0L;
   double drand_num;


   int& valueAtCurLocus=(*sequence)[curLocus];


   //Generate those mutations
   for (int i=0; i<num_mut; ++i) { //Assumes num_mut should be zero  for the MRCA
      //For the moment we do not implement the gamma distributed rates, as we are assigning mutations
      //at each locus separately. We shall assume that the computation of the mutation rate
      //for this locus has been done before

      dejaMut=true;
      ++cur_mut;
      ++hits[curLocus];
      drand_num=ran3(lidum);
      								//Here we implement a 95% transition bias
                           //0 and 1 : A and G
                           //2 and 3 : C and T
      switch (valueAtCurLocus) {
          case 0: if (drand_num<trans_rate)  {
                     valueAtCurLocus=1;
                  }
                  else  {
                   	if (ran3(lidum)<0.5)  valueAtCurLocus=2;
                     else                 valueAtCurLocus=3;
                  }
                   break;
          case 1: if (drand_num<trans_rate)   {
                     valueAtCurLocus=0;
                  }
                  else {
                   	if (ran3(lidum)<0.5) valueAtCurLocus=2;
                     else                 valueAtCurLocus=3;
                  }
                   break;
          case 2: if (drand_num<trans_rate)   {
                     valueAtCurLocus=3;
                  }
                  else {
                   	if (ran3(lidum)<0.5) valueAtCurLocus=0;
                      else                valueAtCurLocus=1;
                  }
                  break;
          case 3: if (drand_num<trans_rate)   {
                     valueAtCurLocus=2;
                  }
                  else {
                   	if (ran3(lidum)<0.5) valueAtCurLocus=0;
                     else                 valueAtCurLocus=1;
                  }
                  break;
      }
   }

   //Continue recursion
   if (desc1 && desc2) {
      //testing if Nodes are true coalescing Nodes
      bool wrongNode=false;
      //TNode* CurrNode;
      //CurrNode=this;
      TNode* tempNode;
      if (desc1) {
         tempNode=desc1;
         wrongNode=true;
         while(wrongNode) {
            if(tempNode->flag[curLocus] ) {
               if(tempNode->desc1 && tempNode->desc2) {
                  if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) {
                     wrongNode=false;
                  }
               }
               else {
                  if (!tempNode->desc1 && !tempNode->desc2) {
                     wrongNode=false;
                  }
               }
            }
            if(wrongNode) {
               if(tempNode->desc1 && tempNode->desc2) {
                  if(tempNode->desc1->flag[curLocus]) {
                     tempNode=tempNode->desc1;
                  }
                  else {
                     tempNode=tempNode->desc2;
                  }
               }
               else {
                  if(!tempNode->desc1) {
                        tempNode=tempNode->desc2;
                  }
                  else {
                     tempNode=tempNode->desc1;
                  }
               }
            }
         }
         readLength+=time-tempNode->time;




         if((readLength > SNP_mut_pos) && !dejaMut) {
            desc_mut= 1;
            SNProot=tempNode;
         }
		          else {
            desc_mut= 0;
         }
		

         //cout << "Curr Node " << this->ID_Node << " ID of desc1 " << desc1->ID_Node<< endl;

         tot_mut+=tempNode->add_SNP_mutation(lidum, desc_mut, curLocus,
                                                //len,
                                                numLinkedLoci,
                                                //mut_type,
    						gamma_par,trans_rate, range_const, this,
                                                SNP_mut_pos, readLength, dejaMut,
                                                SNProot);
      }
      if (desc2) {
         tempNode=desc2;
         wrongNode=true;
         while(wrongNode) {
            if(tempNode->flag[curLocus] ) {
               if(tempNode->desc1 && tempNode->desc2) {
                  if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) {
                     wrongNode=false;
                  }
               }
               else {
                  if (!tempNode->desc1 && !tempNode->desc2) {
                     wrongNode=false;
                  }
               }
            }
            if(wrongNode) {
               if(tempNode->desc1 && tempNode->desc2) {
                  if(tempNode->desc1->flag[curLocus]) {
                     tempNode=tempNode->desc1;
                  }
                  else {
                     tempNode=tempNode->desc2;
                  }
               }
               else {
                  if(!tempNode->desc1) {
                     tempNode=tempNode->desc2;
                  }
                  else {
                     tempNode=tempNode->desc1;
                  }
               }
            }
         }
         readLength+=time-tempNode->time;
		 		
				 // cout << "ReadLength: " << readLength << " snp_mut_pos: " << SNP_mut_pos << endl;



         if((readLength > SNP_mut_pos) && !dejaMut) {

            desc_mut= 1;
            SNProot=tempNode;
         }
		          else {
            desc_mut= 0;
         }
	 


         //cout << "Curr Node " << this->ID_Node << " ID of desc1 " << desc1->ID_Node<< endl;
         tot_mut+=tempNode->add_SNP_mutation(lidum, desc_mut, curLocus,
                                                // len,
                                                numLinkedLoci,
                                                //mut_type,
    						gamma_par,trans_rate, range_const, this,
                                                SNP_mut_pos, readLength, dejaMut,
                                                SNProot);
      }
   }
   tot_mut+=cur_mut;
   return tot_mut;
}




//----------------------------------------------------------------------------
//A recursive routine to print tree topology and branch length
//Loro 27.8.98
void
TNode::print_tree_structure(ostream& os, const tree_output_type& tree_type, const float & mu) {
	if (desc1 && desc2) os << "(";
   if (desc1) desc1->print_tree_structure(os, tree_type, mu);
   if (desc1 && desc2) {
   	os << ", ";
   }
   else {
   	os << node_number << "." << (deme+1);
   }
   if (desc2) desc2->print_tree_structure(os, tree_type, mu);
   if (desc1 && desc2) os << ")";
   if (ancestor) {
   	switch (tree_type)  {
      	case GENERATIONS : 	os << ":" << (ancestor->time-time);
                              break;
      	case MUT_RATE    : 	os << ":" << (ancestor->time-time)*mu;
                              break;
      	case NUM_MUT     : 	os << ":" << num_new_mut;
                              break;
      }
   }
   else os << ";\n";
}
//----------------------------------------------------------------------------
//Recombination: print the coalescent tree for a given loci
void
TNode::print_tree_structure(ostream& os, const tree_output_type& tree_type, const float & mu, 
                                                                            const int curLocus, 
                                                                            const TNode* root,
                                                                            const TNode* ancNode) {
	if (desc1 && desc2) os << "(";
   bool wrongNode=false;
   TNode* tempNode;
   if (desc1) {
      tempNode=desc1;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[curLocus] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }   
      }
      tempNode->print_tree_structure(os, tree_type, mu, curLocus, root, this);  
   }



   
   if (desc1 && desc2) {
   	os << ", ";
   }
   else {
   	os << node_number << "." << (deme+1);
   }
   if (desc2) {
      tempNode=desc2;
      wrongNode=true;
      while(wrongNode) {
         if(tempNode->flag[curLocus] ) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus] && tempNode->desc2->flag[curLocus] ) { 
                  wrongNode=false;
               }
            }
            else {
               if (!tempNode->desc1 && !tempNode->desc2) {
                  wrongNode=false;
               }
            }         
         }
         if(wrongNode) {
            if(tempNode->desc1 && tempNode->desc2) {
               if(tempNode->desc1->flag[curLocus]) {
                  tempNode=tempNode->desc1;
               }
               else {
                  tempNode=tempNode->desc2;
               }
            }
            else {
               if(!tempNode->desc1) {
                  tempNode=tempNode->desc2;
               }
               else {
                  tempNode=tempNode->desc1;
               }
            }
         }    
      }
      tempNode->print_tree_structure(os, tree_type, mu, curLocus, root, this);
   }



   
   if (desc1 && desc2) os << ")";
   
   if (ancNode) {
   	switch (tree_type)  {
      	case GENERATIONS : 	os << ":" << (ancNode->time-time);
                              break;
      	case MUT_RATE    : 	os << ":" << (ancNode->time-time)*mu;
                              break;
      	case NUM_MUT     : 	os << ":" << num_new_mut;
                              break;
      }
   }
   //else os << ";\n";
   if(this==root) { 
      os << ";\n";
   }


   
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int
TListNode::allocate_list(const int& size) {

	if (list) delete[] list;

   try {
      list = new TNode*[size];
      _size=size;
   }
   catch (...) {
   	if (list) delete[] list;
      list=NULL;
      _size=0;
      cout << "TListNode::allocate_list : unable to allocate memory" << endl;
      return 0;
   }
   return 1;
}
//----------------------------------------------------------------------------
TListNode::TListNode(const int& size, TNode * tree) {
	if (tree)
	try {
		list = new TNode*[size];
      for (int i=0; i<size; ++i) {
      	tree[i].time=0;
      	tree[i].desc1=tree[i].desc2=NULL;
         list[i]=tree + i;
      }
      ListSize=size;
   }
   catch (...) {
   	if (list) delete[] list;
      list=NULL;
      _size=0;
      cout << "TListNode::TListNode(const int size, TNode * tree) : unable to allocate memory" << endl;
      return;
   }
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

TTree::TTree(const int& size) :
		chainTree(true) {
	SampleSize=size;
   delete_nodes=1;
   _mean_coalescence_time=0;
   _sd_coalescence_time=0;

   
   /***************************************************************************

   For recombination : Create a chained list of TNodes with size SampleSize

   ****************************************************************************/
   
   try {
   	//tree=new TNode[2*SampleSize];
      for(int cpt=0; cpt<SampleSize; cpt++){
         TNode* curNode=new TNode;
   	   chainTree.Add(curNode);
         curNode->ID_Node=cpt;
         curNode->event=0;
      }
   }
   catch (...) {
   	//if (tree) delete[] tree; tree=NULL;
      if (chainTree.get_beg()) chainTree.Call_deleteChain(); 
      cout << "TTree::TTree(int size) : unable to allocate memory" << endl; 
   }
}
//----------------------------------------------------------------------------
int
TTree::allocate_tree(const int& sampsize) {                           
	SampleSize=sampsize;
   //Create tree structure if needed
	//int array_size=2*SampleSize-1;
   
   //if (tree) delete[] tree;
   //try {
   //	tree=new TNode[array_size];
   //}
   //catch (...) {
   //	if (tree) delete[] tree;
   //    return 0;
   //}

   /***************************************************************************

   For recombination : Create a chained list of TNodes with size SampleSize

   ****************************************************************************/
   
   if (chainTree.get_beg()) chainTree.Call_deleteChain() ;
   try {
      for (int cpt=0; cpt<SampleSize; cpt++){
      	TNode* curNode=new TNode;
   	   chainTree.Add(curNode);
         curNode->ID_Node=cpt;
         curNode->event=0;
      }
   }
   catch (...) {
   	if (chainTree.get_beg()) chainTree.Call_deleteChain();
      return 0;
   }
   
   return 1;
}
//----------------------------------------------------------------------------
//Hudson original method
/***************************************************************************

                           For recombination 

****************************************************************************/
//This method is not possible with recombination and migration
//Since this function is never called in Simcoal I am going to remove some reference to the
//variable "tree" which is not defined with recombination

/**
int
TTree::build_tree() {
       
 	double t;

   Node list initialisation
   TListNode List(SampleSize, tree);    Cancelled by guillaume
   
   

   long lidum=1L;

   Generate the times of the nodes
   t=0.0;
   for (int i=SampleSize; i>1 ; --i) {
   	t+=-2.0*log(1.0-ran3(&lidum)) / ( (double) i*(i-1));
      tree[2*SampleSize-i].time= (long) t;     
      cout <<	(2*SampleSize-i) << "\t" << t << endl;
   }
   Generate the topology of the tree
   for (int i=SampleSize, pick, pos; i>1 ; --i) {
   	  pick= (int) (ran3(&lidum)*i);
      pos=2*SampleSize-i;
      List[pick]->ancestor=tree+pos;  
      tree[pos].desc1=List[pick];     
      List[pick]=List[i-1];
      pick= (int) ((i-1)*ran3(&lidum));
      List[pick]->ancestor=tree+pos;   
      tree[pos].desc2=List[pick];      
      List[pick]=tree+pos;      
   }
   return 1;
}
*/
//----------------------------------------------------------------------------
int
TTree::bottleneck(double time, double factor) {
	for (int i=SampleSize; i<2*SampleSize-1; ++i) {
   /***************************************************************************

      For recombination : use a chained list of TNodes 

   ****************************************************************************/
   // 	if (tree[i].time > time) {
    	if (chainTree(i).time > time) {  
   //   	cout << "\nTime before bottleneck : " << tree[i].time;
         cout << "\nTime before bottleneck : " << chainTree(i).time;
         
   //   	tree[i].time = (long) (factor*(tree[i].time - time) + time);
         chainTree(i).time = (long) (factor*(chainTree(i).time - time) + time);
         
   //   	cout << "\nTime after bottleneck  : " << tree[i].time << "\n";
         chainTree(i).time = (long) (factor*(chainTree(i).time - time) + time);
      }
   }
   return 1;
}

//----------------------------------------------------------------------------
//A routine to compute the first two moments of coalescence times and mean number
//of mutation differences over the whole tree
int
TTree::compute_moments_of_coalescence_times() {

   _mean_coalescence_time=0.0;
   _sd_coalescence_time=0.0;

	double 	sumtime=0.0, sumsquaretime=0.0,
   			num_pairwise=SampleSize*(SampleSize-1)*0.5;
   if (!num_pairwise) return 0;

	//Root the tree
   //TNode& root=tree[2*SampleSize-2];
   //for recombination take the first MRCA
   TNode& root=*MRCA_list[0];
   
   root.count_desc();  
   root.compute_moments_of_pairwise_divergence(sumtime, sumsquaretime);

   
   //_mean_coalescence_time=0.5*sumtime/num_pairwise;

   ofstream ofcoaltime("TMRCA.txt",ios::app);
   _mean_coalescence_time= (double) root.time;
   ofcoaltime <<  _mean_coalescence_time << endl;
   

   
   if (num_pairwise>1)  {
   	_sd_coalescence_time=0.5*sqrt( (sumsquaretime-sumtime*sumtime/num_pairwise)/
      										 (num_pairwise-1) );
   }

   

   
   return 1;
}

//----------------------------------------------------------------------------
//For recombination  USED FOR DEBUG ONLY
//A routine to compute the first two moments of coalescence times and mean number
//of mutation differences over the whole tree
int
TTree::compute_moments_of_coalescence_times(double** list_stat_time) {

   
   _mean_coalescence_time=0.0;
   _sd_coalescence_time=0.0;

	double 	sumtime, sumsquaretime,
   			num_pairwise=SampleSize*(SampleSize-1)*0.5;
   if (!num_pairwise) return 0;

	//Root the tree
   /***************************************************************************

      For recombination : use a chained list of TNodes
     

   ****************************************************************************/
   //TNode& root=tree[2*SampleSize-2];    // old version without recom

   
   TNode* root=MRCA_list[0];
   TNode* tempNode=root;
   for(int cpt=0;cpt<num_linked_loci-1;cpt++) {
      if(tempNode->ID_Node >= MRCA_list[cpt]->ID_Node) {
         //root=tempNode;
      }
      else {
         tempNode=MRCA_list[cpt];
      }
   }

   root->count_desc(0);
   sumtime=0.0;
   sumsquaretime=0.0;
   root->compute_moments_of_pairwise_divergence(sumtime, sumsquaretime, 0);
   _mean_coalescence_time=0.5*sumtime/num_pairwise;
   if (num_pairwise>1)  {
      _sd_coalescence_time=0.5*sqrt( (sumsquaretime-sumtime*sumtime/num_pairwise)/
                                    (num_pairwise-1) );
   }
   else {
      _sd_coalescence_time=0;
   }
     
   int num_linked_loci = root->num_linked_loci;
   double* T_MRCA                   =list_stat_time[0];
   double* within_loci_mean_time    =list_stat_time[1];
   double* within_loci_sd_time      =list_stat_time[2];
   
   double* ProdSumTime              =list_stat_time[5];
   double* sumTime                  =list_stat_time[6];
   //double* THudson                  =list_stat_time[5];
   //double* Between_loci_mean_T      =list_stat_time[4];
   //double* Between_loci_sd_T        =list_stat_time[5];


   
   within_loci_mean_time[0]=0.5*sumtime/num_pairwise;
   if (num_pairwise>1)  {
      within_loci_sd_time[0]=0.5*sqrt( (sumsquaretime-sumtime*sumtime/num_pairwise)/
                                 (num_pairwise-1) );
   }
   else {
      within_loci_sd_time[0]=0;
   }
   


   double BetweenLociMean_Time=0;
   double BetweenLociVar_Time=0;
   double BetweenLociMean_sd_Time=0;
   double BetweenLociVar_sd_Time=0;
   for(int cpt=0;cpt<num_linked_loci;cpt++) { // Compute coalescence time for every linked loci
      TNode* curr0 = (TNode*) chainTree[0];      // Reset the num_desc1 and mum_desc2
      chainTree.resetIterator();
      //cout << chainTree.get_size() << endl;
      for (int cpt1=0;cpt1<chainTree.get_size();) { 
         curr0->num_desc1=0;
         curr0->num_desc2=0;
         if (++cpt1<chainTree.get_size()) curr0=chainTree.next();
      }
      /*if(root==MRCA_list[cpt]){
      }
      else{
        root=MRCA_list[cpt];
      }*/
      if(root!=MRCA_list[cpt]){
        root=MRCA_list[cpt];
      }
      
      root->count_desc(cpt);
      sumtime=0.0;
      sumsquaretime=0.0;
      root->compute_moments_of_pairwise_divergence(sumtime, sumsquaretime, cpt);
      within_loci_mean_time[cpt]=0.5*sumtime/num_pairwise;
      if (num_pairwise>1)  {
         within_loci_sd_time[cpt]=0.5*sqrt( (sumsquaretime-sumtime*sumtime/num_pairwise)/
                                  (num_pairwise-1) );
      }
      else {
         within_loci_sd_time[cpt]=0;
      }
      sumTime[cpt]=sumtime;
      T_MRCA[cpt]=(double) MRCA_list[cpt]->time;
      BetweenLociMean_Time+=within_loci_mean_time[cpt];
      BetweenLociVar_Time+=within_loci_mean_time[cpt]*within_loci_mean_time[cpt];
      BetweenLociMean_sd_Time+=within_loci_sd_time[cpt];
      BetweenLociVar_sd_Time+=within_loci_sd_time[cpt]*within_loci_sd_time[cpt];
   }
   BetweenLociMean_Time=BetweenLociMean_Time/num_linked_loci;
   BetweenLociVar_Time=BetweenLociVar_Time/num_linked_loci;
   BetweenLociVar_Time=BetweenLociVar_Time-(BetweenLociMean_Time*BetweenLociMean_Time);
   _mean_coalescence_time=BetweenLociMean_Time;
   
   BetweenLociMean_sd_Time=BetweenLociMean_sd_Time/num_linked_loci;
   BetweenLociVar_sd_Time=BetweenLociVar_sd_Time/num_linked_loci;
   BetweenLociVar_sd_Time=BetweenLociVar_sd_Time-(BetweenLociMean_sd_Time*BetweenLociMean_sd_Time);
   _sd_coalescence_time=BetweenLociMean_sd_Time;
   
      
   for(int cpt1=0;cpt1<num_linked_loci;cpt1++) {
   ProdSumTime[cpt1]=sumTime[0]*sumTime[cpt1];
   }
   return 1;
   
}
/*
void
TTree::print_nodes(ostream& os, double factor, double time) {

	//Root the tree
	TNode& root=tree[2*SampleSize-2];

   //Create a drawing board
   int 	num_lines=root.time*UNIT_TIME+1,
   		num_col=(SampleSize-1)*WIDTH+2;
   TDrawingBoard DB(num_lines,num_col);


	//Count number of descendents starting from the root of the tree using recursions
//   cout << "\tCounting descendent nodes fron the root\n";
   root.count_desc();
//   cout << "\tDone\n";

   //Printing node information
   os << "\nNode information"
   	<< "\n================\n"
      << "\n Sample size : " << (root.num_desc1+root.num_desc2) << "\n\n";
//   cout << "\tPrinting node information, starting from the root\n";
   root.print_info(os);
//   cout << "\tDone\n";

   //Define how the tree will be drawn
   char node='o', hor_bar='-', left_corner='/', right_corner='\\', vert_bar='|';

   //Get the horizontal position of the root node : Middle of the first horizontal segment
   int node_posx=(int)((float)(root.num_desc1+root.num_desc2-1)*WIDTH/2),
   	 node_posy=0;

   //Call a recursive printing procedure from the root
   root.print_desc_nodes(DB, node_posx, node_posy, node, hor_bar, left_corner,
                          right_corner, vert_bar); 
   os << "\f\nTree structure"
   	<< "\n==============\n"
      << "Sample size : " << (root.num_desc1+root.num_desc2)
      << "\t Expansion factor : " << factor
      << "\t Expansion time : " << time
      << "\n 1 unit of time : " << UNIT_TIME << " characters\n\n";
   os << DB << endl;
}
*/
//----------------------------------------------------------------------------
void
TTree::print_nodes(ostream& os, const tree_output_type& tree_type, const float & mu) {
	//Root the tree
   /***************************************************************************

      For recombination : use a chained list of TNodes 

   ****************************************************************************/
	//TNode& root=tree[2*SampleSize-2];    // old version without recom
   TNode& root=*chainTree.last();
   root.print_tree_structure(os, tree_type, mu);
}
//Recombination: to print coalescent tree for a given linked loci
//----------------------------------------------------------------------------
void
TTree::print_nodes(ostream& os, const tree_output_type& tree_type, const float & mu, const int curLocus) {
	//Root the tree
   /***************************************************************************

      For recombination : use a chained list of TNodes 

   ****************************************************************************/
	//TNode& root=tree[2*SampleSize-2];    // old version without recom
   //TNode& root=*chainTree.last();      
   TNode& root=*MRCA_list[curLocus];      
   const TNode* Root=&root;
   root.print_tree_structure(os, tree_type, mu, curLocus, Root, NULL);
}
//----------------------------------------------------------------------------
void
TTree::print_sequences(ostream& os, const Mut_Type& mut_type) {
   /***************************************************************************

      For recombination : use a chained list of TNodes 

   ****************************************************************************/
   //int len=tree[0].seq_length;
   int len=chainTree(0).seq_length;
   char DNA_letters[4]={'A','G','C','T'};
   if (len) {
   	//Prints first the number of hits per site
   	os << "Hits\t";
		for (int j=0; j<len; ++j) {
      	//os << tree[0].hits[j] << " ";
         os << chainTree(0).hits[j] << " ";
      }
      os << "\n\n";

      //Then prints the sequences
      for (int i=0; i<SampleSize; ++i) {
   		os << (i+1) << "\t\t";
      	//if (tree[i].sequence) {          // modifidations for the recombinations
           if (chainTree(i).sequence) {
            //int* cursite=tree[i].sequence->begin();
         	int* cursite=chainTree(i).sequence->begin();
   			for (int j=0; j<len; ++j) {
            	if (mut_type==DNA)
   					//os << DNA_letters[(*tree[i].sequence)[j]];
                  os << DNA_letters[*cursite++];
               else
   					//os << (*tree[i].sequence)[j];
                  //Loro_1_3_99 : In order to avoid negative numbers
               	if (mut_type==MICROSAT) os << (500+*cursite++) << " ";
               	else os << *cursite++;
            }
      		os << "\n";
      	}
   	}
   }
   os << "\n";
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



	

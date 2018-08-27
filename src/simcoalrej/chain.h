#ifndef _CHAINED_LIST_

#define _CHAINED_LIST_

using namespace std;

//------------------------------------------------------------------------------
template<typename T>
class TSimpleChain {

	struct sTCapsule {
   	T elem;
      sTCapsule * next;
   };

	private:
		sTCapsule beg; //MJJ 10/1/12
		sTCapsule end; //MJJ 10/1/12
		long _size;

      void deleteChain();

   public:
   	TSimpleChain() {
      	beg=end=NULL;
         _size=0;
      }
      void Add(const T& e) {
      	sTCapsule* Caps=new sTCapsule;
         Caps->elem=e;
         Caps->next=NULL;
      	if (!beg) {
         	beg=Caps;
            end=Caps;
            _size=1;
         }
         else {
         	end->next=Caps;
         	end=Caps;
            ++_size;
         }
      }
      ~TSimpleChain(){
      	deleteChain();
      }
      T& last() {return end.elem;}
      T& first() {return beg.elem;}
      T& operator[] (const int& i);

};

template<typename T>
void
TSimpleChain<T>::deleteChain() {
	for (int i=0; i<_size; ++i ) {
   	if (beg) {
      	delete beg;
			if (beg->next) {
         	beg=beg->next;
            beg->next=beg->next->next;
         }
      }
   }
}


template<typename T>
T&
TSimpleChain<T>::operator[](const int&i) {
   sTCapsule cur=beg;
   if (i<_size)
		for (int j=0; j<i; ++j) {cur=cur->next;}
   if (cur) return cur.elem;
}

//------------------------------------------------------------------------------
template<typename T>
class TISimpleChain {

	struct sTCapsule {
   	T* elem;
      sTCapsule * next;
   };

	private:
		sTCapsule* beg;
		sTCapsule* end;
      T* nullElem;

		long _size;
      bool disposeElems;

      void deleteChain();
      void deleteChain(int& begin);

      sTCapsule * curItem;   //Pointer needed for iterator

   public:

      int rank_elem;

   	TISimpleChain() {
      	beg=end=NULL;
         _size=0;
         disposeElems=false;
         nullElem=NULL;
      }
   	TISimpleChain(bool dispose) {
      	beg=end=NULL;
         _size=0;
         disposeElems=dispose;
         nullElem=NULL;
      }
      ~TISimpleChain(){
      	deleteChain();
      }
      const T*& operator[] (const int& i);
      T& operator() (const int& i) {
         sTCapsule* cur=beg;
         if (i<_size) {
		      for (int j=0; j<i; ++j) {
            cur=cur->next;
            //cout << " TT " << cur << endl;
            }
         }
         if (cur) return (*cur->elem);
      }

      T* last() {return end->elem;}
      T* first() {return beg->elem;}

      sTCapsule* CurItem(){return curItem;}

      void Add(T* e) {                  // modifications
      	sTCapsule* Caps=new sTCapsule;
         Caps->elem=e;
         Caps->next=NULL;
      	if (!beg) {
         	beg=Caps;
            end=Caps;
            rank_elem=0;
            _size=1;
         }
         else {
         	end->next=Caps;
         	end=Caps;
            ++rank_elem;
            ++_size;
            //if(disposeElems == false) {
            //   cout << "DONE : Node " << Caps->elem->ID_Node << " is added " << endl ;
            //}
         }
      }

      long get_size() const { return _size; }
      sTCapsule* get_beg(){ return beg; };
      sTCapsule* get_curItem(){ return curItem; };
      void Call_deleteChain() { deleteChain(); }
      void Call_deleteChain(int& begin) { deleteChain(begin); }

      void resetIterator() {curItem=beg;}

      T* next() {
         if(curItem != end){                         // Iterator
      	curItem=curItem->next;
         return curItem->elem;
         }
         else{
            return NULL;
         }
      }
      T* current() {
         return curItem->elem;
      }

      void replace_elem(T* e) {                  // use with iterator
         //sTCapsule *temp=curItem->next;
         if(_size){
            curItem->elem=e;
         }
         else{
            sTCapsule* Caps = new sTCapsule;
            Caps->elem=e;
            Caps->next=NULL;
            beg=Caps;
            end=Caps;
            curItem=beg;
            rank_elem=0;
            _size=1;
         }
         //curItem->next=temp;
      }
      T* remove(const int& n);
      void returnTwoElemsAndReplaceFirstByNewElem(T*& e1, T*& e2, T*& e3, const int& n1,  const int& n2);
      void giveTwoElemsAndReplaceFirstByNewElem(T* e1, T* e2, T* e3);
      void giveTwoElemsAndReplaceSecondByNewElem(T* e1, T* e2, T* e3);
};

template<typename T>
void
TISimpleChain<T>::deleteChain() {
	for (int i=0; i<_size; ++i ) {
   	if (beg) {
      	if (disposeElems) {
         	delete beg->elem;
         }
         sTCapsule* temp=beg->next;
         delete beg;
         beg=temp;
      }
   }
   beg=end=NULL;
   _size=0;
}

template<typename T>
void
TISimpleChain<T>::deleteChain(int& begin) {
   sTCapsule* Curr=beg;
   sTCapsule* tempEnd;
   for(int i=0; i<begin; ++i ) {
      if(i==(begin-1)) {
         tempEnd=Curr;
         //tempEnd->elem=Curr->elem;
         //tempEnd->next=NULL;
      }
      Curr=Curr->next;
   }

   int inttemp=0;
    //cout << _size << endl;
	for (int i=begin; i<_size; ++i ) {
      if (disposeElems) {
      	delete Curr->elem;
      }
      sTCapsule* temp=Curr->next;
      delete Curr;
      Curr=temp;
      ++inttemp;
   }
   _size=_size-inttemp;
    //cout << _size << endl;
   end=tempEnd;
   end->next=NULL;

}

template<typename T>
const T*&
TISimpleChain<T>::operator[](const int&i) {
   sTCapsule* cur=beg;
   if (i<_size)
		for (int j=0; j<i; ++j) {
      cur=cur->next;
      //cout << " TT " << cur << endl;
      }
   /*if (cur) return (const T*&) cur->elem;
   return nullElem;
   */
   //Correction after g++-v3 compilation
   if (cur) {
   	return (const T*&) cur->elem;
   }
   else{
   	return (const T*&) nullElem;
   }
}

template<typename T>
T*
TISimpleChain<T>::remove(const int& n) {
	sTCapsule* cur=beg;
   --_size;
   //Check if we need to remove the first element
   if (!n) {
   	//Establish new links between capsules
    	beg=beg->next;
      T* temp=cur->elem;
      if (disposeElems) {
      	delete cur->elem;
         temp=NULL;
      }
      //Free memory of removed capsule
      delete cur;
      //Return the element that was contained in the capsule
      return temp;
   }

   //Iterate until we find the right capsule to remove
   for (int i=0; i<(n-1); ++i) {
   	cur=cur->next;
   }
   sTCapsule* temp=cur->next;
   T* tempElem=temp->elem;
   if (disposeElems) {
   	delete tempElem;
      tempElem=NULL;
   }
   //Establish new links between capsules
   cur->next=temp->next;
   if(n==_size){
      end=cur;
   }
   //Free memory of removed capsule
   delete temp;
   //Return the element that was contained in the capsule
   return tempElem;
}

template<typename T>
void
TISimpleChain<T>::returnTwoElemsAndReplaceFirstByNewElem(T*& e1, T*& e2, T*& e3,
																		  const int& n1,
                                                        const int& n2) {

   if(n1 >=_size || n2 >= _size){
      cout << n1 << " " << n2 << endl;
      int i;
      cout << "Error: object with a rank higher than the size of the chain list\n";
      cout << "Close the window (advised) or press any key to continue";
      cin >> i;
      return;
   }
	sTCapsule* cur=beg;
   //Check which index is smaller
   int first, second;
   if (n1<n2) {
   	first=n1;
      second=n2;
   }
   else {
   	first=n2;
      second=n1;
   }
   //Find first index
   for (int i=0; i<first;++i) {
   	cur=cur->next;
   }
   //Assign first element
   e1=cur->elem;
   //Replace first elem by third one
   cur->elem=e3;
   //Find second index

   for (int i=first; i<(second-1); ++i) {
   	cur=cur->next;
   }
   sTCapsule* temp=cur->next;
   //Assign second element
   e2=temp->elem;

   //Establish new link
   cur->next=cur->next->next;
   if(second==(_size-1)) {
      //cout << " *********** SIZE " << _size << endl;
      end=cur;
   }
   --_size;
   //cur->next=temp->next;  //??????
   //Destroy second capsule
   if (disposeElems) {
   	delete e2;
      e2=NULL;
   }
   delete temp;

}


template<typename T>
void
TISimpleChain<T>::giveTwoElemsAndReplaceFirstByNewElem(T* e1, T* e2, T* e3) {


	sTCapsule* cur=beg;
   sTCapsule* temp=NULL;
   sTCapsule* prev=NULL;

   //Find first Node
   bool firstNodeFound=false;
   for (int i=0; i<_size;++i) {

      if( cur->elem->ID_Node == e1->ID_Node && !firstNodeFound) { //find first node in the chain
         //Replace first elem by third one
         cur->elem=e3;
         firstNodeFound=true;
      }
      else if ( cur->elem->ID_Node == e2->ID_Node && !firstNodeFound) { //find first node in the chain
         //Replace first elem by third one
         cur->elem=e3;
         firstNodeFound=true;
      }
      else if( cur->elem->ID_Node == e1->ID_Node  && firstNodeFound) { //find second node in the chain
         //Remove capstule
         temp=cur;
         prev->next=cur->next;
         if(i == (_size-1) ) {
            end=prev;
         }
      }
      else if ( cur->elem->ID_Node == e2->ID_Node && firstNodeFound) { //find second node in the chain
         //Remove capstule
         temp=cur;
         prev->next=cur->next;
         if(i == (_size-1) ) {
            end=prev;
         }
      }
      //else {
         prev=cur;
   	   cur=cur->next;
      //}
   }
   --_size;


   //Destroy second capsule
   if (disposeElems) {
   	delete e2;
      e2=NULL;
   }
   delete temp;

}

template<typename T>
void
TISimpleChain<T>::giveTwoElemsAndReplaceSecondByNewElem(T* e1, T* e2, T* e3) {


	sTCapsule* cur=beg;
   sTCapsule* temp=NULL;
   sTCapsule* prev=NULL;

   //Find first Node
   bool firstNodeFound=false;
   for (int i=0; i<_size;++i) {
      if( cur->elem->ID_Node == e1->ID_Node  && !firstNodeFound) { //find second node in the chain
         //Remove capstule
         temp=cur;
         if(i == 0 ) {
            beg=cur->next;
         }
         else {
            prev->next=cur->next;
         }
         firstNodeFound=true;
      }
      else if ( cur->elem->ID_Node == e2->ID_Node && !firstNodeFound) { //find second node in the chain
         //Remove capstule
         temp=cur;
         if(i == 0 ) {
            beg=cur->next;
         }
         else {
            prev->next=cur->next;
         }
         firstNodeFound=true;
      }
      else if( cur->elem->ID_Node == e1->ID_Node && firstNodeFound) { //find first node in the chain
         //Replace first elem by third one
         cur->elem=e3;
      }
      else if ( cur->elem->ID_Node == e2->ID_Node && firstNodeFound) { //find first node in the chain
         //Replace first elem by third one
         cur->elem=e3;
      }
      //else {
         prev=cur;
   	   cur=cur->next;
      //}
   }
   --_size;


   //Destroy second capsule
   if (disposeElems) {
   	delete e2;
      e2=NULL;
   }
   delete temp;

}

#endif

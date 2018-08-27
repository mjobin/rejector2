#ifndef _LOCUS_H_
#define _LOCUS_H_

#include "arrays.h"
#include "cstring.h"
#include "public.h"



typedef MY_TArrayAsVector<my_string>	TLocusData;

//A class to store the results of the simulation of an independent locus
class TLocus {
   public:
      TLocusData data;           //List of alleles at that locus for all individuals in the sample
      int sampSize;
      TLocus() : data(0,10) {
         sampSize       =  0;
      }
      TLocus(const int& size): data(size,10) {sampSize=size;}
      TLocus(const TLocus& L);
};

#endif
 
//------------------------------------------------
//
// REJSTATS
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropology
// Stanford University, 2007
//
//------------------------------------------------

#include <sstream>
#include <sys/time.h>
#include "MersenneTwister.h"
#include "rdeme.h"


//










//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// STAT FXNS
//
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------







//---------------------------------------------
//Function Name:
//Function Definition:
void world::heterozygosity()
{
	
	
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_heterozygosity;
	statresult p = onewayloop("ALL", "ALL", tc);
    stats["Heterozygosity"] = p;
	
	
	
}





//---------------------------------------------
//Function Name:
//Function Definition:
void world::heterozygosityanc()
{
	
	
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_heterozygosity;
	statresult ap = onewayloop("ALL", "ANC", tc);
    stats["Heterozygosity-ANC"] = ap;
	
	
	
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::heterozygositydec()
{
	
	
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_heterozygosity;
	
	
	
	
	statresult dp = onewayloop("ALL", "DEC", tc);
    stats["Heterozygosity-DEC"] = dp;
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::heterozygosityav()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_heterozygosity;
	statresult p = onewayavloop("ALL", "ALL", tc);
    stats["HeterozygosityAv"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::heterozygosityavanc()
{
	
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_heterozygosity;
	statresult ap = onewayavloop("ALL", "ANC", tc);
    stats["HeterozygosityAv-ANC"] = ap;
	
	
	
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::heterozygosityavdec()
{
	
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_heterozygosity;
	
	
	
	
	statresult dp = onewayavloop("ALL", "DEC", tc);
    stats["HeterozygosityAv-DEC"] = dp;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition: returns heterozygosity for all heterozyg calcs
double world::tc_heterozygosity(vector<int>* y, int j)
{
	
	
    double hetz = 0.0;
	double psqrsum = 0.0;
	
	
	
	
	map<string,double> freqs = getdatafreqs(y,j);
	
	
	map<string,double>::const_iterator p = freqs.begin();
	while(p != freqs.end()){
		psqrsum += pow(p->second, 2);
		p++;
	}
	
	hetz = ((double) y->size()*(1-psqrsum))/((double) y->size()-1);//changed 2/09/06 to match what is listed in rejstats guide.. 
	
	
	/*
	 
	 vector<int>::const_iterator p = y->begin();
	 
	 while(p != y->end()){
	 cout << *p << " ";
	 string entry = getdata(*p, j);
	 cout << entry << endl;
	 
	 
	 
	 p++;
	 }
	 */
	
	
	/*
	 map<string, double>::const_iterator fp = loci[j].freqs.begin();
	 while(fp != loci[j].freqs.end()){
	 psqrsum += pow(fp->second, 2); // count up number of these alleles
	 fp++;
	 }
	 
	 hetz = (loci[j].lmap.size()*(1-psqrsum))/(loci[j].lmap.size()-1);//changed 2/09/06 to match what is listed in rejstats guide.. 
	 
	 */
	
	
	
	
	/*
	 double psqrsum = 0.0;
	 vector<multimap<string, segment*, lociCmp > >::const_iterator loci = a->lmap.begin();
	 for(int i=0; i < loc; i++){ //theres got to be a better way to do this!!!
	 loci++;			//gets the iterators pointing to the right locus
	 }
	 
	 
	 multimap<string, segment*, lociCmp >::const_iterator locusi = loci->begin();
	 while(locusi != loci->end()){ // iterate through all alleles
	 psqrsum += pow(((double)loci->count(locusi->first)/loci->size()), 2); // count up number of these alleles
	 locusi = loci->upper_bound(locusi->first); //advance to next allele
	 }
	 hetz = (a->lmap[loc].size()*(1-psqrsum))/(a->lmap[loc].size()-1);//changed 2/09/06 to match what is listed in rejstats guide.. 
	 */
	
	return hetz;
}


/*

//---------------------------------------------
//Function Name:
//Function Definition:
void world::mismatch()
{
	
	
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_mismatch;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["Mismatch"] = p;
	
	
	
}


double world::tc_mismatch(string dname){
    
    //sum of squares.. wil need to call another fxn
    

    map<int, int> mismatches;
    //first build the mismatch distn
    vector<int>::const_iterator thisun = y->begin();
	while(thisun != --y->end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator nextun = thisun;
		nextun++;
        
        double nomiss = 0.0;
		double nummatch = 0.0;
		bool miss = false;
		bool match = true;
		
		for(unsigned int i = 0; i<loci.size(); i++){
			
			for(int j=i; j<(i+loci[i].length); j++){
				
                
                if(data[*thisun][j] == 'X' || data[*nextun][j] == 'X') {
                    miss = true;
                    break;
                }
                
				if(data[*thisun][j] !=  data[*nextun][j]){ 
					match = false;	
					break;
				}
				
			}
			
			if (!miss) nomiss++;
			if (match) nummatch++;
        } 
        
        if(mismatch.find(nummatch) == mismatch.end()){
            mismatch.insert ( pair<int,int> (nummatch, 1));
        }
        else{//increment
            int theval = mismatches[nummatch];
            theval++;
            mismatches[nummatch] = theval;
        }
        thisun++;
    }
    
    //now turn this into distn

}

*/




//---------------------------------------------
//Function Name:
//Function Definition:
void world::hamming()
{
	
	
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_hamming;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["Hamming"] = p;
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition: hamming distance
double world::tc_hamming(vector<int>* y)
{
	double hamdist = 0.0;
	
	vector<double> hams;
	
	vector<int>::const_iterator thisun = y->begin();
	while(thisun != --y->end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator nextun = thisun;
		nextun++;
		
		
		double nomiss = 0.0;
		double nummatch = 0.0;
		bool miss = false;
		bool match = true;
		
		for(unsigned int i = 0; i<loci.size(); i++){
			
			for(int j=i; j<(i+loci[i].length); j++){
				
				 
				 if(data[*thisun][j] == 'X' || data[*nextun][j] == 'X') {
				 miss = true;
				 break;
				 }
				 
				if(data[*thisun][j] !=  data[*nextun][j]){ 
					match = false;	
					break;
				}
				
			}
			
			if (!miss) nomiss++;
			if (match) nummatch++;
			
			/*
			string thisdata = getdata(*thisun, i);
			
			string nextdata = getdata(*nextun, i);
			if(thisdata[0] != 'X' && nextdata[0] != 'X') nomiss++;
			
			if(thisdata == nextdata) nummatch++;; //homozygosity detected, thgus not all loci heterozygous
			 */
			
		}
		
		
		
		
		hams.push_back(nummatch/nomiss);
		
		
		
		thisun++;
		
		
	}
	
	double hamtot = 0.0;
	vector<double>::const_iterator hp = hams.begin();
	while(hp != hams.end()){
		hamtot += *hp;
		hp++;
	}
	
	hamdist = hamtot/(double)hams.size();
	
	return hamdist;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::hammingbtwn()
{
	double (world::*tc) (string, string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_hammingbtwn;
	statresult p = bytwodemeloop("ALL", "ALL", tc);
    stats["HammingBtwn"] = p;
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition: hamming distance
double world::tc_hammingbtwn(string dnamex, string dnamey)
{
	double hamdist = 0.0;
	
	vector<double> hams;
	
	vector<int>::const_iterator thisun = dmap[dnamex].begin();
	while(thisun != dmap[dnamex].end()){ //end at one before acvtual end as nextun will be the last one

		vector<int>::const_iterator otherun = dmap[dnamey].begin();
		while(otherun != dmap[dnamey].end()){

			
			double nomiss = 0.0;
			double nummatch = 0.0;
			bool miss = false;
			bool match = true;
			for(unsigned int i = 0; i<loci.size(); i++){
				
				//miss = false; 
				match = true;
				

				
				for(int j=i; j<(i+loci[i].length); j++){
					if(data[*thisun][j] == 'X' || data[*otherun][j] == 'X') {
						miss = true;
						break;
					}
					
					if(data[*thisun][j] !=  data[*otherun][j]){ 
						match = false;	
						break;
					}

				}
				//cout << endl;
				
				if (!miss) nomiss++;
				//nomiss++;
				if (match) nummatch++;
				

				
			}
			
			
			
			
			hams.push_back(nummatch/nomiss);
			
			
			otherun++;
		}
		
		thisun++;

		
		
	}
	
	double hamtot = 0.0;
	vector<double>::const_iterator hp = hams.begin();
	while(hp != hams.end()){
		hamtot += *hp;
		hp++;
	}
	
	hamdist = hamtot/(double)hams.size();
	
	return hamdist;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::numhaploblocks()
{
	
	if(dblocks.size() == 0) rsquaredblocksbydeme();
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_numhaploblocks;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["NumHaploBlocks"] = p;
	
	
	
}

double world::tc_numhaploblocks(string dname)
{
	
	return dblocks[dname].size();
	
	
}



void world::sharedhaploblocks(){
	
	if(pblocks.size() == 0) {
		vector <int> all;
		for(int i=0; i < datalength; i++){
			all.push_back(i);
		}
		
		//get haplotypes for ALL samples
		pblocks = rsquaredblocks(all);
	}
	
	
	double (world::*tc) (string, string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_sharedhaplo;
	statresult p = bytwodemeloop("ALL", "ALL", tc);
    stats["SharedHaplotypeBlocks"] = p;
	
	
	
}

double world::tc_sharedhaplo(string dx, string dy) {
	
	double sharedh = 0.0;
	
	vector<int>::const_iterator thisun = dmap[dx].begin();
	while(thisun != dmap[dx].end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator otherun = dmap[dy].begin();
		while(otherun != dmap[dy].end()){
			
			
			map<int, int>::const_iterator bp = pblocks.begin();
			while(bp != pblocks.end()){
				
				//cout << bp->first << " : " << bp->second << endl;
				
				string thisuns = makehaplo(*thisun, bp->first, bp->second);
				string otheruns = makehaplo(*otherun, bp->first, bp->second);
				
				if(thisuns == otheruns) sharedh++;
				
				bp++;
			}
			
			
			otherun++;
		}
		thisun++;
	}
	
	
	
	
	return  sharedh/(2.0 * (double) dmap[dx].size() * (double) dmap[dy].size());
}


void world::longsharedbetween(){
	
	
	double (world::*tc) (string, string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_longsharedbetween;
	statresult p = bytwodemeloop("ALL", "ALL", tc);
    stats["LongSharedBetween"] = p;
	
	
	
}

	
void world::longsharedwithin(){
	
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_longsharedwithin;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["LongSharedWithin"] = p;
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition: hamming distance
double world::tc_longsharedwithin(string dname)
{
	
	double overs = 0.0;
	double tested = 0.0;
	
	vector<int>::const_iterator thisun = dmap[dname].begin();
	while(thisun != dmap[dname].end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator otherun = dmap[dname].begin();
		while(otherun != dmap[dname].end()){
			
			//cout << endl << endl << "Testing " << *thisun << " and " << *otherun << endl;
			
			if (*thisun != *otherun) { //dont compre from same sequence
				
			
			
			
			
			int longestrun = 0;
			
			bool overthreshold = false;
				
			bool sharedrun = false;	
			
			for(unsigned int i = 0; i<loci.size(); i++){
				
				
				
				bool shared = true;
				
				for(int j=i; j<(i+loci[i].length); j++){
					if(data[*thisun][j] == 'X' || data[*otherun][j] == 'X') {
						shared = false;
					}
					
					if(data[*thisun][j] !=  data[*otherun][j]){ 
						shared = false;
					}
					
							
					
				}
				
				
				if(shared){
					
					
					if (i > 0 && sharedrun) {//dont check the first one
							if(loci[(i-1)].phdist < MAX_LOCI_JUMP){//do not include jumps over centromeres and to nbew chromosomes
								longestrun += loci[(i-1)].phdist; //keepadding
								if(longestrun >= LONGBLOCK_THRESHOLD)overthreshold = true;
								//if(longestrun > 3000000)cout << dname << ": " << *thisun << " and " << *otherun << ".   Run of " << longestrun << " by " << i << endl;
							}
						else { //if jumped over centromere/chromsoome, start over
							
							
							longestrun = 0;
							sharedrun = false;
						}

		

					}
					else sharedrun = true; //if shared on the first one, just pass on that a run has started
				}
				else{
					
					longestrun = 0;
					sharedrun = false;
				}
					
				

				if(overthreshold) break;
			}
			
			tested++;
			
			if(overthreshold) overs++;
			
				//end if they dont' match
			}
			
			
			otherun++;
		
		

		
		}
		
		thisun++;
		
		
		
	}
	
	return overs/tested;
}


//---------------------------------------------
//Function Name:
//Function Definition: hamming distance
double world::tc_longsharedbetween(string dnamex, string dnamey)
{

	double overs = 0.0;
	double tested = 0.0;
	
	vector<int>::const_iterator thisun = dmap[dnamex].begin();
	while(thisun != dmap[dnamex].end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator otherun = dmap[dnamey].begin();
		while(otherun != dmap[dnamey].end()){
			
			//cout << endl << endl << "Testing " << *thisun << " and " << *otherun << endl;

			int longestrun = 0;
			
			bool overthreshold = false;
			
			bool sharedrun = false;	
			
			for(unsigned int i = 0; i<loci.size(); i++){
				
				
				
				bool shared = true;
				
				for(int j=i; j<(i+loci[i].length); j++){
					if(data[*thisun][j] == 'X' || data[*otherun][j] == 'X') {
						shared = false;
					}
					
					if(data[*thisun][j] !=  data[*otherun][j]){ 
						shared = false;
					}
					
				
					

					
				}
				
				if(shared){
					
					
					if (i > 0 && sharedrun) {//dont check the first one
						if(loci[(i-1)].phdist < 500000){//do not include jumps over centromeres and to nbew chromosomes
							longestrun += loci[(i-1)].phdist; //keepadding
							if(longestrun >= LONGBLOCK_THRESHOLD)overthreshold = true;
							//if(longestrun > 3000000)cout << dname << ": " << *thisun << " and " << *otherun << ".   Run of " << longestrun << " by " << i << endl;
						}
						else { //if jumped over centromere/chromsoome, start over
							
							
							longestrun = 0;
							sharedrun = false;
						}
						
						
						
					}
					else sharedrun = true; //if shared on the first one, just pass on that a run has started
				}
				else{
					
					longestrun = 0;
					sharedrun = false;
				}
				
				
				
				if(overthreshold) break;
				
				
			}
			
			tested++;
			
			if(overthreshold) overs++;

			
			
			otherun++;
		}
		
		thisun++;
		
		
		
	}
	
	return overs/tested;
}

/*
double world::tc_longshared(string dx, string dy) {
	
	double sharedlong = 0.0;
	double tested = 0.0;
	
	vector<int>::const_iterator thisun = dmap[dx].begin();
	while(thisun != dmap[dx].end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator otherun = dmap[dy].begin();
		while(otherun != dmap[dy].end()){
			
			tested++;
			
			map<int, int>::const_iterator bp = pblocks.begin();
			while(bp != pblocks.end()){
				
				cout << bp->first << " : " << bp->second << endl;
				
				string thisuns = makehaplo(*thisun, bp->first, bp->second);
				string otheruns = makehaplo(*otherun, bp->first, bp->second);
				
				if(thisuns == otheruns){
					
					double totlength = 0.0;

						for(int i=bp->first; i < bp->second; i++){ //inclusive
							totlength += loci[i].phdist;
							//end for i	
						}

					if(totlength >= LONGBLOCK_THRESHOLD) {
						sharedlong++;
						break;
					}
					
					
					
					
				}
				
				bp++;
			}
			
			
			otherun++;
		}
		thisun++;
	}
	
	
	
	
	return  sharedlong/tested;
}
*/


//---------------------------------------------
//Function Name:
//Function Definition:
void world::numroh()
{
	
	if(dslideroh.size() == 0) rohslidebydeme();
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_numroh;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["NumROH"] = p;
	
	
	
}

double world::tc_numroh(string dname)
{

	
	double tot=0.0;
	
	map <int, vector<bool> >::const_iterator cp = dslideroh[dname].begin();
	while(cp != dslideroh[dname].end()){
		
		tot+=((double)cp->second.size()); //How many rohs are in that individual

		cp++;
	}
	
	
	return tot/(double) dslideroh[dname].size();
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::froh()
{
	
	if(dslideroh.size() == 0) rohslidebydeme();
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_froh;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["FROH"] = p;
	
	
	
}

double world::tc_froh(string dname)
{
	
	
	
	if(dslideroh[dname].size() == 0) return -999;
	
	double tot = 0.0;
	
	
	map <int, vector<bool> >::const_iterator dp = dslideroh[dname].begin();
	
	while(dp != dslideroh[dname].end()){
		
		//cout << dp->first << ":: " << endl;
		unsigned int i =0;
		double lroh = 0.0; 
		double lauto = 0.0;
		
		if(dp->second.size() != loci.size()) return -999;
		
		vector<bool>::const_iterator dpp = dp->second.begin();
		while(dpp != dp->second.end()){
			
			if(*dpp) {
				
				lroh += loci[i].phdist;
				
			}
			lauto += loci[i].phdist;
			i++;
			dpp++;
		}
		
		tot += lroh/lauto; //going to average all the froh for a deme
		//cout << endl;
		dp++;	
	}
	
	
	return tot/(double)dslideroh[dname].size();
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::frohvar()
{
	
	if(dslideroh.size() == 0) rohslidebydeme();
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_frohvar;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["FROHVar"] = p;
	
	
	
}

double world::tc_frohvar(string dname)
{
	
	
	
	if(dslideroh[dname].size() == 0) return -999;
	
	double tot = 0.0;
	
	vector<double> frohs;
	
	map <int, vector<bool> >::const_iterator dp = dslideroh[dname].begin();
	
	while(dp != dslideroh[dname].end()){
		
		//cout << dp->first << ":: " << endl;
		unsigned int i =0;
		double lroh = 0.0; 
		double lauto = 0.0;
		
		if(dp->second.size() != loci.size()) return -999;
		
		vector<bool>::const_iterator dpp = dp->second.begin();
		while(dpp != dp->second.end()){
			
			if(*dpp) {
				
				lroh += loci[i].phdist;
				
			}
			lauto += loci[i].phdist;
			i++;
			dpp++;
		}
		
		frohs.push_back(lroh/lauto);
		tot += lroh/lauto;
		//cout << endl;
		dp++;	
	}
	
	//mean
	
	double dmean = tot/(double)dslideroh[dname].size();
	//cout << "Mean: " << dmean << endl;
	
	double sumsq = 0.0;
	
	vector<double>::const_iterator fp = frohs.begin();
	while(fp != frohs.end()){
	//	cout << *fp << endl;
		sumsq += pow((*fp - dmean),2.0);
		fp++;
	}
	
	
	return(sumsq/(frohs.size()-1));
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::varrohlength()
{
	
	if(dslideroh.size() == 0) rohslidebydeme();
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_varrohlength;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["VarROHLength"] = p;
	
	
	
}

double world::tc_varrohlength(string dname)
{
 
    if(dslideroh[dname].size() == 0) return -999;
    
    double tot = 0.0;
    vector<double> lrohs;
    
    
    map <int, vector<bool> >::const_iterator dp = dslideroh[dname].begin();
	
	while(dp != dslideroh[dname].end()){
        
        unsigned int i =0;
        double lroh = 0.0; 
        
        if(dp->second.size() != loci.size()) return -999;
        
        vector<bool>::const_iterator dpp = dp->second.begin();
		while(dpp != dp->second.end()){
            
            if(*dpp) {
				lroh += loci[i].phdist;
				
			}
            
            dpp++;
            i++;
            
        }
        
        lrohs.push_back(lroh);
        tot += lroh;
        
        dp++;
        
    }
    
    //mean
	
	double dmean = tot/(double)dslideroh[dname].size();
	//cout << "Mean: " << dmean << endl;
	
	double sumsq = 0.0;
	
	vector<double>::const_iterator lp = lrohs.begin();
	while(lp != lrohs.end()){
        //	cout << *fp << endl;
		sumsq += pow((*lp - dmean),2.0);
		lp++;
	}
	
	
	return(sumsq/(lrohs.size()-1));
    
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::blocklength()
{
	
	if(dblocks.size() == 0) rsquaredblocksbydeme();
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_blocklength;
	statresult p = bydemeloop("ALL", "ALL", tc);
    stats["BlockLength"] = p;
	
	
	
}

double world::tc_blocklength(string dname)
{
	
	int totlength = 0;
	
	
	map<int, int>::const_iterator dp = dblocks[dname].begin();
	
	if(dp == dblocks[dname].end()) return -999; //send an NA if the block list is empty
	
	while(dp != dblocks[dname].end()){
		for(int i=dp->first; i < dp->second; i++){ //inclusive
			totlength += loci[i].phdist;
			//end for i	
		}
		dp++;
	}
	
	
	
	return totlength/(double)dblocks[dname].size();
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::blockboundary()
{
	
	if(dblocks.size() == 0) rsquaredblocksbydeme();
	
	double (world::*tc) (string, string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_blockboundary;
	statresult p = bytwodemeloop("ALL", "ALL", tc);
    stats["BlockBoundary"] = p;
	
	
	
}

double world::tc_blockboundary(string dnamex, string dnamey)
{
	
	double corr = 0.0;
	
	
	vector<bool> x;
	
	for(unsigned int i = 0; i<loci.size(); i++){
		bool inablock = false;
		map<int, int>::iterator bp = dblocks[dnamex].begin(); //begin at first place its less than or equal to sp second
		while(bp != dblocks[dnamex].end()){
			if(i >= bp->first && i <= bp->second){ //it is inside the block, so no new block need be created
				inablock = true;
				break;
				//end if nin a block
			}
			
			
			bp++;
		}
		if(inablock) x.push_back(true);
		else x.push_back(false);
	}
	
	
	
	
	
	vector<bool> y;
	
	for(unsigned int j = 0; j<loci.size(); j++){
		bool inablock = false;
		map<int, int>::iterator bp = dblocks[dnamey].begin(); //begin at first place its less than or equal to sp second
		while(bp != dblocks[dnamey].end()){
			if(j >= bp->first && j <= bp->second){ //it is inside the block, so no new block need be created
				inablock = true;
				break;
				//end if nin a block
			}
			
			
			bp++;
		}
		if(inablock) y.push_back(true);
		else y.push_back(false);
	}
	
	if(x.size() != y.size()) return -999;
	
	double matches = 0.0;
	
	vector<bool>::const_iterator xp = x.begin();
	vector<bool>::const_iterator yp = y.begin();
	while (xp != x.end()&& yp != y.end()) {
		if (*xp == *yp) matches++;
		xp++;
		yp++;
	}
	
	
	corr = matches/(double)x.size();
	
	return corr;
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
/*
 void world::blockdensitycorr()
 {
 
 if(dblocks.size() == 0) rsquaredblocksbydeme();
 
 double (world::*tc) (string, string) = NULL; //fxn pointer with list of indivs and locus to use
 tc = &world::tc_blockdensitycorr;
 statresult p = bytwodemeloop("ALL", "ALL", tc);
 stats["BlockDensity"] = p;
 
 
 
 }
 
 double world::tc_blockdensitycorr(string dnamex, string dnamey)
 
 {
 double corr = 0.0;
 
 double xnum = (double)dblocks[dnamex].size();
 
 
 retrun corr;
 
 }
 */

// BELOW DOES NOT WORK WITH SLIDING WINDOW ROH'S YET
/*
//---------------------------------------------
//Function Name:
//Function Definition:
void world::rohcorr()
{
	
	if(droh.size() == 0) rohbydeme();
	
	double (world::*tc) (string, string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_rohcorr;
	statresult p = bytwodemeloop("ALL", "ALL", tc);
    stats["ROHCorr"] = p;
	
	
	
}

double world::tc_rohcorr(string dx, string dy)
{
	
	double corr = -999;
	if(droh[dx].size() != droh[dy].size()) return -999;
	
	double matches = 0.0;
	
	vector<bool>::const_iterator xp = droh[dx].begin();
	vector<bool>::const_iterator yp = droh[dy].begin();
	while (xp != droh[dx].end()&& yp != droh[dy].end()) {
		if (*xp == *yp) matches++;
		xp++;
		yp++;
	}
	
	
	corr = matches/(double)droh[dx].size();
	
	return corr;
	
}

*/

//---------------------------------------------
//Function Name: Number of segregating sites
//Function Definition:
void world::numsegsites(){
	
	
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_numsegsites;
	statresult p = bydemeloop("ALL", "ALL", tc);
     stats["NumSegSites"] = p;
	

}

//---------------------------------------------
//Function Name:
//Function Definition: returns heterozygosity for all heterozyg calcs
double world::tc_numsegsites(string dname){
	
	return (double) countsegsites(&dmap[dname]);
	
		
		
		
	
}


//---------------------------------------------
//Function Name: Number of segregating sites
//Function Definition:
void world::tajimasd(){
	
	double (world::*tc) (string) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_tajimasd;
	statresult p = bydemeloop("ALL", "ALL", tc);
	stats["TajimasD"] = p;
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition: returns heterozygosity for all heterozyg calcs
double world::tc_tajimasd(string dname){
	
	//a la Statistical Methodfor Testing the Neutral Mutation Hypothesis by DNA Polymorphism, Tajima 1989
	
	double khat = meanpairwise(&dmap[dname]);
	
	
	double caps = countsegsites(&dmap[dname]);
	
	
	double aone = 0.0;
	double atwo = 0.0;
	
	for(unsigned int i=0; i < (dmap[dname].size()-1); i++){
		aone +=  (1/((double)i+1));
		atwo += pow((1/((double)i+1)), 2);
	}
	
	
	double ssz = (double)dmap[dname].size();
	
	double bone = (ssz+1.0)/(3.0*(ssz-1.0));
	
	double btwo = (2.0*(pow(ssz,2.0)+ssz+3.0))/((9.0*ssz) * (ssz-1.0));
	
	double cone = bone - (1.0/aone);
	
	double ctwo = btwo - ((ssz+2.0)/(aone*ssz)) + atwo/(pow(aone,2.0));
	
	double eone = cone/aone;
	
	double etwo = ctwo/(pow(aone,2.0)+atwo);
	 
	
	
	double dee = khat - (caps/aone);
	
	double veedee = pow( (eone*caps)+((etwo*caps)*(caps-1.0)), 0.5);
	
	double capd = dee/veedee;
	
	
	return capd;
	

	
	
	
	
}




//---------------------------------------------
//Function Name:
//Function Definition:
void world::allheterozygotes()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_allhet;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["CompleteHeterozygotes"] = p;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_allhet(vector<int>* y)
{
	int hetcount = 0;
	
	
	//go through and grab all double entries. i.e. two chromosomes in segs listing with same name
	vector<int>::const_iterator thisun = y->begin();
	while(thisun != --y->end()){ //end at one before acvtual end as nextun will be the last one
		bool hetstate = false;
		vector<int>::const_iterator nextun = thisun;
		nextun++;
		
		
		//cout << "thisun:" << thisun->second->name << " nextun:" << nextun->second->name << endl;
		
		
		
		if(nmap[*thisun] == nmap[*nextun]){//they have same names thus are two chromosomes from same individual
			hetstate = true;
			
			for(unsigned int i = 0; i<loci.size(); i++){
				
				
	
				
				
				string thisdata = getdata(*thisun, i);
				if(thisdata[0] == 'X') hetstate = false;//assume that missing data anywhere will invalidate test and assume can't SHOW its allhet 2/9/06
				string nextdata = getdata(*nextun, i);
				if(nextdata[0] == 'X') hetstate = false;//assume that missing data anywhere will invalidate test and assume can't SHOW its allhet 2/9/06
				if(thisdata == nextdata) hetstate = false; //homozygosity detected, thgus not all loci heterozygous
				
			}
			
			
			
		}
		
		
		//cout << "State (0 is false): " << thisunstate << endl << endl;
		if(hetstate) hetcount++;
		
		thisun++;
		
		
	}
	
	
	return (double) hetcount;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::allhomozygotes()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_allhom;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["CompleteHomozygotes"] = p;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_allhom(vector<int>* y)
{
	int homcount = 0;
	
	
	//go through and grab all double entries. i.e. two chromosomes in segs listing with same name
	vector<int>::const_iterator thisun = y->begin();
	while(thisun != --y->end()){ //end at one before acvtual end as nextun will be the last one
		bool homstate = false;
		vector<int>::const_iterator nextun = thisun;
		nextun++;
		
		
		//cout << "thisun:" << thisun->second->name << " nextun:" << nextun->second->name << endl;
		
		
		
		if(nmap[*thisun] == nmap[*nextun]){//they have same names thus are two chromosomes from same individual
			homstate = true;
			
			for(unsigned int i = 0; i<loci.size(); i++){
				string thisdata = getdata(*thisun, i);
				if(thisdata[0] == 'X') homstate = false;//assume that missing data anywhere will invalidate test and assume can't SHOW its allhet 2/9/06
				string nextdata = getdata(*nextun, i);
				if(nextdata[0] == 'X') homstate = false;//assume that missing data anywhere will invalidate test and assume can't SHOW its allhet 2/9/06
				if(thisdata != nextdata) homstate = false; //homozygosity detected, thgus not all loci heterozygous
				
			}
			
			
			
		}
		
		
		//cout << "State (0 is false): " << thisunstate << endl << endl;
		if(homstate) homcount++;
		
		thisun++;
		
		
	}
	
	
	return (double) homcount;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::haplotypes()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_haplo;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["Haplotypes"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_haplo(vector<int>* y)
{
	int gh = hcounter(y, 0);
	
	
	return (double) gh;
	
}

//---------------------------------------------
//Function Name:
//Function Definition: recursive function
int world::hcounter(vector<int>* y, int j)
{
	unsigned int gh = 0;
	
	if((j+1) == loci.size()){ 
		map<string, int> counts = getdatacounts(y, j);
		gh += counts.size();
		
	}
	
	
	else{ //if not make a subsystem and go recursive
		
		set<string> alleles = getalleles(y,j);
		set<string>::const_iterator ap = alleles.begin();
		while(ap != alleles.end()){
			
			vector<int> z = getsubset(y, j, *ap);
			gh += hcounter(&z, (j+1));
			
			ap++;
		}
		
	}
	
	return gh;
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::numdiffalleles()
{
   	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_ndalleles;
	statresult p = onewayloop("ALL", "ALL", tc);
    stats["NumDiffAlleles"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::numdiffallelesav()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_ndalleles;
	statresult p = onewayavloop("ALL", "ALL", tc);
    stats["NumDiffAllelesAv"] = p;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_ndalleles(vector<int>* y, int j)
{
	
	set<string> alleles = getalleles(y,j);
	
	return (double) alleles.size();
	
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::samplesize()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_sampsz;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["SampleSize"] = p;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::samplesizeanc()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_sampsz;
	statresult ap = persystemloop("ALL", "ANC", tc);
    stats["SampleSize-ANC"] = ap;
	
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::samplesizedec()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_sampsz;
	
	
	statresult dp = persystemloop("ALL", "DEC", tc);
    stats["SampleSize-DEC"] = dp;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::meanlength()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_meanl;
	statresult p = onewayloop("ALL", "ALL", tc);
    stats["MeanLength"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_meanl(vector<int>* y, int j) 
{
	return(locmean(y,j));
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::minlength()
{
   	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_minl;
	statresult p = onewayloop("MICROSAT", "ALL", tc);
    stats["MinLength"] = p;
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_minl(vector<int>* y, int j) 
{
	
	
	double min = -999;
	
	vector<int>::const_iterator p = y->begin();
	
	while(p != y->end()){
		double thisun = getdataasdbl(*p, j);
		if(min == -999){
			if(thisun != -999) min = thisun;
		}
		else{
			if(thisun != -999 && thisun < min) min = thisun;
		}
		p++;
	}
	
	return min;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::maxlength()
{
   	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_maxl;
	statresult p = onewayloop("MICROSAT", "ALL", tc);
    stats["MaxLength"] = p;
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_maxl (vector<int>* y, int j) 
{
	
	
	double max = -999;
	
	vector<int>::const_iterator p = y->begin();
	
	while(p != y->end()){
		double thisun = getdataasdbl(*p, j);
		if(max == -999){
			if(thisun != -999) max = thisun;
		}
		else{
			if(thisun != -999 && thisun > max) max = thisun;
		}
		p++;
	}
	
	return max;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::rangeoflocus()
{
    double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_range;
	statresult p = onewayloop("MICROSAT", "ALL", tc);
    stats["RangeofLocus"] = p;
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::rangeoflocusav()
{
    double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_range;
	statresult p = onewayavloop("MICROSAT", "ALL", tc);
    stats["RangeofLocusAv"] = p;
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::nucleotidediversity()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_nucdiv;
	statresult p = onewayloop("DNA", "ALL", tc);
    stats["NucleotideDiversity"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::nucleotidediversityav()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_nucdiv;
	statresult p = onewayavloop("DNA", "ALL", tc);
    stats["NucleotideDiversityAv"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition: calculates nucelotide diveristy a la Molecular Evo Genetics by Nei p. 256
double world::tc_nucdiv(vector<int>* y, int j) 
{
	
	double tot = 0.0;
	map<string, double> freqs =getdatafreqs(y, j);
	
	map<string, double>::const_iterator li = freqs.begin();
	map<string, double>::const_iterator lj;
	while(li != freqs.end()){
		map<string, double>::const_iterator lj = freqs.begin();
		while(lj != freqs.end()){
			tot += (li->second * lj->second * difnuc(li->first, lj->first)); 
			
			lj++;
		}
		
		
		li++;
	}
	
	double pihat = (((double)y->size()/((double)y->size() - 1.0))) * tot;
	
	return pihat;
	
}

//---------------------------------------------
//Function Name:
//Function Definition: counts the number of differences between two sequences and returns the proportion of nucelotides that are different
double world::difnuc(string w, string x)
{
	if(w.size() != x.size()){
		cerr << "Error in difnuc function! DNA sequence " << w << " is not equal in length to sequence " << x << endl;
		abort();
	}
	
	unsigned int numdif = 0;
	
	unsigned int strsz = w.size();
	
	string::const_iterator pa = w.begin();
	string::const_iterator pb = x.begin();
	
	while(pa != w.end() && pb != x.end()){
		
		//cout << "difnuc: " << *pa << "   " << *pb << endl;
		
		if(*pa != *pb) numdif++;
		
		pa++;
		pb++;
	}
	
	double difprop = (double)numdif/(double)strsz;
	return difprop;
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::Fstdist()
{
	double (world::*tc) (int) = NULL;
	tc = &world::tc_fst;
	statresult p = perlocusxpoploop("ALL", "ALL", tc);
    stats["Fst"] = p;
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::Fstdistav()
{
	double (world::*tc) (int) = NULL;
	tc = &world::tc_fst;
	statresult p = perlocusxpopavloop("ALL", "ALL", tc);
    stats["FstAv"] = p;
}				




//---------------------------------------------
//Function Name:
//Function Definition:
double world::averageheterozygosity(int j)
{
	double tot = 0.0;
	map<string, vector<int> >::iterator dp = dmap.begin();
	while(dp != dmap.end()){
		vector<int> dl = getintersection(&dp->second, &loci[j].lmap);
		tot += tc_heterozygosity(&dl,j);
		
		dp++;
	}
	
	return tot/(double) dmap.size();
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_fst(int j) 
{
	
	
	
	
	double alldemehet = tc_heterozygosity(&loci[j].lmap, j);
	double hetz = ((alldemehet-averageheterozygosity(j))/alldemehet);
	if(hetz<0) return 0; //floor of 0
	return hetz;
	
	
	
	
	
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::betaimbalance()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_beta;
	statresult p = onewayloop("MICROSAT", "ALL", tc);
    stats["Beta"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::betaimbalanceanc()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_beta;
	statresult ap = onewayloop("MICROSAT", "ANC", tc);
    stats["Beta-ANC"] = ap;
	
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::betaimbalancedec()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_beta;
	
	
	statresult dp = onewayloop("MICROSAT", "DEC", tc);
    stats["Beta-DEC"] = dp;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::betaimbalanceav()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_beta;
	statresult p = onewayavloop("MICROSAT", "ALL", tc);
    stats["BetaAv"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::betaimbalanceavanc()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_beta;
	statresult ap = onewayavloop("MICROSAT", "ANC", tc);
    stats["BetaAv-ANC"] = ap;
	
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::betaimbalanceavdec()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_beta;
	
	
	statresult dp = onewayavloop("MICROSAT", "DEC", tc);
    stats["BetaAv-DEC"] = dp;
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_beta(vector<int>* y, int j)
{
	double Vt = 0.0;
	double P0t = 0.0;
	double beta = 0.0;
	
	double dememean = locmean(y, j);
	
	
	vector<int>::const_iterator yp = y->begin();
	while(yp != y->end()){
		double xi = atof(getdata(*yp, j).c_str());
		Vt += (pow((xi - dememean),2)) / (y->size() - 1);
		yp++;
	}
	
	/*
	 map<string, double> freqs = getdatafreqs(y,j);
	 map<string, double>::const_iterator ap = freqs.begin();
	 while(ap != freqs.end()){
	 cout << "Freq: " << ap->first << endl;
	 double xi = atof(ap->first.c_str());
	 Vt += (pow((xi - dememean),2)) / (y->size() - 1);
	 ap++;
	 }
	 */
	
	
	/*
	 vector<multimap<string, segment*, lociCmp > >::const_iterator locip = a->lmap.begin();
	 for(int i=0; i < loc; i++){ //theres got to be a better way to do this!!!
	 locip++;			//gets the iterators pointing to the right locus
	 }
	 
	 multimap<string, segment*, lociCmp >::const_iterator locusi = locip->begin();
	 while(locusi != locip->end()){ // iterate through all alleles
	 double xi = atof(locusi->second->loci[loc].c_str());
	 Vt += (pow((xi - dememean),2)) / (a->lmap[loc].size() - 1);
	 locusi++;
	 }
	 */
	
	
	
	Vt = Vt*2;
	double sumpk = 0.0;
	
	map<string, double> freqs = getdatafreqs(y,j);
	map<string, double>::const_iterator ap = freqs.begin();
	while(ap != freqs.end()){
		sumpk += pow(ap->second, 2); // count up number of these alleles
		ap++;
	}
	
	/*
	 locusi = locip->begin();
	 while(locusi != locip->end()){ // iterate through all alleles
	 
	 sumpk += pow(((double)locip->count(locusi->first)/a->lmap[loc].size()), 2); // count up number of these alleles
	 
	 locusi = locip->upper_bound(locusi->first); //advance to next allele
	 }
	 */
	
	sumpk = (sumpk * y->size()) -1.0;
	
	
	
	//sumpk = (sumpk * a->lmap[loc].size()) -1;
	
	P0t = sumpk/(y->size() - 1.0);
	
	//P0t = sumpk/(a->lmap[loc].size() - 1);
	
	beta = Vt/(((1/(pow(P0t,2))) - 1)/2);
	
	
	
	return beta;
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::microsatvariance()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_micsatvar;
	statresult p = onewayloop("MICROSAT", "ALL", tc);
    stats["MicrosatVariance"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::microsatvarianceav()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_micsatvar;
	statresult p = onewayavloop("MICROSAT", "ALL", tc);
    stats["MicrosatVarianceAv"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_micsatvar(vector<int>* y, int j)
{
	
	
	return locvar(y,j);
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::derivedfraction()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_derfrac;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["DerivedFraction"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_derfrac(vector<int>* y)
{
	
	vector<int> decvec  = getancdecvec(y, false);
	return((double)decvec.size()/(double)y->size());
	
	
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::minorallelefreq()
{
	
	double (world::*tc) (vector<int>*) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_maf;
	statresult p = persystemloop("ALL", "ALL", tc);
    stats["MinorAlleleFrequency"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_maf(vector<int>* y)
{
	
	
	double maf;
	
	for(unsigned int i=0; i<loci.size(); i++){
		map<string, int> counts = getdatacounts(y, i);
		int locount = 0;
		bool firsttime = true;
		int totcount = 0;
		if(counts.size() > 1){
			map<string, int>::const_iterator cp = counts.begin();
			while (cp != counts.end()) {
				totcount += cp->second;
				if (firsttime) {
					locount = cp->second;
				}
				else {
					if(cp->second < locount) locount = cp->second;
				}
				firsttime = false;
				cp++;
			}
			
			maf += (double) locount/(double)totcount;
		}
		//cout << (double) locount/(double)totcount << endl;
		
	}
	
	
	
	return maf/(double)loci.size(); //avergaed by loci
}




//---------------------------------------------
//Function Name:
//Function Definition:
void world::jstatistic()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_jstat;
	statresult p = onewayloop("ALL", "ALL", tc);
    stats["Jstatistic"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::jstatisticav()
{
	
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_jstat;
	statresult p = onewayavloop("ALL", "ALL", tc);
    stats["JstatisticAv"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::modjstatisticav()
{
    
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_jstat;
	statresult p = modonewayavloop("ALL", "ALL", tc);
    stats["ModJstatisticAv"] = p;
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_jstat(vector<int>* y, int j)
{
	
	vector<int> decvec  = getancdecvec(y, false);
	if(decvec.empty()) return -999;
	double lvar = locvar(&decvec,j);
	
	
	
	
	//double lvar = dec.locvar(loc);
	
	//fraction of derived snps GOING WITH SEGS.SIZE BECAUSE MISSING DATA AT MICROSAT WILL STILL BE COUNTED AT the SNP FOR DETERMINING IF ANCESTRAL!!
	double frac = (double)decvec.size()/(double)y->size();
	
	return lvar/frac;
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::istatistic()
{
    
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_istat;
	statresult p = onewayloop("DNA", "ALL", tc);
    stats["Istatistic"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::istatisticav()
{
    double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_istat;
	statresult p = onewayavloop("DNA", "ALL", tc);
    stats["IstatisticAv"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::modistatisticav()
{
	double (world::*tc) (vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_istat;
	statresult p = modonewayavloop("DNA", "ALL", tc);
    stats["ModIstatisticAv"] = p;
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_istat(vector<int>* y, int j)
{
	
	
	// variance of str on derived background
	
	vector<int> decvec  = getancdecvec(y, false);
	if(decvec.empty()) return -999;
	
	
	
	double phat = tc_nucdiv( y, j) ;
	
	//fraction of derived snps GOING WITH SEGS.SIZE BECAUSE MISSING DATA AT MICROSAT WILL STILL BE COUNTED AT the SNP FOR DETERMINING IF ANCESTRAL!!
	//double frac = ((double)dec.segs.size()) / ((double)a->segs.size());
	double frac = (double)decvec.size()/(double)y->size();
	
	return phat/frac;
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::delmusquared()
{
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_delmusqrd;
	statresult p = twowayloop("MICROSAT", "ALL", tc);
    stats["DeltaMuSquared"] = p;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::delmusquaredanc()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_delmusqrd;
	statresult ap = twowayloop("MICROSAT", "ANC", tc);
    stats["DeltaMuSquared-ANC"] = ap;
	
	
	
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::delmusquareddec()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_delmusqrd;
	
	
	statresult dp = twowayloop("MICROSAT", "DEC", tc);
    stats["DeltaMuSquared-DEC"] = dp;
	
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::ubdelmusquared()
{
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_ubdelmusqrd;
	statresult p = twowayloop("MICROSAT", "ALL", tc);
    stats["UBDeltaMuSquared"] = p;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::ubdelmusquaredanc()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_ubdelmusqrd;
	statresult ap = twowayloop("MICROSAT", "ANC", tc);
    stats["UBDeltaMuSquared-ANC"] = ap;
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::ubdelmusquareddec()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_ubdelmusqrd;
	
	
	statresult dp = twowayloop("MICROSAT", "DEC", tc);
    stats["UBDeltaMuSquared-DEC"] = dp;
	
	
}





//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_delmusqrd(vector<int>* a, vector<int>* b, int j)
{   
	
    return (double) pow((locmean(a, j) - locmean(b, j)), 2);
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_ubdelmusqrd(vector<int>* a, vector<int>* b, int j)
{
	return (double) (asd(a, b, j) - vx(a, j) - vx(b, j));
}



//---------------------------------------------
//Function Name:
//Function Definition: microsats only!
double world::asd(vector<int>* a, vector<int>* b, int j)
{
	double asd = 0.0;
	
	map<string, double>  aal = getdatafreqs(a, j);
	map<string, double>  bal = getdatafreqs(b, j);
	
	
	map<string, double> ::const_iterator ap = aal.begin();
	while(ap != aal.end()){
		map<string, double> ::const_iterator bp = bal.begin();
		while(bp != bal.end()){
			asd += pow(atof(ap->first.c_str()) - atof(bp->first.c_str()), 2) * (ap->second * bp->second);
			bp++;
		}
		ap++;
	}
	
	/*
	 vector<multimap<string, segment*, lociCmp > >::const_iterator loci = a->lmap.begin();
	 vector<multimap<string, segment*, lociCmp > >::const_iterator locj = b->lmap.begin();
	 
	 for(int i=0; i < loc; i++){ 
	 loci++;	
	 locj++;
	 }
	 
	 
	 multimap<string, segment*, lociCmp >::const_iterator locusi = loci->begin();
	 multimap<string, segment*, lociCmp >::const_iterator locusj;
	 double asd = 0.0;
	 while(locusi != loci->end()){ // iterate through all i alleles
	 locusj = locj->begin(); //set j allel to beginning of list
	 while(locusj != locj->end()){	//iterqatee through all j alleles for each i allele
	 //cout << "LENGTH: " << locusi->first << "  " << "Count: " << loci->count(locusi->first) << " out of " << loci->size() << " LENGTH: " << locusj->first << " Count: " << locj->count(locusj->first) << " out of " << locj->size() << endl;
	 //asd += ((double)(loci->count(locusi->first))/loci->size()) * ((double)(locj->count(locusj->first))/locj->size()) * abs(((atoi(locusi->first.c_str())) - (atoi(locusj->first.c_str())))); // Dw calc
	 asd += pow(atof(locusi->first.c_str()) - atof(locusj->first.c_str()), 2) * ((double)(loci->count(locusi->first))/loci->size()) * ((double)(locj->count(locusj->first))/locj->size());
	 
	 locusj = locj->upper_bound(locusj->first); //advance to next allele
	 }
	 locusi = loci->upper_bound(locusi->first); //advance to next allele
	 }
	 */
	//cout << "ASD: " << asd << endl;
	return asd;
}


//---------------------------------------------
//Function Name:
//Function Definition: microsats only!
double world::vx(vector<int>* y,int j)
{
	
	/*
	 vector<multimap<string, segment*, lociCmp > >::const_iterator loci = a->lmap.begin();
	 
	 for(int i=0; i < loc; i++){ 
	 loci++;	
	 }
	 
	 multimap<string, segment*, lociCmp >::const_iterator locusi = loci->begin();
	 */
	
	double vx = 0.0;
	
	double lmean = locmean(y,j);
	
	map<string, double> freqs = getdatafreqs(y, j);
	
	map<string, double> ::const_iterator p = freqs.begin();
	while(p != freqs.end()){
		vx += pow((atof(p->first.c_str()) - lmean), 2) * p->second;
		p++;
	}
	
	/*
	 while(locusi != loci->end()){ // iterate through all i alleles
	 vx += pow((atof(locusi->first.c_str()) - lmean), 2) * (loci->count(locusi->first))/loci->size();
	 locusi = loci->upper_bound(locusi->first); //advance to next allele
	 }
	 */
	
	double nover = (double)y->size()/(double)(y->size()-1);
	
	//double nover = (double)loci->size()/(double)(loci->size()-1);
	
	//cout << "pre-nover Vx " << vx << endl;
	
	//cout << "nover: " << nover << endl;
	
	//cout << "Vx: " << (nover * vx) << endl;
	return (nover*vx);
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::t_d()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_td;
	statresult p = twowayloop("MICROSAT", "ALL", tc);
    stats["Td"] = p;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::t_danc()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_td;
	statresult ap = twowayloop("MICROSAT", "ANC", tc);
    stats["Td-ANC"] = ap;
	
	
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::t_ddec()
{
	
	double (world::*tc) (vector<int>*, vector<int>*, int) = NULL; //fxn pointer with list of indivs and locus to use
	tc = &world::tc_td;
	
	
	statresult dp = twowayloop("MICROSAT", "DEC", tc);
    stats["Td-DEC"] = dp;
	
	
}





//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_td(vector<int>* a, vector<int>* b, int j)
{
	
    return (double) ((pow(locmean(a,j) - locmean(b,j), 2))/(2.0*(loci[j].mutrate)))  - (Vnaught/loci[j].mutrate);
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::NeiDist()
{
	if(dmap.size() > 1){
		double (world::*tc) (vector<int>*, vector<int>*, string) = NULL; //fxn pointer with list of indivs and locus to use
		tc = &world::tc_nei;
		statresult p = twowayxlociloop("ALL", "ALL", tc);
		stats["Nei"] = p;
	}
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_nei(vector<int>* a, vector<int>* b,  string locustype)
{
	double nei;
	double dxx = 0.0, dyy = 0.0, dxy = 0.0;
	double denom = 0.0;
	
	//nei 1987 p 219
	
	
	
	
	
	
	
	for(unsigned int k=0; k < loci.size(); k++){
		
		vector<int> y,z;
		vector<int>::iterator lp, lq;
		unsigned int vecsz = 0;
		lp = loci[k].lmap.begin();
		lq = loci[k].lmap.end();
		vecsz = loci[k].lmap.size();
		set_intersection(a->begin(), a->end(), lp, lq, std::back_inserter(y));
		set_intersection(b->begin(), b->end(), lp, lq, std::back_inserter(z));
		
		if((loci[k].type == locustype || locustype == "ALL") && (y.size() >= DEME_CUTOFF) && (z.size() >= DEME_CUTOFF)){
			
			
			
			
			map<string, double> afreqs = getdatafreqs(&y, k);
			map<string, double> bfreqs = getdatafreqs(&z, k);
			
			dxy += capd(&afreqs, &bfreqs);
			
			dxx += capd(&afreqs, &afreqs);
			dyy += capd(&bfreqs, &bfreqs);
			
			
			//cout << "dxy: " << dxy << "  dxx: " << dxx << "  dyy: " << dyy << endl;
			
			denom++;
			
			
			
		}
		
		
		//end k
	}
	//end if sysj	
	
	
	//end if testsystem
	
	
	
	dxy = dxy/denom;
	dxx = dxx/denom;
	dyy = dyy/denom;
	
	
	nei = dxy - ((dxx+dyy)/2);
	
	
	
    return nei;
}




//---------------------------------------------
//Function Name:
//Function Definition:
double world::capd(map<string, double>* a, map<string, double>* b)
{
    //for calculating Dxw, simply make loci and locj point tot he same place
    map<string, double>::const_iterator ap = a->begin();
    map<string, double>::const_iterator bp;
    double tot = 0.0;
	
    while(ap != a->end()){ // iterate through all i alleles
        bp = b->begin(); //set j allel to beginning of list
        while(bp != b->end()){	//iterqatee through all j alleles for each i allele
			
            if(ap->first != bp->first){ // for all i != j
				//something about the way frequencies are done
				//cout << ap->first << " " << bp->first << " " << (ap->second * bp->second) << endl;
				tot += (ap->second * bp->second); 
				
            }
			
			bp++;
        }
		
        ap++;
    }
    return tot;
}

//---------------------------------------------
//Function Name:
//Function Definition:
void world::DeeSW()
{
	if(dmap.size() > 1){
		double (world::*tc) (vector<int>*, vector<int>*, string) = NULL; //fxn pointer with list of indivs and locus to use
		tc = &world::tc_dsw;
		statresult p = twowayxlociloop("MICROSAT", "ALL", tc);
		stats["Dsw"] = p;
	}
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::tc_dsw(vector<int>* a, vector<int>* b,  string locustype)
{
	
	double totdsw = 0.0, dxw = 0.0, dyw = 0.0, dxyw = 0.0;
	double denom = 0.0;
	
	//nei 1987 p 219
	
	
	for(unsigned int k=0; k < loci.size(); k++){
		if((loci[k].type == locustype || locustype == "ALL") && (a->size() >= DEME_CUTOFF) && (b->size() >= DEME_CUTOFF)){
			map<string, double> afreqs = getdatafreqs(a, k);
			map<string, double> bfreqs = getdatafreqs(b, k);
			dxw = dw(&afreqs, &afreqs);
			dyw = dw(&bfreqs, &bfreqs);
			dxyw = dw(&afreqs, &bfreqs);
			
			
			
			denom++;
			
			
			
		}
		
		
		//end k
	}
	
	
	dxw = dxw/denom;
	dyw = dyw/denom;
	dxyw = dxyw/denom;
	
	
	totdsw = (double) dxyw - ((dxw+dyw)/2.0);
	
	
	
    return totdsw;
}



//---------------------------------------------
//Function Name:
//Function Definition:
double world::dw(map<string, double>* a, map<string, double>* b)
{
    //for calculating Dxw, simply make loci and locj point tot he same place
    map<string, double>::const_iterator ap = a->begin();
    map<string, double>::const_iterator bp;
    double tot = 0.0;
	
    while(ap != a->end()){ // iterate through all i alleles
        bp = b->begin(); //set j allel to beginning of list
        while(bp != b->end()){	//iterqatee through all j alleles for each i allele
			
            if(ap->first != bp->first){ // for all i != j
				tot += ((ap->second * bp->second) * abs(((atoi(ap->first.c_str())) - (atoi(bp->first.c_str()))))); // Dw calc
				
			}
			
			bp++;
        }
		
		ap++;
    }
    return tot;
}



//---------------------------------------------
//Function Name:
//Function Definition:
void world::multalleleLDp()
{
	if(loci.size() >1){
		double (world::*tc) (vector<int>*, int, int) = NULL; //fxn pointer with list of indivs and locus to use
		tc = &world::tc_multldp;
		statresult p = twolociloop("ALL", "ALL", tc);
		stats["LDChiSquarePvalue"] = p;
	}
	
}

//---------------------------------------------
//Class Name 
//Class Definition: Chi Square LD
double world::tc_multldp(vector<int>* y, int j, int k)
{
	
	//map<string, int> ycounts = getdatacounts(y,j);
	//map<string, int> zcounts = getdatacounts(z,k);
	//if(ycounts.size() >= 3 || zcounts.size() >= 3) return -999;
	
	//first check. if either one or the other locus has more than two alleles, just return an NA
	//if(a->allelethreshold(loc0, 3, true) && a->allelethreshold(loc1,3, true)) return -999; //the test is not applicable if there is more tah two alleles in both loci
	int numchiover = 0;
	
	int usediterations = LDP_ITERATIONS; 
	
	//vector<vector<string> > vecs = a->vecsfromlocinomiss(loc0, loc1);
	
	
	
	vector<int> locs;
	locs.push_back(j);
	locs.push_back(k);
	
	vector<int> nomiss = getnomissing(y, locs);
	
	vector<string> yal = getallelesasvec(&nomiss,j);
	vector<string> zal = getallelesasvec(&nomiss,k);
	
	
	
	
	//vector<string> vec0 = a->vecfromloci(loc0);
	//vector<string> vec1 = a->vecfromloci(loc1);
	//map<string, map <string, double> > counts = twowaycount(vec0, vec1);
	//map<string, map <string, double> > counts = twowaycount(vecs[0], vecs[1]);
	map<string, map <string, double> > counts = twowaycount(yal, zal);
	double X2 = chisqrld(counts);
	
	if(X2 == -999){// No comparisons possible
		return -999;
	}
	
	//randomize lost of first 
	for(int iti=0; iti < LDP_ITERATIONS; iti++){
		
		vector<string> vec0r = yal;
		vector<string> vec1r = zal;
		
		//randomize the r vec. Only necessary to randomize one to make the reltive positions random
		random_shuffle(vec0r.begin(), vec0r.end());
		map<string, map <string, double> > countsr = twowaycount(vec0r, vec1r);
		double X2r = chisqrld(countsr);
		if(X2r == -999){
			usediterations--;
		}
		
		else{
			//cout << "X2: " << X2 << " X2r: " << X2r << endl;
			if(X2r > X2) numchiover++;
		}
		//end iti	
	}
	
	//cout << endl << endl << endl;
	
	if(usediterations < 1){
		return -999;
	}
	
	return ((double)numchiover/(double)usediterations);
}

//---------------------------------------------
//Function Name:
//Function Definition: creates, from two vectors known to be linked position-by-position, a 2d table of counts indexed by their names
map<string, map <string, double> > world::twowaycount(vector<string> a, vector<string> b)
{
	map<string, map <string, double> > o;
	set<string> as;
	set<string> bs;
	
	
	//first test sizes, they must match
	if(a.size() != b.size()){
		cerr << "Error in test::twowaycount. Vector sizes do not match!";
		abort();
	}
	
	for (unsigned int i=0; i<a.size(); i++){
		as.insert(a[i]);
		bs.insert(b[i]);
		map <string, double> inside = o[a[i]];//should create one if emptry so far
		map <string, double>::const_iterator inp = inside.find(b[i]);
		if(inp == inside.end()){ //it was empty i.e. count is 0
			inside[b[i]] = 1;
		}
		else{//add 1
			//double cur = inside[b[i]];
			//cur++;
			//inside[b[i]] = cur;
			inside[b[i]]++;
		}
		o[a[i]] = inside;
	}
	
	//fill empty space in 2xN with zeroes
	set<string>::const_iterator asp = as.begin();
	while(asp != as.end()){
		//cout << *asp << endl;
		map<string, map <string, double> >::iterator outa = o.find(*asp);
		if(outa == o.end()){
			cerr << "Error in twowaycount. Expected at least one entry of " << *asp << "in one of the vectors.";
			abort();
		}
		set<string>::const_iterator bsp = bs.begin();
		while(bsp != bs.end()){
			map<string, double>::iterator inb = outa->second.find(*bsp);
			if(inb ==  outa->second.end()){
				//create a zero entry!
				outa->second[*bsp] = 0;
			}
			bsp++;
		}
		asp++;
	}
	
	
	
	
	
	
	
	return o;	
}

//---------------------------------------------
//Function Name:
//Function Definition: returns an LD X2 result from a 2d map of counts created by twolociu:twowaycount
//double twoloci::chisqrld(map<string, map<string, double> > a)
double world::chisqrld (map<string, map<string, double> > a)
{
	//first turn into a 2d map of frequencies by dividing through by the total
	//get total
	
	int numhaplo = 0;
	
	double tot = 0.0;
	map<string, map<string, double> >::const_iterator cpi = a.begin();
	while(cpi != a.end()){
		//cout << "cpi: " << cpi->first << endl;
		map<string, double>::const_iterator cpj = cpi->second.begin();
		while(cpj != cpi->second.end()){
			//cout << "   cpj: " << cpj->first << endl;
			tot += cpj->second;
			numhaplo++;
			cpj++;
		}
		cpi++;
	}
	
	
	//IF there is a single haplotype, there is NO variation in the entire deme
	if(numhaplo == 1){
		return -999;
	}
	
	map<string, double> toti;
	map<string, double> totj;
	
	map<string, map<string, double> >::iterator pi = a.begin();
	while(pi != a.end()){
		//cout << "pi: " << pi->first << endl;
		map<string, double>::iterator pj = pi->second.begin();
		while(pj != pi->second.end()){
			//cout << "   pj: " << pj->first << "    " << pj->second << endl;
			pj->second = pj->second/tot;
			//cout << "pj-second/tot:   " << pj->second << endl;
			
			//add to column totals outside loop first
			map <string, double>::const_iterator totip = toti.find(pi->first);
			if(totip == toti.end()){ //it was empty i.e. count is 0
				toti[pi->first] = pj->second;
			}
			else{//add current
				//double cur = toti[pi->first];
				//cur += pj->second;
				//toti[pi->first] = cur;
				toti[pi->first] += pj->second;
			}
			
			//then inside loop
			map <string, double>::const_iterator totjp = totj.find(pj->first);
			if(totjp == totj.end()){ //it was empty i.e. count is 0
				totj[pj->first] = pj->second;
			}
			else{//add current
				//double cur = totj[pj->first];
				//cur += pj->second;
				//totj[pj->first] = cur;
				totj[pj->first] += pj->second;
			}
			
			
			pj++;
		}
		pi++;
	}
	
	
	// now calc the X2
	double X2 = 0.0;
	
	
	cpi = a.begin();
	while(cpi != a.end()){
		map<string, double>::const_iterator cpj = cpi->second.begin();
		while(cpj != cpi->second.end()){
			double pipj = ((toti[cpi->first])*(totj[cpj->first]));
			double Dij = cpj->second - pipj;
			//cout << "DIJ: " << Dij << endl;
			X2 += (tot * pow(Dij,2.0))/pipj;
			cpj++;
		}
		cpi++;
	}
	
	
	return X2;
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::multalleleLD()
{
	if(loci.size() >1){
		double (world::*tc) (vector<int>*, int, int) = NULL; //fxn pointer with list of indivs and locus to use
		tc = &world::tc_multld;
		statresult p = twolociloop("ALL", "ALL", tc);
		stats["LDChiSquare"] = p;
	}
	
}


//---------------------------------------------
//Class Name:
//Class Definition: Chi Square LD
double world::tc_multld(vector<int>* y, int j, int k)
{
	//first check. if either one or the other locus has more than two alleles, just return an NA
	//if(a->allelethreshold(loc0, 3, true) && a->allelethreshold(loc1,3, true)) return -999; //the test is not applicable if there is more tah two alleles in both loci
	
	/*
	 map<string, int> ycounts = getdatacounts(y,j);
	 map<string, int>::const_iterator yy = ycounts.begin();
	 while(yy != ycounts.end()){
	 cout << yy->first << " " << yy->second << endl;
	 yy++;
	 }
	 
	 map<string, int> zcounts = getdatacounts(z,k);
	 map<string, int>::const_iterator zz = zcounts.begin();
	 while(zz != zcounts.end()){
	 cout << zz->first << " " << zz->second << endl;
	 zz++;
	 }
	 */
	
	//cout << endl << endl;
	
	//if(ycounts.size() >= 3 || zcounts.size() >= 3) return -999;
	
	vector<int> locs;
	locs.push_back(j);
	locs.push_back(k);
	
	vector<int> nomiss = getnomissing(y, locs);
	
	vector<string> yal = getallelesasvec(&nomiss,j);
	vector<string> zal = getallelesasvec(&nomiss,k);
	
	//vector<vector<string> > vecs = a->vecsfromlocinomiss(loc0, loc1);
	//vector<string> vec0 = a->vecfromloci(loc0);
	//vector<string> vec1 = a->vecfromloci(loc1);
	//map<string, map <string, double> > counts = twowaycount(vecs[0], vecs[1]);
	map<string, map <string, double> > counts = twowaycount(yal, zal);
	
	/*
	 map<string, map <string, double> >::const_iterator cp = counts.begin();
	 while(cp != counts.end()){
	 cout << "----- "<< cp-> first << endl;
	 map <string, double>::const_iterator cpp = cp->second.begin();
	 while(cpp != cp->second.end()){
	 cout << cpp->first << " " << cpp->second << endl;
	 cpp++;
	 }
	 
	 
	 
	 
	 cp++;
	 }
	 
	 */
	
	
	//map<string, map <string, double> > counts = twowaycount(vec0, vec1);
	double X2 = chisqrld(counts);
	
	
	return X2;
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::LD()
{
	
	if(loci.size() >1){
		double (world::*tc) (vector<int>*, int, int) = NULL; //fxn pointer with list of indivs and locus to use
		tc = &world::tc_ld;
		statresult p = twolociloop("ALL", "ALL", tc);
		stats["LD"] = p;
	}
	
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
void world::LDadj()
{
	
	if(loci.size() >1){
		double (world::*tc) (vector<int>*, int, int) = NULL; //fxn pointer with list of indivs and locus to use
		tc = &world::tc_ld;
		statresult p = twolociadjloop("ALL", "ALL", tc);
		stats["LDadj"] = p;
	}
	
	
}


//---------------------------------------------
//Class Name:
//Class Definition: LD
double world::tc_ld(vector<int>* y, int j, int k)
{
	//first check. if either one or the other locus has more than two alleles, just return an NA
	map<string, int> ycounts = getdatacounts(y,j);
	map<string, int> zcounts = getdatacounts(y,k);
    //cout << "ycounts: " << ycounts.size() << " zcounts: " << zcounts.size() << endl;
	if(ycounts.size() >= 3 || zcounts.size() >= 3 || ycounts.size() <= 0 || zcounts.size() <= 0) return -999;
	
	
	
	vector<vector<double> > cvec;
	cvec.push_back(vector<double> (2,0.0));
	cvec.push_back(vector<double> (2,0.0));
	
	
	
	
	vector<int> locs;
	locs.push_back(j);
	locs.push_back(k);
	
	vector<int> nomiss = getnomissing(y, locs);
	
	
	
	string yfirst = getdata(nomiss[0], j);
	string zfirst = getdata(nomiss[0], k);
	
	
	
	double tot = 0.0;
	
	//for (unsigned int i = 0; i < a->segs.size(); i++){
	for (unsigned int i = 0; i < nomiss.size(); i++){
		string thej = getdata(nomiss[i], j);
		string thek = getdata(nomiss[i], k);
		
		//cout << "nomiss i: " << nomiss[i]  << " locj:" << thej << " lock:" << thek << endl;
		if((thej != "NA") && (thek != "NA")){
			tot++;
			
			if(thej == yfirst){
				if(thek == zfirst) cvec[0][0]++;
				else cvec[0][1]++;
			}
			else{
				if(thek == zfirst) cvec[1][0]++;
				else cvec[1][1]++;
			}
			
		}
	}
	
	double x11 = cvec[0][0]/tot;
	double x12 = cvec[0][1]/tot;
	double x21 = cvec[1][0]/tot;
	double x22 = cvec[1][1]/tot;
	
	double p1 = x11+x12;
	double p2 = x21+x22;
	double q1 = x11+x21;
	double q2 = x12+x22;
	
	//Gimme a D
	double D = (cvec[0][0] / tot) - (p1 * q1);
	
	double r2 = (D*D)/(p1*p2*q1*q2);
	
	//Normalize
	
	double r2c = (r2-(1/tot))/(1-(1/tot));
	
	return r2c;
	
	
}




//---------------------------------------------
//Function Name:
//Function Definition: Return true if the only entry in the list is "ALL"
bool priordist::testlistall(){
	map<string,double>::const_iterator p = alphalist.begin();
	
	if(p->first =="ALL" && alphalist.size() == 1) return true;
	
	return false;
}



//---------------------------------------------
//Function Name:
//Function Definition: Perform a limited list of tests
void world::dotests(priordist* pt, bool timetests){
	
	//time_t teststart, pretest, posttest;
	//time (&teststart);
	
	//time setup
	timeval teststart, pretest, posttest;
	//Stamp time
	gettimeofday(&teststart, NULL);
	double starttime = teststart.tv_sec + (teststart.tv_usec/1000000.0);
	
	map<string,double>::const_iterator p = pt->alphalist.begin();
	
	while(p != pt->alphalist.end()){
		
		//time(&pretest);
		
		gettimeofday(&pretest, NULL);
		double pretime = pretest.tv_sec + (pretest.tv_usec/1000000.0);
		
		//cout << "Test:" << p->first << endl;
		
		if(p->first == "Heterozygosity") heterozygosity(); 
		//else if(p->first == "HeterozygositySub") heterozygositysub();
		else if(p->first == "Heterozygosity-ANC") heterozygosityanc();
		else if(p->first == "Heterozygosity-DEC") heterozygositydec();
		else if(p->first == "HeterozygosityAv") heterozygosityav();	
		//else if(p->first == "HeterozygosityAvSub") heterozygosityavsub();
		else if(p->first == "HeterozygosityAv-ANC") heterozygosityavanc();
		else if(p->first == "HeterozygosityAv-DEC") heterozygosityavdec();
		else if(p->first == "DeltaMuSquared") delmusquared();
		//else if(p->first == "DeltaMuSquaredSub") delmusquaredsub();
		else if(p->first == "DeltaMuSquared-ANC") delmusquaredanc();
		else if(p->first == "DeltaMuSquared-DEC") delmusquareddec();
		else if(p->first == "Td") t_d();
		//else if(p->first == "TdSub") t_dsub();
		else if(p->first == "Td-ANC") t_danc();
		else if(p->first == "Td-DEC") t_ddec();
		else if(p->first == "Dsw") DeeSW();
		else if(p->first == "Nei") NeiDist();
		else if(p->first == "Beta") betaimbalance();
		else if(p->first == "BetaAv") betaimbalanceav();
		//else if(p->first == "BetaSub") betaimbalancesub();
		else if(p->first == "Beta-ANC") betaimbalanceanc();
		else if(p->first == "Beta-DEC") betaimbalancedec();
		//else if(p->first == "BetaAvSub") betaimbalanceavsub();
		else if(p->first == "BetaAv-ANC") betaimbalanceavanc();
		else if(p->first == "BetaAv-DEC") betaimbalanceavdec();
		else if(p->first == "LD") LD();
		else if(p->first == "LDadj") LDadj();
		else if(p->first == "LDChiSquare") multalleleLD();
		else if(p->first == "LDChiSquarePvalue") multalleleLDp();
		else if(p->first == "RangeofLocus") rangeoflocus();
		else if(p->first == "RangeofLocusAv") rangeoflocusav();
		else if(p->first == "SampleSize") samplesize();
		//else if(p->first == "SampleSizeSub") samplesizesub();
		else if(p->first == "SampleSize-ANC") samplesizeanc();
		else if(p->first == "SampleSize-DEC") samplesizedec();
		else if(p->first == "MaxLength") maxlength();
		else if(p->first == "MinLength") minlength();
		else if(p->first == "MeanLength") meanlength();
		else if(p->first == "NumDiffAlleles") numdiffalleles();
		else if(p->first == "NumDiffAllelesAv") numdiffallelesav();
		else if(p->first == "CompleteHeterozygotes") allheterozygotes();
		else if(p->first == "CompleteHomozygotes") allhomozygotes();
		else if(p->first == "Haplotypes") haplotypes();
		else if(p->first == "Jstatistic") jstatistic();
		else if(p->first == "JstatisticAv") jstatisticav();
		else if(p->first == "ModJstatisticAv") modjstatisticav();
		else if(p->first == "Istatistic") istatistic();
		else if(p->first == "IstatisticAv") istatisticav();
		else if(p->first == "ModIstatisticAv") modistatisticav();
		else if(p->first == "MicrosatVariance") microsatvariance();
		else if(p->first == "MicrosatVarianceAv") microsatvarianceav();
		else if(p->first == "DerivedFraction") derivedfraction();
		else if(p->first == "NucleotideDiversity") nucleotidediversity();
		else if(p->first == "NucleotideDiversityAv") nucleotidediversityav();
		else if(p->first == "Fst") Fstdist();
		else if(p->first == "FstAv") Fstdistav();
		else if(p->first == "UBDeltaMuSquared") ubdelmusquared();
		//	else if(p->first == "UBDeltaMuSquaredSub") ubdelmusquaredsub();
		else if(p->first == "UBDeltaMuSquared-ANC") ubdelmusquaredanc();
		else if(p->first == "UBDeltaMuSquared-DEC") ubdelmusquareddec();
		else if(p->first == "Hamming") hamming();
		else if(p->first == "HammingBtwn") hammingbtwn();
		else if(p->first == "NumHaploBlocks") numhaploblocks();
		else if(p->first == "BlockLength") blocklength();
		else if(p->first == "BlockBoundary") 	blockboundary();
		else if(p->first == "NumROH") 	numroh();
	//	else if(p->first == "ROHCorr") 	rohcorr();
		else if(p->first == "FROH") 	froh();
		else if(p->first == "FROHVar") frohvar();
		else if(p->first == "SharedHaplotypeBlocks") sharedhaploblocks();
		else if(p->first == "MinorAlleleFrequency") minorallelefreq();
		else if(p->first == "TajimasD") tajimasd();
		else if(p->first == "NumSegSites") numsegsites();
		else if(p->first == "LongSharedBetween") longsharedbetween();
		else if(p->first == "LongSharedWithin") longsharedwithin();
        else if(p->first == "VarROHLength") varrohlength();

		
		
		else {
			cerr << "UNRECOGNIZED TEST CALL " << p->first << endl;
			abort();
		}
		
		//time (&posttest);
		//double postteststart = difftime(posttest, teststart);
		//double posttestdiff = difftime(posttest, pretest);
		
		gettimeofday(&posttest, NULL);
		double posttime = posttest.tv_sec + (posttest.tv_usec/1000000.0);
		
		if(timetests) cout << "Time elapsed: " << posttime-starttime << "  Test " << p->first << " took " << posttime-pretime << " seconds." << endl;
		
		p++;
	}
	
	
	
	
}




//---------------------------------------------
//Function Name:
//Function Definition: Perform all tests
void world::doalltests(){
	heterozygosity();
	heterozygosityav();
	heterozygosityanc();
	heterozygositydec();
	heterozygosityavanc();
	heterozygosityavdec();
	delmusquared();
	delmusquaredanc();
	delmusquareddec();
	t_d();
	t_danc();
	t_ddec();
	DeeSW();
	NeiDist();
	nucleotidediversity();
	betaimbalance();
	betaimbalanceav();
	betaimbalanceanc();
	betaimbalancedec();
	betaimbalanceavanc();
	betaimbalanceavdec();
	LD();
	LDadj();
	multalleleLD();
	multalleleLDp();
	rangeoflocus();
	rangeoflocusav();
	maxlength();
	minlength();
	meanlength();
	samplesize();
	samplesizeanc();
	samplesizedec();
	numdiffalleles();
	numdiffallelesav();
	allhomozygotes();
	allheterozygotes();
	haplotypes();
	jstatistic();
	jstatisticav();
	modjstatisticav();
	istatistic();
	istatisticav();
	modistatisticav();
	microsatvariance();
	microsatvarianceav();
	derivedfraction();
	Fstdist();
	Fstdistav();
	nucleotidediversity();
	nucleotidediversityav();
	ubdelmusquared();
	ubdelmusquaredanc();
	ubdelmusquareddec();
	hamming();
	numhaploblocks();
	blocklength();
	blockboundary();
	numroh();
	//rohcorr();
	froh();
	frohvar();
	hammingbtwn();
	sharedhaploblocks();
	minorallelefreq();
	tajimasd();
	numsegsites();
	longsharedbetween();
	longsharedwithin();
    varrohlength();

	
}



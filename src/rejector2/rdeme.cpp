//------------------------------------------------
//
// REJSTATS
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropology
// Stanford University, 2007
//
//------------------------------------------------

#include <cstring>
#include <sstream>
#include <set>
#include "rdeme.h"
#include "public.h"



int world::locinfosizechk(ostream &out){
	int totsize =0;	
	int i = 0;
	vector<locinfo>::iterator lp = loci.begin();
	while(lp != loci.end()){
		out << i;
		//int lsc = lp->locinfosizechk(out);
		
		
		//
		
		int locsize =0;	
		
		locsize += sizeof(lp->type);
		locsize += sizeof(lp->mutrate);
		
		out << " Type: " << sizeof(lp->type) << " Mut: " << sizeof(lp->mutrate);
		
		int anc = 0;
		vector<string>:: const_iterator ancp = lp->ancstate.begin();
		while(ancp != lp->ancstate.end()){
			anc += sizeof(*ancp);
			ancp++;
		}
		
		locsize += anc;
		locsize += sizeof(lp->length);
		locsize += sizeof(lp->offset);
		locsize += sizeof(lp->recrate);
		
		out << " Anc: " << anc << " Length: " << sizeof(lp->length) << " Offset: " << sizeof(lp->offset) << " Rec: " << sizeof(lp->recrate);
		
		int lm = 0;
		vector<int>::const_iterator lmp = lp->lmap.begin();
		while(lmp != lp->lmap.end()){
			lm += sizeof(*lmp);
			lmp++;
		}
		out << " Lmap: " << lm;
		
		locsize += lm;
		
		
		//
		
		totsize += locsize;
		
		
		//totsize += lsc;
		out << " Tot: " << locsize;
		out << " Tot from Outside: " << sizeof(*lp) << endl;
		i++;
		lp++;
	}
	
	return totsize;
}


int world::worldsizechk(ostream &out){
	int sizetot = 0;
	
	int locisize = locinfosizechk(out);
	sizetot += locisize;
	
	int datasize =0;
	for(int i = 0 ; i < datalength ; i++ ){
		for(int j=0; j<datawidth; j++){
			datasize += sizeof(char);
		}
	}
	sizetot += datasize;
	
	
	int nmapsize = 0;
	vector<string>:: const_iterator nmapp = nmap.begin();
	while(nmapp != nmap.end()){
		nmapsize += sizeof(*nmapp);
		nmapp++;
	}
	
	sizetot += nmapsize;
	
	
	int dmapsize = 0;
	map<string, vector<int> >::const_iterator dp  = dmap.begin();
	while(dp != dmap.end()){
		vector<int>::const_iterator dpp = dp->second.begin();
		while(dpp != dp->second.end()){
			dmapsize += sizeof(*dpp);
			
			dpp++;
		}
		dp++;
		
		
		
	}
	
	sizetot += dmapsize;
	
	
	return sizetot;
}

//---------------------------------------------
//Function Name:
//Function Definition: returns true if all parts of a locus match for two entries. Can turn on check for missing data too.
bool world::datamatch(int x, int y, int i, bool missingchk){
	for(int j=i; j<(i+loci[i].length); j++){
		
		if(missingchk){
			if(data[x][j] == 'X' || data[y][j] == 'X') {
				return false;
			}
		}
		
		if(data[x][j] !=  data[y][j]){ 
			return false;	
		}
		
	}
	return true;
	
}


//---------------------------------------------
//Function Name:
//Function Definition: tests if any of a string is a digit
bool world::anydigits(const char* Str) const
{
    const int len = strlen(Str);
    for(int i = 0;i < len;i++) if (isdigit(Str[i])) return true;
    return false;
}

//---------------------------------------------
//Function Name:
//Function Definition:
string world::getdata(int i, int lp){
	string entry;
	for(unsigned int j = loci[lp].offset; j < (loci[lp].offset+loci[lp].length); j++){
		
		entry += data[i][j];
		
	}
	return entry;
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::getdataasdbl(int i, int lp){
	string entry;
	for(unsigned int j = loci[lp].offset; j < (loci[lp].offset+loci[lp].length); j++){
		
		entry += data[i][j];
		
	}
	return atof(entry.c_str());
}


//---------------------------------------------
//Function Name:
//Function Definition: This means no duplicates! Just a list of all the alleles
set<string> world::getalleles(vector<int>* y, int j){
	
	set<string> alleles;
	
	vector<int>::const_iterator p = y->begin();
	
	while(p != y->end()){
		
		string entry = getdata(*p, j);
		
		string::const_iterator ep = entry.begin();
		
		if(*ep  != 'X') alleles.insert(entry);
		
		p++;
	}
	
	
	return alleles;	
}

//---------------------------------------------
//Function Name:
//Function Definition:
vector<string> world::getallelesasvec(vector<int>* y, int j){
	
	vector<string> alleles;
	
	vector<int>::const_iterator p = y->begin();
	
	while(p != y->end()){
		
		string entry = getdata(*p, j);
		
		string::const_iterator ep = entry.begin();
		
		if(*ep  != 'X') alleles.push_back(entry);
		
		p++;
	}
	
	
	return alleles;	
}


//---------------------------------------------
//Function Name:
//Function Definition:
vector<int> world::getintersection(vector<int>* y, vector<int>* z){
	
	vector<int> a;
	
	set_intersection(y->begin(), y->end(), z->begin(), z->end(), std::back_inserter(a));
	
	return a;
	
}

//---------------------------------------------
//Function Name:
//Function Definition:
vector<int> world::getsubset(vector<int>* y, int j, string match){
	
	vector<int> sub;
	
	vector<int>::const_iterator p = y->begin();
	
	
	while(p != y->end()){
		
		string entry = getdata(*p, j);
		string::const_iterator ep = entry.begin();
		
		
		if(*ep  != 'X' && entry == match) sub.push_back(*p); 
		
		
		
		
		p++;
	}
	
	
	return sub;	
}


//---------------------------------------------
//Function Name:
//Function Definition:
bool world::missingchk(int match, vector <int> locs, int j){
	
	vector<int>::iterator location;
	location = std::find( loci[locs[j]].lmap.begin(), loci[locs[j]].lmap.end(), match );
	
	if ( location != loci[locs[j]].lmap.end() ) {
		if((j+1) == locs.size()) return true;
		else{
			if(missingchk(match, locs, (j+1))) return true;
			else return false;
		}
	}
	
	
	
	
	return false;
}

//---------------------------------------------
//Function Name:
//Function Definition:
vector<int> world::getnomissing(vector<int>* y, vector <int> locs){
	
	vector<int> nomiss;
	
	vector<int>::const_iterator p = y->begin();
	
	
	while(p != y->end()){
		
		if(missingchk(*p, locs, 0)) nomiss.push_back(*p);
		
		p++;
	}
	
	
	return nomiss;	
}




//---------------------------------------------
//Function Name:
//Function Definition:
map<string, int> world::getdatacounts(vector<int>* y, int j){
	
	map<string, int> counts;
	
	vector<int>::const_iterator p = y->begin();
	
	while(p != y->end()){
		
		string entry = getdata(*p, j);
		
		string::const_iterator ep = entry.begin();
		//cout << entry << endl;
		
		if(*ep  != 'X'){
			
			if(counts.find(entry) == counts.end()){//if it does not have an entry make one
				counts[entry] = 1;
				
			}
			else{ //increment
				counts.find(entry)->second++;	
			}
		}
		
		
		
		
		p++;
	}
	
	
	return counts;	
}



//---------------------------------------------
//Function Name:
//Function Definition:
map<string, double> world::getdatafreqs(vector<int>* y, int j){
	
	map<string, double> freqs;
	
	map<string, int> counts = getdatacounts(y, j);
	
	map<string, int>::const_iterator p = counts.begin();
	
	while(p != counts.end()){
		//cout << p->first << " Freq: " << (double) p->second/(double) y->size() << " Count: " << p->second << " Size: " << y->size()	<< endl;
		freqs[p->first] = ((double) p->second/(double) y->size());
		p++;
	}
	
	
	return freqs;	
}


//---------------------------------------------
//Function Name:
//Function Definition: calculates total at a specified locus. use only on numerical sites like microsats
double world::loctot(vector<int>* y, int j){
    double tot = 0.0;
	
	
	vector<int>::const_iterator p = y->begin();
	
	while(p != y->end()){
		tot += getdataasdbl(*p, j);
		p++;
	}
	
    return tot;
}

//---------------------------------------------
//Function Name:
//Function Definition: computes the variance at a (numerical) vector
double world::locvar(vector<int>* y, int j){
	
	double sumsq = 0.0;
	double lmean = locmean(y,j);
	
	
	vector<int>::const_iterator yp = y->begin();
	while(yp != y->end()){
		double tval = atof(getdata(*yp, j).c_str());
		sumsq += pow((tval - lmean),2.0);
		yp++;
	}
	
	
	return(sumsq/(y->size()-1));
}














//---------------------------------------------
//Function Name:
//Function Definition: ensures all parts of a string are integers
bool StringIntTest(const char* Str)
{
    const int len = strlen(Str);
    for(int i = 0;i < len;i++) if (!isdigit(Str[i])) return false;
    return true;
}










//---------------------------------------------
//Function Name:
//Function Definition: simple boolean to determine if a world has no demes
bool world::isempty()
{
	//if the dmee map is empty, the world is empty
	if(dmap.empty()) return true;
	return false;
}



//---------------------------------------------
//Function Name:
//Function Definition: constructor creates empty world object
world::world(){
	
}

//---------------------------------------------
//Function Name:
//Function Definition: destructor.. calls function to delete all data
world::~world()
{
	
	
	//DELETE DATA
	for(int i = 0 ; i < datalength ; i++ ) delete[] data[i];
	delete[] data;
	
	/*
	 vector<locinfo>::iterator dell =  loci.begin();
	 while(dell != loci.end())
	 {
	 
	 delete *dell;
	 dell++;
	 }
	 
	 */
	
	
	
	/*
	 //DELETE LOCI
	 for(int j=0; j < loci.size(); j++)
	 {
	 
	 delete loci[j];
	 }
	 */
	
	
	
}






//---------------------------------------------
//Function Name:
//Function Definition: ensures all parts of a string are integers
bool lociCmp::strIntTest(const char* Str) const
{
    const int len = strlen(Str);
    for(int i = 0;i < len;i++) if (!isdigit(Str[i])) return false;
    return true;
}


//---------------------------------------------
//Function Name:
//Function Definition: checks to ensure that the number of alleles is at or over a threshold number if over is true, or at or under if false. Used mainly to ensure that calculations are not doen on demes with no SNP diversity.
bool world::allelethreshold(int loc, int thresh, bool over){ //
	
	set<string> alleles;
	
	vector<int>::const_iterator lp = loci[loc].lmap.begin();
	while(lp != loci[loc].lmap.end()){
		string thisdata = getdata(*lp , loc);
		alleles.insert(thisdata);
		if(over){
			if(alleles.size() > thresh) return true;
		}
		
		lp++;	
	}
	
	if(!over){
		if(alleles.size() < thresh) return true;
	}
	
	return false;
	
	
	
}




//---------------------------------------------
//Function Name:
//Function Definition: determines whether any in a set of alleles matches the ancestral state
bool world::allelematch(set<string> alleles, string match){
	set<string>::const_iterator ap = alleles.begin();
	while(ap != alleles.end()){
		if(*ap == match) return true;
		ap++;
	}
	return false;
}


//---------------------------------------------
//Function Name:
//Function Definition: calls allelethreshold for all systems
bool world::checksnpdiversity(int thresh, bool over){
	
	for(int i=0; i<loci.size(); i++){
		if(loci[i].type == "SNP"){
			if(over){
				if(allelethreshold(i, thresh, true)){
					cerr << "Error: Locus" << i << " has  more than " << thresh << " SNP alleles." << endl;
					return false;
				}
			}
			else{
				if(allelethreshold(i, thresh, false)){
					cerr << "Error: Locus" << i << " has  more than " << thresh << " SNP alleles." << endl;
					return false;
				}
				
				
			}
			
			
		}
		
	}
	
	
	
	
	
	
	
	return true;
}







//---------------------------------------------
//Function Name:
//Function Definition: check to make sure no name is present more than twice
bool world::allelenumcheck(){
	
	vector<string>::const_iterator np = nmap.begin(); //Referenced by the vec #
	map<string, int> namecounts;
	while(np != nmap.end()){
		if(namecounts.find(*np) == namecounts.end()){
			namecounts[*np] = 1;
		}
		else if(namecounts[*np] > 1) return false;//and thereforew adding 1 will make it greater than 2
		else namecounts[*np]++;
		
		
		
		np++;
	}
	return true;
	
	
}

//---------------------------------------------
//Function Name:
//Function Definition: controller function to simply outpuit all test data
void world::printalltests(ostream &out, string inFile){
	out << "Rejstats Output" << endl;
	out << "Infile: " << inFile << endl;
	
	time_t rawtime = time(0);
	struct tm * timeinfo = localtime ( &rawtime );
	
	out << "Printed " << asctime (timeinfo) << endl << endl;
	
	
	
	
	
	
	out << "Vnaught\t" << Vnaught << endl;
	
	
	out << endl << endl;
	
	
	map<string, statresult>::iterator t = stats.begin();
	
	
	
	
	
	
	
	while(t != stats.end()){
		out << "/---" << endl;
		out << t->first << endl;
		//cout << t->first << " ";
		statresult * sr = &t->second;
		printstatresult(sr, out);
		
		//t->second->printtest(out);
		out << endl << endl;
		
		t++;
	}
	
	
	out << "/---" << endl;
	out << "END TESTS" << endl << "Run complete. ";
}


//---------------------------------------------
//Function Name:
//Function Definition:
int world::ancdecchk(int i) {
	
	bool allna = true;
	
	
	int k = 0;
	for(unsigned int lp=0; lp < loci.size(); lp++){
		
		vector<string>::const_iterator ac = loci[lp].ancstate.begin();
		while(ac != loci[lp].ancstate.begin()){
			
			bool anc = false;
			
			if(*ac != "NA"){
				
				if(loci[lp].ancstate.size() > 1){
					cerr << "Error in ancdecchk function. Locus " << k << " contains an NA for ancestral info, but there is more than one entryu for this locus. This is not allowed.";
					abort();
				}
				
				anc = true; //NA's will not prevent matches unless all are loci NA
				
			}
			
			else{ //not an NA, so test
				
				allna = false;
				
				string entry;
				for(unsigned int j = loci[lp].offset; j < (loci[lp].offset+loci[lp].length); j++) entry += data[i][j];
				
				if(*ac == entry) {
					anc = true;
					break;
				}
				
			}
			
			
			if(!anc) return -1; //if any locus has no anc match, return that its dec
			
			
			
			
			ac++;
		}
		
		k++;
		lp++;
		
		//end step through loci
	}
	
	if(allna) return 0; //return that its neither
	return 1; //return that its anc
}


//---------------------------------------------
//Function Name:
//Function Definition:
vector<int> world::getancdecvec(vector<int>* orig, bool anc){
	vector<int> outvec;
	
	
	vector<int>::const_iterator p = orig->begin();
	
	while(p != orig->end()){
		int adchk = ancdecchk(*p);
		if(anc){
			if(adchk == 1){
				outvec.push_back(*p);
			}
			
		}
		else{//it's dec we're loking for then
			if(adchk == -1){
				outvec.push_back(*p);
			}
		}
		
		p++;
	}
	
	return outvec;
}










//---------------------------------------------
//Function Name: dv8
//Function Definition: returns specified deviate
double dv8(double v1, double v2, string dtype, bool posit, bool zeroup)
{
	double dv8val = -9999999.99;
	
	while(true){
		
		
		long lidum = 1L;
		if(dtype == "uniform") dv8val = ((v2-v1)*ran3(&lidum)+v1); // v1 min, v2 max
		else if (dtype == "gaussian") dv8val =((v2*gasdev(&lidum))+v1); //v1 mean, v2 standard deviation
		else if (dtype == "gamma") dv8val = (gamma_dev(v2)*v1); //v1 mean, v2 shape param
		else {
			cerr << "UNKNOWN DEVIATE TYPE";
			abort();
		}
		
		
		
		
		if(posit){
			if(dv8val > 0.0) return dv8val; //if it must be positive
		}
		else if(zeroup){
			if(dv8val >= 0.0) return dv8val; // if it must be zero or positive
		}
		else return dv8val; //otherwise just return
		
	}
	
	
}







//---------------------------------------------
//Function Name:
//Function Definition: calculates r-squared from input data
map<int, int> world::rsquaredblocks(vector<int> invec){
	
	
	map<int, int> localblocks;
	
	
	//sorted vector of LD's for input vector (usually the individuals in a deme)
	multimap<double, int> sorted;
	
	
	for(unsigned int i =0; i<(loci.size()-1); i++){
		
		double rsq =  tc_ld(&invec, i, i+1);
		
		//cout << i << ": " << rsq << endl;
		loci[i].rsquared = rsq;
		
		sorted.insert(make_pair(rsq, i));
		
		
		
		
		//end for i to loci
	}
	
	
	//The last one will have no right-side LD
	loci[loci.size()-1].rsquared = -999;
	//cout << loci.size()-1 << ": " << loci[loci.size()-1].rsquared << endl;
	
	
	multimap<double, int>::const_iterator sp = sorted.end();
	sp--;
	while(sp != sorted.begin()){
		
		//	cout << sp->first << ": " << sp->second << endl;
		
		if(sp->first > 0.4 && !isnan(sp->first)){
			
			bool inablock = false;
			map<int, int>::iterator bp = localblocks.begin(); //begin at first place its less than or equal to sp second
			while(bp != localblocks.end()){
				if(sp->second >= bp->first && sp->second <= bp->second){ //it is inside the block, so no new block need be created
					inablock = true;
					//end if nin a block
				}
				
				bp++;
			}
			
			int lower = sp->second;
			int upper = sp->second+1; //Since its rsq to the right, the next block is automatically in
			
			
			
			
			if(!inablock){
				
				//map<double, int> rsqset;
				
				
				
				
				//Upper
				//cout << "up: ";
				for(unsigned int up = upper; up <(loci.size()-1); up++){
					//cout << loci[up].rsquared << " ";
					
					//check that pairwises with new one are high enough - only need to check with the newest as by now all other combos were already checked
					bool blockpw = true;
					for(int pwchk = lower; pwchk <= upper; pwchk++){
						double pwrsq =  tc_ld(&loci[pwchk].lmap, pwchk, up);
						if(pwrsq < 0.1){
							blockpw = false;
							break;
						}
					}
					
					
					if(loci[up].rsquared > 0.3 && blockpw && !isnan(loci[up].rsquared)){//if there's enough LD, extend the block
						upper++;
						
						//end up >0.3	
					}
					else break;
					//end for up
				}
				
				//cout << endl;
				//Lower
				//cout << "down: ";
				if(lower >0){
					for(unsigned int low = (lower-1); low >0; low--){
						//cout << loci[low].rsquared << " ";
						
						
						//check that pairwises with neww one are high enough - only need to check with the newest as by now all other combos were already checked
						bool blockpw = true;
						for(int pwchk = lower; pwchk <= upper; pwchk++){
							double pwrsq =  tc_ld(&loci[pwchk].lmap, pwchk, low);
							if(pwrsq < 0.1){
								blockpw = false;
								break;
							}
						}
						
						if(loci[low].rsquared > 0.3 && blockpw && !isnan(loci[low].rsquared)){//if there's enough LD, extend the block
							lower--;
							
							//end up >0.3	
						}
						else break;
						//end for up
					}
				}
				//cout << endl;
				
				//end if !inablock
				
			}
			
			
			//cout << "New block between " << lower << " and " << upper << endl;
			localblocks.insert(make_pair(lower, upper));
			
			
			
			//end if sp->first ios less than alpha
		}
		
		
		sp--;
	}		
	
	
	//Merge overlapping blocks
	map<int, int> finalblocks;
	
	
	
	if (localblocks.size()>0) {
		
		map<int, int>::iterator bx = localblocks.begin(); //begin at first place its less than or equal to sp second
		
		int low = bx->first;
		int high = bx->second;
		//cout << endl << endl << "LOCAL" << endl;
		while(bx != localblocks.end()){
			
			
			
			//cout << bx->first << " : " << bx->second << endl;
			if (bx->first <= high) {//if the low end touches or is lower than the current block high end, incorporate the block
				high = bx->second; //keep growing block
			}
			else{ //new block to finalblocks
				finalblocks.insert(make_pair(low, high));
				low = bx->first;
				high = bx->second;
				
				
			}
			
			bx++;
		}
		
		//insert the last one
		finalblocks.insert(make_pair(low, high));
		//end isempty localblocks
	}
	
	/*
	 cout << endl << endl << "FINAL" << endl;
	 map<int, int>::iterator fx = finalblocks.begin(); //begin at first place its less than or equal to sp second
	 while(fx != finalblocks.end()){
	 cout << fx->first << " : " << fx->second << endl;
	 
	 fx++;
	 }
	 */
	
	
	
	return finalblocks;
}


//---------------------------------------------
//Function Name:
//Function Definition: calculates r-squared from input data
void world::rsquaredblocksbydeme(){
	
	
	map<string, vector<int> >::iterator dp = dmap.begin();
	while (dp != dmap.end()) {
		//cout << "*********************************" << endl << endl << dp->first << endl << endl;
		map<int, int> curblock = rsquaredblocks(dp->second);
		dblocks[dp->first] = curblock;
		
		
		/*	
		 map<int, int>::iterator bp = curblock.begin();
		 while(bp != curblock.end()){
		 cout << bp->first << " : "<< bp->second << endl;
		 bp++;
		 }
		 
		 cout << "*********************************" << endl << endl;
		 */
		
		dp++;
	}
	
	
	
	//cout << "Size of dblocks: " << dblocks.size() << endl;
	
}


//---------------------------------------------
//Function Name:
//Function Definition: finds runs of homozygosity for the given samples
vector<bool> world::runsofhomozyg(vector<int>* y){
	
	unsigned int i;
	
	
	
	vector<double> runbysnp;
	
	for(i = 0; i<loci.size(); i++){//quick and dirty
		runbysnp.push_back(0.0);
	}
	
	
	
	double numindivs = 0.0;
	
	//go through and grab all double entries. i.e. two chromosomes in segs listing with same name
	vector<int>::const_iterator thisun = y->begin();
	while(thisun != --y->end()){ //end at one before acvtual end as nextun will be the last one
		bool homstate = false;
		vector<int>::const_iterator nextun = thisun;
		nextun++;
		
		
		//cout << "thisun:" << nmap[*thisun] << " nextun:" << nmap[*nextun];
		
		//int numrunchk = 0;
		//int pdchk = 0;
		
		
		if(nmap[*thisun] == nmap[*nextun]){//they have same names thus are two chromosomes from same individual
			vector<bool> homvec;
			numindivs++;
			
			
			for(unsigned int i = 0; i<loci.size(); i++){
				homstate = true;
				string thisdata = getdata(*thisun, i);
				if(thisdata[0] == 'X') homstate = false;//assume that missing data anywhere will invalidate test and assume can't SHOW its allhet 2/9/06
				string nextdata = getdata(*nextun, i);
				if(nextdata[0] == 'X') homstate = false;//assume that missing data anywhere will invalidate test and assume can't SHOW its allhet 2/9/06
				if(thisdata != nextdata) homstate = false; 
				homvec.push_back(homstate);
				
			}
			
			//Now examine for runs and ONLY insert as true if it is over a user-defined minimum run sized..
			vector<bool> runvec;
			
			int truecount = 0;
			int chk = 0;
			
			unsigned int j = 0;
			int physrun = 0;
			vector<bool>::const_iterator hp = homvec.begin();
			while(hp != homvec.end()){
				if(*hp == false){
					
					if(truecount >= MIN_HOMOZYGOSITY_RUN && physrun >= MIN_HOMOZYGOSITY_LENGTH){
						//numrunchk += truecount;
						//pdchk += physrun;
						//cout << "Run of " << truecount << " SNPs and " << physrun << " distance" << endl; 
						for(int i =0; i < truecount; i++){
							runvec.push_back(true);	
							chk++;
						}
						
					}
					else{
						for(int i =0; i < truecount; i++){
							runvec.push_back(false);
							chk++;
						}
						
					}
					
					truecount = 0;
					physrun = 0;
					//add on the false at the end
					runvec.push_back(false);
					chk++;
					
				}
				
				else {
					truecount++;
					physrun += loci[j].phdist;
				}
				j++;
				hp++;	
			}
			
			
			
			
			//If it does not end on a false, read out the final trues
			if(truecount >= MIN_HOMOZYGOSITY_RUN && physrun >= MIN_HOMOZYGOSITY_LENGTH){
				//numrunchk += truecount;
				//pdchk += physrun;
				//cout << "Run of " << truecount << " SNPs and " << physrun << " distance" << endl; 
				for(int i =0; i < truecount; i++){
					runvec.push_back(true);	
					chk++;
				}
				
				
				
			}
			else{
				for(int i =0; i < truecount; i++){
					runvec.push_back(false);	
					chk++;
				}
				
				
			}
			
			
			
			//cout << "homvecsize: " << homvec.size() << "    runvec: " << runvec.size() << " runbysnpsize: " << runbysnp.size() << endl;
			
			for(unsigned int j = 0; j<loci.size(); j++){
				if (runvec[j]) {
					runbysnp[j]++;
				}
				
			}
			
		}
		
		//cout << " had " << numrunchk << " SNP in ROH's covering " << pdchk << "  bp" << endl;
		
		thisun++;
		
		
	}
	
	
	
	//cout << "Percs:" << endl;
	for(i = 0; i<loci.size(); i++){//quick and dirty
		if(numindivs < 1) runbysnp[i] = -999;
		else runbysnp[i] = runbysnp[i]/numindivs;
		
		//cout << runbysnp[i] << " ";
	}
	
	
	vector<bool> inaroh;
	
	for(i = 0; i<loci.size(); i++){//quick and dirty
		if(runbysnp[i] >= MIN_HOMOZYGOSITY_SAMPLE_PROP) inaroh.push_back(true);
		else inaroh.push_back(false);
	}
	
	
	return inaroh;
}

//---------------------------------------------
//Function Name:
//Function Definition: runs of homozygosity
void world::rohbydeme(){
	
	
	map<string, vector<int> >::iterator dp = dmap.begin();
	while (dp != dmap.end()) {
		//cout << "*********************************" << endl << endl << "DEME " << dp->first << endl << endl;
		vector<bool> curroh = runsofhomozyg(&dp->second);
		droh[dp->first] = curroh;
		
		
		/*
		 vector<bool> ::iterator bp = curroh.begin();
		 while(bp != curroh.end()){
		 cout << *bp << " ";
		 bp++;
		 }
		 
		 cout << endl << "*********************************" << endl;
		 */
		
		dp++;
	}
	
	//cout << "Size of dblocks: " << dblocks.size() << endl;
	
}

//---------------------------------------------
//Function Name:
//Function Definition: finds runs of homozygosity for the given samples
map <int, vector<bool> > world::rohsliding(vector<int>* y){
	int numhet = 0;
	int nummiss = 0;
	bool homstate = true;
	
	
	map <int, vector<bool> > slideroh;
	
	
	
	vector<slidingwindow> swindows;
	//for each locus
	for(unsigned int i = 0; i<loci.size(); i++){
		//generate a sliding window
		int windowlength = 0;
		int windowsnps = 0;
		for(unsigned int j=i; j<loci.size(); j++){
			if(loci[j].phdist >= MAX_LOCI_JUMP) break; //do not allow a ROH to ho over a centromere or off the end of a chromosome
			
			windowlength += loci[j].phdist;
			windowsnps++;
			if (windowlength >= MIN_SLIDING_WINDOW_LENGTH || windowsnps >= MIN_SNPS_IN_WINDOW)break;
		}
		
		
		if(windowsnps >= MIN_SNPS_IN_WINDOW){
			slidingwindow curwind;
			curwind.start = i;
			curwind.stop = i+windowsnps;
			swindows.push_back(curwind);
		}
		
		//windows.insert(make_pair(i,(i+windowsnps))); // all of this will only be done for rows of snps longer than specified minimu
		
		
		//end for i<size
	}
	
	
	
	
	
	
	
	//go through and grab all double entries. i.e. two chromosomes in segs listing with same name
	vector<int>::const_iterator thisun = y->begin();
	while(thisun != --y->end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator nextun = thisun;
		nextun++;
		
		
		
		if(nmap[*thisun] == nmap[*nextun]){//they have same names thus are two chromosomes from same individual
			
			//Make vector of doubles for homozyg, and for all loci in a widnow
			vector<double> homs (loci.size(), 0.0);
			vector<double> wins (loci.size(), 0.0);
			
			vector<slidingwindow>::iterator wp = swindows.begin();
			while(wp != swindows.end()){
				
				//hom/het state redone for every diploid indiv
				numhet = 0;
				nummiss = 0;
				homstate = true;
				for(int k=wp->start; k<wp->stop; k++){
					
					for(int l=k; l<(k+loci[k].length); l++){
						
						
						if(data[*thisun][l] == 'X' || data[*nextun][l] == 'X') {
							nummiss++;
							break;
						}
						if(data[*thisun][k] !=  data[*nextun][k]) {
							numhet++;
							break;
						}
					}
					
					if(numhet > MAX_HETZ_IN_WINDOW || nummiss > MAX_MISS_IN_WINDOW){
						homstate = false;
						break;
					}
					
					
					//end for k
					
				}
				if(homstate) wp->hom = true;
				else wp->hom = false;
				wp++;
			}
			
			wp = swindows.begin();
			while(wp != swindows.end()){
				//cout << "WiNDOW: ";
				for(int m=wp->start; m<wp->stop; m++){
					//cout << m << " ";
					wins[m]++;
					if (wp->hom)homs[m]++;			
				}	
				//cout << endl;
				wp++;
			}
			
			
			
			
			vector<bool> homvec;
			
			for(unsigned int i = 0; i<loci.size(); i++){
				if(wins[i] < 1.0) homvec.push_back(false); //no ROH posible if no window over it
				else {
					if(homs[i]/wins[i] > MIN_HOM_WIND_THRESH) homvec.push_back(true);
					else homvec.push_back(false);
				}
				
			}
			
			
			int truecount = 0;
			//int chk = 0;
			vector<bool> runvec;
			unsigned int j = 0;
			int physrun = 0;
			vector<bool>::const_iterator hp = homvec.begin();
			while(hp != homvec.end()){
				if(*hp == false){
					
					if(truecount >= MIN_HOMOZYGOSITY_RUN && physrun >= MIN_HOMOZYGOSITY_LENGTH){
						//numrunchk += truecount;
						//pdchk += physrun;
						//cout << "Run of " << truecount << " SNPs and " << physrun << " distance" << endl; 
						for(int i =0; i < truecount; i++){
							runvec.push_back(true);	
							//chk++;
						}
						
					}
					else{
						for(int i =0; i < truecount; i++){
							runvec.push_back(false);
							//chk++;
						}
						
					}
					
					truecount = 0;
					physrun = 0;
					//add on the false at the end
					runvec.push_back(false);
					//chk++;
					
				}
				
				else {
					truecount++;
					physrun += loci[j].phdist;
				}
				j++;
				hp++;	
			}
			
			
			
			
			//If it does not end on a false, read out the final trues
			if(truecount >= MIN_HOMOZYGOSITY_RUN && physrun >= MIN_HOMOZYGOSITY_LENGTH){
				//numrunchk += truecount;
				//pdchk += physrun;
				//cout << "Run of " << truecount << " SNPs and " << physrun << " distance" << endl; 
				for(int i =0; i < truecount; i++){
					runvec.push_back(true);	
					//chk++;
				}
				
				
				
			}
			else{
				for(int i =0; i < truecount; i++){
					runvec.push_back(false);	
					//chk++;
				}
				
				
			}
			
			
			
			//cout << "homvecsize: " << homvec.size() << "    runvec: " << runvec.size() << " runbysnpsize: " << runbysnp.size() << endl;
			
			slideroh.insert(make_pair(*thisun, runvec));
			
			
			//end if nextun == thisun
		}
		
		thisun++;
		
		
	}
	
	
	
	
	return slideroh;
}





//---------------------------------------------
//Function Name:
//Function Definition: runs of homozygosity
void world::rohslidebydeme(){
	
	
	map<string, vector<int> >::iterator dp = dmap.begin();
	while (dp != dmap.end()) {
		//cout << "*********************************" << endl << endl << "DEME " << dp->first << endl << endl;
		map <int, vector<bool> > curroh = rohsliding(&dp->second);
		dslideroh[dp->first] = curroh;
		

		dp++;
	}
	
	//cout << "Size of dblocks: " << dblocks.size() << endl;
	
}




//---------------------------------------------
//Function Name:
//Function Definition: glues together a haplotype from the data
string world::makehaplo(int indiv, int low, int high){
	string haplo;
	for(int lp = low; lp < high; lp++){
		haplo += getdata(indiv, lp);
	}
	return haplo;
}


//---------------------------------------------
//Function Name:
//Function Definition: returns locations segregating sites ina  given group
set<int> world::segsites(vector<int>* y){
	set<int> ss;
	
	//From Efficient Algorithms for Counting and Reporting Segregating Sites in Genomic Sequences, Christodoulakis et al.
	
	vector<int>::const_iterator p = y->begin();
	vector<int>::const_iterator b = y->begin();
	
	p++; //don't need to check the first one
	
	while(p != y->end()){
		
		for(int j=0; j < loci.size(); j++)
		{
			if(ss.find(j)==ss.end() && !datamatch(*p, *b, j, true)) ss.insert(j); //if the locus is not already present in  the set AND does not match the first, it is a segsite
			
			
		}
		
		
		
		p++;
	}
	
	
	return ss;
}


//---------------------------------------------
//Function Name:
//Function Definition: Counts number of segregating sites ina  given group
int world::countsegsites(vector<int>* y){
	
	set<int> segsiteset = segsites(y);
	
	return segsiteset.size();
	
}


//---------------------------------------------
//Function Name:
//Function Definition:
double world::countpairwise(int x, int y){
	
	double pairwisediffs = 0.0;
	
	for(unsigned int i = 0; i<loci.size(); i++){
		if(!datamatch(x, y,i, true))  pairwisediffs++;
	}
	return pairwisediffs;
}


//---------------------------------------------
//Function Name:
//Function Definition: Counts number of airwise differences
double world::meanpairwise(vector<int>* y){
	
	double tot = 0.0;
	double over = 0.0;
	
	//cout << y->size() << endl;
	
	vector<int>::const_iterator p = y->begin();
	while(p != --y->end()){ //end at one before acvtual end as nextun will be the last one
		
		vector<int>::const_iterator q  = p;
		q++;
		while(q != y->end()){ 
			tot += countpairwise(*p, *q);
			
			over++;
			q++;
		}
		
		p++;
	}
	

	
	return tot/over;
	
}




//---------------------------------------------
//Function Name:
//Function Definition: Outputs the line for recombination in msHOt
//just makes hotspots for very long, discontinuous jumps - otherwise averages across the distance
string world::recrateandhotspots(double n0){
	
	long totwithin = 0.0;
	
	map<long, long> hotspots;
	
	
	for(unsigned int i = 0; i<loci.size(); i++){
		if(loci[i].phdist >= MAX_LOCI_JUMP){
			hotspots[totwithin] = loci[i].phdist; //don't add thge actiual length, just make a very big hotspot
			
			
		}
		

		else{
			totwithin += loci[i].phdist;
		}	
		
	}

	
	
	
	double rho = 4.0 * n0 * ADJ_RECOMB_PROB * totwithin; 
	
	if (totwithin <= 0.0) {
		cerr << "WARNING: Rejector2 requires physical distances between sites to calculate recombination (msHOT's -r switch). " << endl;
		cerr << "Please include physical distances in your Rejector2 input file!" << endl;
		abort();
	}
	
	//recombination
	stringstream rline;
	rline << "-r ";
	rline << rho;
	rline << " ";
	rline << totwithin;
	
	
	
	//hotspots calculate crossover rate 
	double backrate = rho/totwithin;
	
	
	
	int lasthotspot = -1;
	// merge
	map<long, long> mergedhotspots;
	map<long, long>::iterator hp = hotspots.begin();
	while(hp != hotspots.end()){
		long mergespot = hp->first;
		long mergedist = hp->second;
		
		
		
		
		if(mergespot <= lasthotspot){
			
			
			long thisdist = hp->second;
			
			map<long, long>::iterator ihp = hp;
			ihp--;
			
			long thispos = ihp->first;
			
			long newdist = thisdist + (ihp->second);
			mergedhotspots[thispos]=newdist;
			
			
			
			
		}
		
		else mergedhotspots[mergespot]=mergedist;
		
		lasthotspot = (mergespot+1); //assumes all hotspots take up one space
		hp++;
	}
	
	if(mergedhotspots.size() > 0){ //only mention hotspots
		rline << " -v ";
		rline << mergedhotspots.size();
		rline << " ";
		
		
		
		map<long, long>::const_iterator sp = mergedhotspots.begin();
		while(sp != mergedhotspots.end()){

			
			
			double jumprate = sp->second * backrate;
			if (jumprate > 0.5) jumprate = 0.5;
			rline << (sp->first) << " " << ((sp->first)+1) << " " << jumprate << " ";

			sp++;
		}
	}

	
	return rline.str();
}


//---------------------------------------------
//Function Name:
//Function Definition: Spits out the line for recombination in msHOt
//just makes hotspots for very long, discontinuous jumps - otherwise averages across the distance
/*SCALED by 1000 version
string world::recrateandhotspots(double n0){
	
	long totwithin = 0.0;
	
	map<long, long> hotspots;
	

	for(unsigned int i = 0; i<loci.size(); i++){
		if(loci[i].phdist >= MAX_LOCI_JUMP){
			hotspots[totwithin] = loci[i].phdist; //don't add thge actiual length, just make a very big hotspot
			
			
		}

	
			
		else{
			totwithin += loci[i].phdist;
			//cout << loci[i].phdist << " " << totwithin << endl;
		}	

	}
	//USING LARGE SEUENCES, SO SCALING BY 1000
	long scaledtotwithin = totwithin/1000;
	
	
	
	double rho = 4.0 * n0 * ADJ_RECOMB_PROB * scaledtotwithin *1000; //scaled up by 1000 to account for scale-down above
	
	//recombination
	stringstream rline;
	rline << "-r ";
	rline << rho;
	rline << " ";
	rline << scaledtotwithin;
	

	
	//hotspots calculate crossover rate 
	double backrate = rho/scaledtotwithin;
	
	
	
	int lasthotspot = -1;
	//scale and merge
	map<long, long> schotspots;
	map<long, long>::iterator hp = hotspots.begin();
	while(hp != hotspots.end()){
		long scaledspot = hp->first/1000;
		long scaleddist = hp->second/1000;
		
		
		cout << scaledspot << " " << scaleddist << endl;
		
		if(scaledspot <= lasthotspot){
			
			
			long thisdist = hp->second/1000;
			
			map<long, long>::iterator ihp = hp;
			ihp--;
			
			long thispos = ihp->first/1000;
			
			long newdist = thisdist + (ihp->second/1000);
			schotspots[thispos]=newdist;
			
			
			
			
		}
		
		else schotspots[scaledspot]=scaleddist;
		
		lasthotspot = (scaledspot+1); //assumes all hotspots take up one space
		hp++;
	}
	
	
	
	
	
	
	
	
	
	
	int wha = 1;
	
	
	
	rline << " -v ";
	rline << schotspots.size();
	rline << " ";
	
	
	
	map<long, long>::const_iterator sp = schotspots.begin();
	while(sp != schotspots.end()){
		
		
		
		double jumprate = sp->second * backrate;
		if (jumprate > 0.5) jumprate = 0.5;
		

		

		
		rline << (sp->first) << " " << ((sp->first)+1) << " " << jumprate << " ";
		
		
		wha++;
		sp++;
	}
	
	
	return rline.str();
}
*/

priordist world::blankpriordist(){
    priordist pt;
    int i;
    
    pt.numpops = dmap.size();
    
    for(i=0; i < pt.numpops; i++){
		pt.popsizes.push_back(distn());

		pt.popsizes[i].name = "uniform";
		pt.popsizes[i].v1 = 0;
		pt.popsizes[i].v2 = 0;
	}
	

	
	for(i=0; i < pt.numpops; i++){
		pt.growthrates.push_back(distn());
		pt.growthrates[i].name = "uniform";
		pt.growthrates[i].v1 = 0;
		pt.growthrates[i].v2 = 0;
	}
    
  
    
    
	pt.msrc.name = "uniform";
	pt.msrc.v1 = 0;
	pt.msrc.v2 = 0;
    
    //just one test
    pt.alphalist["SampleSize"] = 0.1;
    
    return pt;
}

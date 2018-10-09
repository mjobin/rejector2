//------------------------------------------------
//
// REJECTOR2
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropological Sciences
// Stanford University, 2009
//
//------------------------------------------------

#include <fstream>
#include <cctype>
#include <sstream>
#include "rejector.h"
#include "rdeme.h"
#include <cstddef>
#include <utility>




//---------------------------------------------
//Function Name: outfiletop
//Function Definition: prints outfile
void outfiletop(ostream &rout, string inFile, int runlimit){
	
	//Print top
	rout << "Rejector Output" << endl;
	rout << "Infile: " << inFile << endl;
	rout << "Run Limit: " << runlimit << endl;
	time_t fstart;
	time (&fstart);
	struct tm * timeinfo = localtime ( &fstart );
	rout << "Printing started " << asctime (timeinfo) << endl << endl;
	
}



//---------------------------------------------
//Function Name:
//Function Definition: for reading in data NB SINCE SIMCOAL OLY WORKLS WITH ONE SYSTEM AT A TIME SO DOES THIS NOW
priordist world::readinfile(ifstream &in)
{
	string buf, abuf, cbuf, name, pop, systemname, pbuf;
	int i, j, k;
	vector <string> insysTypes;
	vector <double> inmutRates;
	vector <vector <string> > inancStates;
	vector <double> inrecRates;
	vector <int> inphysDist;
	vector <int> inNumLoci, inLength;
	
	
	
	priordist pt;
	
	
	
	
	bool hitdata = false;
	
	in >> buf; // Get the first entry to figure out which kind of infile it is.
	if(buf != "REJECTOR"){
		cerr << "Rejstats: Does not seem to be a Rejector input file.";
		abort();
	}
	
	while(in){//Loop won't end unless it hits Tag
		in >> buf;
		if(buf == "Tag") {
			hitdata = true;
			break;
		}
	}
	
	in >> buf;
	if(buf != "Population") {
		cerr << "Rejstats:Infile Alignment Error. Expecting Population";
		abort();
	}
	if(!hitdata)  {
		cerr << "Rejstats:Infile Error Checking Bounds. Tag never found";
		abort();
	}
	
	getline (in,buf); //rest of Tag pop line
	
	int lc = 0; //line count
	
	
	while(in){
		
		getline (in,buf);
		
		if(!buf.empty()){
			
			//cout << "Line:" << buf << endl;
			lc++;
		}
	}
	
	in.clear();
	
	in.seekg (0, ios::beg); //jump back to the beginning of the file
	
	
	summaryprint = false; //don't print summary unless it's requested
	
	in >> buf; // Get the first entry to figure out which kind of infile it is.
	
	
	if(buf != "REJECTOR"){
		cerr << "Rejector: Error in input file alignment! Expecting: REJECTOR but read " << buf << endl;
		abort();
	}
	
	in >> buf;  //Number
	in >> buf;  //of
	in >> buf;  //Populations
	
	in >> buf; //numpops
	pt.numpops = atoi(buf.c_str());
	
	in >> buf;
	if(buf != "/--Priors") {
		cerr << "Rejector: Error in input file alignment: Expecting: /--Priors but read " << buf << endl;
		abort();
	}
	
	
	//NOW READ IN PRIORS
	
	
	
	int nummatrices, matsize, numhistev;
	vector<double> indVec;
	
	in >> buf; //Population
	in >> buf; //Sizes
	
	
	for(i=0; i < pt.numpops; i++){
		pt.popsizes.push_back(distn());
		in >> buf;
		pt.popsizes[i].name = buf;
		in >> buf;
		pt.popsizes[i].v1 = atoi(buf.c_str());
		in >> buf;
		pt.popsizes[i].v2 = atoi(buf.c_str());
	}
	
	in >> buf;
	if (buf != "Growth"){
		cerr << "Rejector: Error in input file alignment: Expecting: Growth but read " << buf << ". Do you have the same number of populations listed here as stated in Number of popualtions above?" << endl;
		abort();
	}
	
	
	in >> buf; //Rates
	
	for(i=0; i < pt.numpops; i++){
		pt.growthrates.push_back(distn());
		in >> buf;
		pt.growthrates[i].name = buf;
		in >> buf;
		pt.growthrates[i].v1 = atof(buf.c_str());
		in >> buf;
		pt.growthrates[i].v2 = atof(buf.c_str());
	}
	
	in >> buf;
	if (buf != "Number"){
		cerr << "Rejector: Error in input file alignment: Expecting: Number but read " << buf << ". Do you have the same number of populations listed here as stated in Number of popualtions above?" << endl;
		
		abort();
	}
	
	while(buf != "Matrices") in >> buf;
	
	in >> buf;
	nummatrices = atoi(buf.c_str());
	
	for(i=0; i < nummatrices; i++){
		pt.migrmat.push_back(distnmm());
		in >> buf;
		pt.migrmat[i].name = buf;
		
		
		matsize = pt.numpops; 
		
		for(j=0; j < matsize; j++){
			indVec.clear();
			for(k=0; k < matsize; k++){
				in >> buf;
				indVec.push_back(atof(buf.c_str()));
			}
			pt.migrmat[i].v1.push_back(indVec);
		}
		for(j=0; j < matsize; j++){
			indVec.clear();
			for(k=0; k < matsize; k++){
				in >> buf;
				indVec.push_back(atof(buf.c_str()));
			}
			pt.migrmat[i].v2.push_back(indVec);
		}
	}
	
	while(buf != "Events") in >> buf;
	
	in >> buf;
	numhistev = atoi(buf.c_str());
	for(i=0; i < numhistev; i++){
		pt.histev.push_back(distnhe());
		in >> buf;
		pt.histev[i].name = buf;
		for(j=0; j<7; j++){
			in >> buf;
			pt.histev[i].v1.push_back(atof(buf.c_str()));
		}
		for(j=0; j<7; j++){
			in >> buf;
			pt.histev[i].v2.push_back(atof(buf.c_str()));
		}
	}
	
	
	in >> buf;  //Microsat
	in >> buf;  //Range
	in >> buf;
	if (buf != "Constraint"){
		cerr << "Rejector: Error in input file alignment: Expecting: Constraint but read " << buf<< endl; 
		abort();
	}
	
	
	
	
	in >> buf;
	pt.msrc.name = buf;
	in >> buf;
	pt.msrc.v1 = atoi(buf.c_str());
	in >> buf;
	pt.msrc.v2 = atoi(buf.c_str());
	
	
	in >> buf; //testlist
	if (buf != "/--Statlist"){
		cerr << "Rejector: Error in input file alignment: Expecting: /--Statlist but read " << buf<< endl; 
		abort();
	}
	
	
	while(true){
		in >> buf;
		in >> abuf;
		
		if(buf == "/--Data") break;
		if (abuf == "/--Data"){
			cerr << "Rejector: Error in input file alignment: Read end of Statlist where an alpha value should be. All stats have a name and an alpha value on the same line. " << endl; 
			abort();
		}
		
		
		//CREATE ALPHALIST
		if(buf.find("Sub") != string::npos){ //subdivided tests need two entries
			acdc = true; //set global so ancestral and descendant tests done
			// pt.alphalist[buf] = atof(abuf.c_str());
			
			string subbase = buf;
			subbase.replace(subbase.find("Sub"),3,"");
			string subanc = subbase;
			subanc += "-ANC";
			pt.alphalist[subanc] = atof(abuf.c_str());
			string subdec = subbase;
			subdec += "-DEC";
			pt.alphalist[subdec] = atof(abuf.c_str());
			
			
		}
		else pt.alphalist[buf] = atof(abuf.c_str());
		
	}
	
	
	
	
	
	
	//READ IN DATA
	if(buf != "/--Data") {
		cerr << "Rejector: Error in input file alignment: Expecting: /--Data but read " << buf << endl; 
		abort();
	}
	
	if(abuf != "Vnaught") {
		cerr << "Rejstats: Expecting Vnaught";
		abort();
	}
	
	//Vnaught
	in >> buf;
	Vnaught = atof(buf.c_str());
	
	
	in >> buf;
	if (buf != "Loci"){	
		cerr << "Rejstats: Loci expected.";
		abort();
	}
	
	int numloci = 0; //counts number of loci
	
	//LOCI TYPES
	insysTypes.clear();
	
	while(true){
		in >> buf;
		if(buf == "MutRates") break;
		
		numloci++;
		
		if(buf == "STR"){
			insysTypes.push_back("MICROSAT"); //needed to keep consistent with simcoal
		}
		else {
			insysTypes.push_back(buf);
		}
		
	}
	
	
	
	//mut rates now using numloci to make sure all lines up
	inmutRates.clear();
	for(i=0; i < numloci; i++){
		in >> buf;
		inmutRates.push_back(atof(buf.c_str()));
	}
	
	
	in >> buf;
	if(buf != "Ancestral") {
		cerr << "Rejstats:Infile Alignment Error. Expecting Ancestral";
		abort();
	}
	
	//anc state
	inancStates.clear();
	for(i=0; i < numloci; i++){
		in >> buf;
		//remove pesky "'s
		
		while(buf.find_first_of("\"") != string::npos){
			buf.replace(buf.find("\""), 1, "");
		}
		
		vector<string> anclist;
		
		while(!buf.empty()){
			string ancsub = buf.substr(0, buf.find_first_of(","));
			anclist.push_back(ancsub);
			if(buf.find_first_of(",") == string::npos) buf.clear();
			else buf.erase(0, (buf.find_first_of(",")+1)); 
			
		}
		inancStates.push_back(anclist);
		
	}
	
	
	
	
	
	//TEST FOR RECOMB RATE
	
	char pchar;
	// READ TO END OF LINE
	bool rpread = true;
	while(rpread){
		pchar=in.peek();
		
		if((pchar == '\n') || (pchar == '\r')){ //either kind of pesky newline
			
			while((pchar == '\n') || (pchar == '\r')){//read all newlines out
				pchar=in.peek();
				if((pchar == '\n') || (pchar == '\r')) in.get(pchar); //if its not an extra carriage return, put back in buffer to be read on next line
			}
			
			rpread = false;
			
		}
		
		else {
			in.get(pchar);
		}
		
		
		//end while rpread
	}
	
	
	
	
	char rectest = in.peek();
	
	if (rectest == 'R'){ //OPTIONAL recomb rate in put file 10/27/08
		
		inrecRates.clear();
		//Will use system name of above
		
		in >> buf;
		
		if(buf != "RecombRt") {
			cerr << "Recombination Rate Specified but 'RecombRt' entry not found.\n";
			abort();
		}
		
		for(i=0; i < numloci; i++){
			in >> buf;
			inrecRates.push_back(atof(buf.c_str()));
		}
		
		//if recomb rate specified, then put in dummy values for physical distance
		inphysDist.push_back(-999);
		
		
		//END Recomb Rate
	}
	
	
	
	
	if (rectest == 'P'){ //OPTIONAL phys dist put in 4/9/10
		
		inphysDist.clear();
		//Will use system name of above
		
		in >> buf;
		
		if(buf != "PhysDist") {
			cerr << "Physical Disatnce Specified but 'PhysDist' entry not found.\n";
			abort();
		}
		
		for(i=0; i < numloci; i++){
			in >> buf;
			int ipd = atoi(buf.c_str());
			inphysDist.push_back(ipd);
			
			//and convert based on 1cM/Mb
			
			inrecRates.push_back((double)ipd/100000000);
			
		}
		
		
		//END PhysDist
	}
	
	
	
	
	in >> buf;
	if(buf != "NumLoci") {
		cerr << "Rejstats:Infile Alignment Error. Expecting NumLoci. You cannot spoecifyt Rec Rates AND PhysDist!";
		abort();
	}		
	
	for(i=0; i < numloci; i++){
		in >> buf;
		inNumLoci.push_back(atoi(buf.c_str()));
	}
	
	
	in >> buf;
	if(buf != "Length") {
		cerr << "Rejstats:Infile Alignment Error. Expecting Length";
		abort();
	}
	for(i=0; i < numloci; i++){
		in >> buf;
		inLength.push_back(atoi(buf.c_str()));
	}
	
	//Generate Loci vector
	
	int inOffset = 0;
	
	for(unsigned int it = 0; it<insysTypes.size(); it++){
		
		for(unsigned int in = 0; in < inNumLoci[it]; in++){
			//cout << "Offset: " << inOffset << endl;
			
			locinfo newlocinfo;
			newlocinfo.type = insysTypes[it];
			newlocinfo.mutrate = inmutRates[it];
			newlocinfo.ancstate = inancStates[it];
			newlocinfo.recrate = inrecRates[it];
			newlocinfo.phdist = inphysDist[it];

			
			
			//IGNORE Length listing for all but DNA. STR length MUST be 5 for Simcoal
			if(insysTypes[it] == "SNP" || insysTypes[it] == "UEP"){
				newlocinfo.length = 1;
			}
			else if (insysTypes[it] == "STR" || insysTypes[it] == "MICROSAT"){
				newlocinfo.length = 5;
			}
			else if (insysTypes[it] == "DNA"){
				newlocinfo.length = inLength[it];
			}
			newlocinfo.offset = inOffset;
			inOffset += newlocinfo.length;
			loci.push_back(newlocinfo);
			
			
			
			
			
			
		}
	}
	
	//cout << "locilength: " << loci.size() << endl;
	
	
	
	
	in >> buf; //MUST READ Tag!
	
	if(buf != "Tag") {
		cerr << "Rejstats:Infile Alignment Error. Expecting Tag.";
		abort();
	}	
	
	
	in >> buf; //Population
	
	//CREATE THE DATA
	int nl = 0;
	//find total data length
	for(unsigned int lp=0; lp < loci.size(); lp++) {
		//cout << lp << " " << loci[lp].type << " " << loci[lp].length << endl;
		nl += loci[lp].length;
	}
	
	cout << "Total size: " << nl << endl;
	
	datalength = lc;
	datawidth = nl;
	
	
	
	data = new char *[datalength] ;
	for(i = 0 ; i < datalength ; i++ ) data[i] = new char[datawidth];
	
	
	
	int di = 0;
	
	bool missing = false;
	
	while(in){
		in >> name;
		//cout << "Name: " << name << "  ";
		//nmap.insert(pair<string, int>(name, di));
		nmap.push_back(name);
		in >> pop;
		//cout << "Pop: " << pop << ": " << di << endl;
		
		if(dmap.find(pop) == dmap.end()){
			vector<int> newvec;
			dmap[pop] = newvec;
		}
		map<string, vector<int> >::iterator dp = dmap.find(pop);
		dp->second.push_back(di);
		
		
		//dmap.insert(pair<string, int>(pop, di));
		
		
		
		unsigned int dk = 0;
		int lwidth = 0;
		//cout << " Lengths: ";
		for(unsigned int dj=0; dj < loci.size(); dj++){ 
			missing = false;
			
			//cout << loci[dj].length << "  ";
			in >> buf;
			
			//Test lmap
			
			
			
			if(atoi(buf.c_str()) < 0) missing = true;
			if(buf == "N") missing = true;
			if(buf == "-") missing = true;
			if(buf == ".") missing = true;
			if(buf == "NA") missing = true;
			
			/*
			 if(atoi(buf.c_str()) < 0) buf = "X"; //all negative numbers are assumed to mean missing data
			 if(buf == "N") buf = "X"; //using a single missing data code internally
			 if(buf == "-") buf = "X"; //using a single missing data code internally
			 if(buf == ".") buf = "X"; //using a single missing data code internally
			 if(buf == "NA") buf = "X"; //using a single missing data code internally
			 */
			
			
			
			
			
			if(dk > datawidth){
				cerr << "Error reading data at line " << di << ". Data in line longer than allotted space of " << datawidth << "characters." << endl;
				abort();
			}
			
			if(missing){
				for(i=0; i<loci[dj].length; i++){
					//cout << "Entering: " << 'X' << endl;
					data[di][dk] = 'X';
					dk++;
				}
			}
			
			else{
				
				
				
				
				if(loci[dj].type == "MICROSAT"){
					//If the char length isnt 1 there is bad juju so check
					if(buf.length() > loci[dj].length){
						cerr << "Error: For Tag " << name << " length of locus " << (dj+1) << " too large! Read length: " << buf.length() << " Loci length is given in infile as " << loci[dj].length << endl;
						abort();
					}
					int zeroes = loci[dj].length - buf.length();
					
					for(i=0; i<zeroes; i++) {
						//cout << "Entering: " << '0' << endl;
						data[di][dk] = '0';
						lwidth++;
						dk++;
					}
					
					
					lwidth += buf.length();
					if(lwidth > datawidth){
						cerr << "Error: For Tag " << name << " length of locus " << (dj+1) << " too large! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
						abort();
					}
					
					string::iterator sit;
					for ( sit=buf.begin() ; sit < buf.end(); sit++ ){
						//cout << "Entering a: " << *sit << endl;
						data[di][dk] = *sit;
						dk++;
					}
					
					
					
					
				}
				else {
					
					//check if length matches expected length
					if(buf.length() != loci[dj].length){
						cerr << "Error: For Tag " << name << " length of locus! Read: " << buf.length() << " but expected " << loci[dj].length << endl;
						
						abort();
					}
					
					
					lwidth += buf.length();
					if(lwidth > datawidth){
						cerr << "Error: For Tag " << name << " length of locus " << (dj+1) << " too large! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
						abort();
					}
					
					string::iterator sit;
					for ( sit=buf.begin() ; sit < buf.end(); sit++ ){
						//cout << "Entering: " << *sit << endl;
						data[di][dk] = *sit;
						dk++;
					}
					
				}
				
				//end else (not missing)
			}
			
			
			//cout <<loci[dj].type << "\t" << cbuf << endl;
			/*
			 for(unsigned int dk=0; dk < loci[dj].length; dk++){ 
			 
			 in.get(c);
			 
			 cout << c << endl;
			 
			 }
			 */
			
			
			/*
			 
			 */
			
			
			
			
		}
        
        
		
		di++;
		
		
		
		
		//newline and cr eater
		char c;
		while(in.get(c)){
			if(c != '\n' && c != '\r' && c != '\t' && c != ' '){
				in.putback(c);
				break;
			}
		}
	}
	
	
	
	
	
	//check
	/*
	int newcount = 0;
	for(int j=0; j < loci.size(); j++){
		//		cout << loci[j].type << " " << loci[j].length << "| ";
		newcount += loci[j].length;
	}
	//	cout << endl;
	
	cout << "New Count: " << newcount << endl;
	*/
	
	
	for(unsigned int lp=0; lp < loci.size(); lp++) {
		//cout << "Type: " << loci[lp].type << " Offset: " << loci[lp].offset << " Length " << loci[lp].length << endl;
		
		set<string, int> sorted;
		for(i = 0; i < lc; i++){
			
			string entry = getdata(i,lp);
			if(entry[0] != 'X'){
				//loci[lp].lmap.insert(i);
				loci[lp].lmap.push_back(i);
			}
			
		}
		
		
		sort(loci[lp].lmap.begin(), loci[lp].lmap.end());
		
		
		
		/*
		 multimap<string, int> sorted;
		 for(i = 0; i < lc; i++){
		 //cout  << i << ": ";
		 string entry;
		 for(unsigned int j = loci[lp].offset; j < (loci[lp].offset+loci[lp].length); j++){
		 //cout << data[i][j];
		 entry += data[i][j];
		 
		 
		 
		 }
		 
		 //cout << endl;
		 
		 
		 
		 
		 //If first char is "x" don't put it in list... ignore missing
		 string::const_iterator ep = entry.begin();
		 //cout << entry << endl;
		 
		 if(*ep  != 'X'){
		 sorted.insert(make_pair(entry, i));
		 //Add to freq list
		 if(loci[lp].counts.find(entry) == loci[lp].counts.end()){//if it does not have an entry make one
		 loci[lp].counts[entry] = 1;
		 
		 }
		 else{ //increment
		 loci[lp].counts.find(entry)->second++;
		 
		 }
		 
		 
		 }
		 //if(*ep  != 'X') sortvec.push_back(make_pair(entry, i));
		 }
		 
		 
		 
		 
		 
		 //Multimap
		 multimap<string, int>::const_iterator sp = sorted.begin();
		 while(sp != sorted.end()){
		 loci[lp].lmap.push_back(sp->second);
		 //cout << sp->first << ": " << sp->second << endl;
		 
		 
		 
		 sp++;
		 }
		 
		 //Turn Count into freqs
		 map<string, int>::const_iterator fr = loci[lp].counts.begin();
		 
		 while(fr != loci[lp].counts.end()){
		 loci[lp].freqs[fr->first] = (double) fr->second/loci[lp].lmap.size();
		 
		 
		 fr++;
		 }
		 
		 */
		
	}
	
	
	
	
	
	//SORT  dmap and add sample sizes to the priortest
	
	
	map<string, vector<int> >::iterator dp = dmap.begin();
	while(dp != dmap.end()){
		sort(dp->second.begin(), dp->second.end());
		pt.demesampsizes[dp->first] = dp->second.size();
		
		dp++;
	}
	
	
	/*
	 for(unsigned int lp=0; lp < loci.size(); lp++) {
	 cout << "For Locus: " << lp << " Size: " << loci[lp].lmap.size() << endl;
	 
	 map<string, double>::const_iterator fr = loci[lp].freqs.begin();
	 map<string, int>::const_iterator cs = loci[lp].counts.begin();
	 
	 while(fr != loci[lp].freqs.end()){
	 cout << fr->first << " " << cs->first << " " << cs->second << " "<< fr->second << endl;
	 
	 
	 fr++;
	 cs++;
	 }
	 cout << endl;
	 
	 
	 }
	 */
	
	
	/*
	 //VEC MUCH SLOWER than MULYIMAP!!
	 //vector<pair<string, int> >::iterator sv;
	 sort(sortvec.begin(), sortvec.end());
	 vector<pair<string, int> >::iterator sv = sortvec.begin();
	 while(sv != sortvec.end()){
	 loci[lp].lmap.push_back(sv->second);
	 cout << sv->first << ": " << sv->second << endl;
	 
	 sv++;
	 }
	 */
	
	
	/*
	 for(unsigned int lp=0; lp < loci.size(); lp++) {
	 cout << "For Locus: " << lp <<": ";
	 
	 vector<int>::const_iterator lm = loci[lp].lmap.begin();
	 while(lm != loci[lp].lmap.end()){
	 cout << *lm << " ";
	 
	 
	 lm++;
	 }
	 cout << endl;
	 
	 
	 }
	 */
	
	
	
	/*
	 for(i = 0;i < lc ; i++){
	 
	 
	 //need offset here
	 
	 for(int j=0; j < loci.size(); j++){
	 for(unsigned int k = loci[j].offset; j < (loci[j].offset+loci[j].length); j++){
	 
	 cout << data[i][k];
	 
	 }
	 cout << " | ";
	 }
	 cout << endl;
	 }
	 
	 */
	
	
	
	
	
	
	
	
	
	
	
	return pt; //return priors and tests
}


int world::readarpfile(ifstream &in)
{
    string buf;
    int i;
    bool genotypic = false;
    string datatype;
    char achar;
    bool aread = true;
    
    
    //---------------------
	//then read in arp
	in >> buf;
	if(buf != "[Profile]") {
		cerr << "Rejstats:Either the infile is not an .arp file or there is a read error.";
		abort();
	}
    getline (in,buf);
    getline (in,buf);
    getline (in,buf);
    
    
   
    while(true){
        in.get(achar);
        if(achar == '=') break;
    }
    in >> buf;
    if(buf == "1") genotypic = true;
    getline (in,buf);
    
    //If there is genotypic data repeat with the same name 
    int hapgeno = 1; //assume this is haplotypic data unless specified
    if(genotypic) hapgeno = 2;
    
    
    while(true){ //gameticphase
        in.get(achar);
        if(achar == '=') break;
    }
    getline (in,buf);
    
    while(true){//data type
        in.get(achar);
        if(achar == '=') break;
    }
    in >> datatype;    
    getline (in,buf);
    
    while(buf != "[[Samples]]"){
        in >> buf;   
    }
    
   // while(buf != "{") in >> buf;
    
    string demename;
    int demesampsize;
    bool firsttime = true;
    int numsamp = 0;
    int itemsinline = 0;
    unsigned int di = 0; //line counter
    
    Vnaught = 0;
    
    while(buf != "[[Structure]]"){
        
        
       
        while(true){
            in.get(achar);
            if(achar == '=') break;
        }
        in >> buf;
        demename = buf;
        getline (in,buf);

        
        while(true){
            in.get(achar);
            if(achar == '=') break;
        }
        in >> buf;
        demesampsize = atoi(buf.c_str());
        getline (in,buf);
        
        getline (in,buf);
        

       
        string name;
        for(i=0; i<demesampsize; i++){ //by line
            
            if(firsttime){// must calibrate data length on the first pass
                firsttime = false;
                in >> name;
                in >> buf; //number after name
                
                aread = true;
				
				while(aread){
					achar=in.peek();
					if((achar == '\n') || (achar == '\r')){ //either kind of pesky newline
						
						while((achar == '\n') || (achar == '\r')){//read all newlines out
							achar=in.peek();
							if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
						}
						
						aread = false;
						
						//end if
					}
                    
                    else if((achar == '\t') || (achar == ' ')){
                        while(in.get(achar)){
                            if(achar != '\t' && achar != ' '){
                                in.putback(achar);
                                break;
                            }
                        }
                        
                        
                    }
					
					else{
						in >> buf;
                        itemsinline++;
					}
					//end while	
                    
				}
                numsamp++;
                
                if(genotypic) getline (in,buf);
                
                numsamp++;
                
            }
            else{
               
                for(int dh=0; dh<hapgeno; dh++){
                     int itemsinthisline = 0;
                    bool aread = true;
                    if (dh==0){
                        in >> name;
                        in >> buf; //number after name
                    }
                    

                    
                    //cout << numsamp << " " << name << " " << demename << endl;
                    
                    while(aread){
                        achar=in.peek();
                        if((achar == '\n') || (achar == '\r')){ //either kind of pesky newline
                            
                            while((achar == '\n') || (achar == '\r')){//read all newlines out
                                achar=in.peek();
                                if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
                            }
                            
                            aread = false;
                            
                            //end if
                        }
                        
                        else if((achar == '\t') || (achar == ' ')){
                            while(in.get(achar)){
                                if(achar != '\t' && achar != ' '){
                                    in.putback(achar);
                                    break;
                                }
                            }
                            
                            
                        }
                        
                        else{
                           in >> buf;
                        itemsinthisline++;
                        }
                        //end while	
                        
                        
                    }
                    if(itemsinthisline != itemsinline){
                        cerr << "Error! Number of entries in arlequin file line does  not match the first line!" << endl;
                        abort();
                    }
                    //end for hapgeno
                    numsamp++;
                     di++;
                }
                
                
                
            }
            
            
            
        }
        
        in >> buf; // }
        
 
        aread = true;
        while(aread){
            achar=in.peek();
            if((achar == '\n') || (achar == '\r')){ //either kind of pesky newline
                
                while((achar == '\n') || (achar == '\r')){//read all newlines out
                    achar=in.peek();
                    if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
                }
                
                

            }
            else {
                
                aread = false;
                achar=in.peek();
            }
            
            
        }
        if(achar == '[') break;
        
        
       
    
    }
    
   int inOffset = 0;
    for(i=0; i<itemsinline; i++){ //each itme in the line is a locus
        
        locinfo newlocinfo;
        newlocinfo.mutrate = 0;
        string anc = "NA";
        vector<string> vecanc;
        vecanc.push_back(anc);
        newlocinfo.ancstate = vecanc;
        newlocinfo.recrate = 0;
        newlocinfo.phdist = 0;
        
        if(datatype == "STANDARD"){ //Assuming its SNP
            newlocinfo.type = "SNP";
            newlocinfo.length = 1;

        }
        else if (datatype == "MICROSAT"){
            newlocinfo.type = "MICROSAT";
            newlocinfo.length = 5;
            
        }
        else{
            cerr << "ERROR: in arpreadfile: unrecognized data type in .arp file" << endl;
            return 1;
        }
        newlocinfo.offset = inOffset;
        inOffset += newlocinfo.length;
        loci.push_back(newlocinfo);
    
    }
    
    datalength = numsamp;
    int nl = 0;
	//find total data length
	for(unsigned int lp=0; lp < loci.size(); lp++) {
		//cout << lp << " " << loci[lp].type << " " << loci[lp].length << endl;
		nl += loci[lp].length;
	}
	
	//cout << "Total size: " << nl << endl;
	datawidth = nl;
	
	
	
	data = new char *[datalength] ;
	for(i = 0 ; i < datalength ; i++ ) data[i] = new char[datawidth];
    

    
   
    

    // GO TO TOP OF FILE AND READ DATA
    //IF STANDARD DATA IS NOT SNP MAKE IT ABORT
    
    in.seekg(0, ios::beg);
    map<string, vector<int> >::iterator dp = dmap.begin();
    

     dp = dmap.begin();
    di = 0;
    
    //---------------------
	//then read in arp
	in >> buf;
	if(buf != "[Profile]") {
		cerr << "Rejstats:Either the infile is not an .arp file or there is a read error.";
		abort();
	}
    while(buf != "[[Samples]]"){
        in >> buf;   
    }
    
    
    while(buf != "[[Structure]]"){

        while(true){
            in.get(achar);
            if(achar == '=') break;
        }
        in >> buf;
        demename = buf;
        getline (in,buf);
        //cout << dp->first << " " << demename << endl;
        vector<int> newvec;
        dmap[demename] = newvec;
        

        
        while(true){
            in.get(achar);
            if(achar == '=') break;
        }
        in >> buf;
        demesampsize = atoi(buf.c_str());
        getline (in,buf);
        
        getline (in,buf);
        
        string name;
        for(i=0; i<demesampsize; i++){ //by line
            
            
            
            for(int dh=0; dh<hapgeno; dh++){
                
                unsigned int dk = 0; //position in data row
                bool aread = true;
                if (dh==0){
                    in >> name;
                    in >> buf; //number after name
                }
                nmap.push_back(name);
                dmap[demename].push_back(di);

                
                //cout << di << " " << name << " " << demename << " ";
                
                while(aread){
                    achar=in.peek();
                    if((achar == '\n') || (achar == '\r')){ //either kind of pesky newline
                        
                        while((achar == '\n') || (achar == '\r')){//read all newlines out
                            achar=in.peek();
                            if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
                        }
                        
                        aread = false;
                        
                        //end if
                    }
                    
                    else if((achar == '\t') || (achar == ' ')){
                        while(in.get(achar)){
                            if(achar != '\t' && achar != ' '){
                                in.putback(achar);
                                break;
                            }
                        }
                        
                        
                    }
                    
                    else{
                        in >> buf;
                        if(di >= numsamp) {
                            cerr << "ERROR in readarpfile, more lines counted than samples!" << di << " " << numsamp << endl;;
                            abort();
                        }
                        if(dk >= itemsinline){
                            cerr << "ERROR in readarpfile, more entries in line " << di << " than in first line" << dk << " " << itemsinline << endl;
                            abort();
                        }
                        
                        
                        //insert data
                        if(datatype == "STANDARD"){ //Assuming its SNP
                            if(buf.length() > 1) {
                                cerr << "ERROR: in arpreadfile: STANDARD Arlequin data can only be interpeteed as SNPs (length 1 per locus)" << endl;
                                return 1;
                            }
                            
                            if(buf == "0") data[di][dk] = 'N';
                            else if(buf == "1") data[di][dk] = 'A';
                            else if(buf == "2") data[di][dk] = 'C';
                            else if(buf == "3") data[di][dk] = 'G';
                            else if(buf == "4") data[di][dk] = 'T';
                            else {
                                cerr << "ERROR: in arpreadfile: Unreadable SNP data" << endl;
                                return 1;
                            }
                           //cout << data[di][dk] << " ";
                            dk++;
                            


                            
                        }
                        else if (datatype == "MICROSAT"){

                            
                        }
                        else{
                            cerr << "ERROR: in arpreadfile: unrecognized data type in .arp file" << endl;
                            return 1;
                        }
                        
                        
                        
                    }
                    //end while	

                    
                }
 
                
                
                //newline, cr and whitesapce eater for line end
                char c;
                while(in.get(c)){
                    if(c != '\n' && c != '\r' && c != '\t' && c != ' '){
                        in.putback(c);
                        break;
                    }
                }
                
              //  cout << "last" << data[di][itemsinline-1] << endl;
            
                di++; //EACH line of genotypic gets an entry
                //cout << endl;
            
            
            
            
            }
            
            
            
           
            
            
        }

        
        in >> buf; // }
        
        aread = true;
        while(aread){
            achar=in.peek();
            if((achar == '\n') || (achar == '\r')){ //either kind of pesky newline
                
                while((achar == '\n') || (achar == '\r')){//read all newlines out
                    achar=in.peek();
                    if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
                }
                
                
                
            }
            else {
                
                aread = false;
                achar=in.peek();
            }
            
            
        }
        dp++;
        if(achar == '[') break;
        
     //end while buf != structure   
    
    }
    

    //end of data in file
        
        //sort loci
        
        for(unsigned int lp=0; lp < loci.size(); lp++) {
            
            
            set<string, int> sorted;
            for(i = 0; i < datalength; i++){
                
                string entry = getdata(i,lp);
                
                if(entry[0] != 'X'){
                    loci[lp].lmap.push_back(i);
                }
                
            }
            
            
            sort(loci[lp].lmap.begin(), loci[lp].lmap.end());
            
        }
        
        
        
        //Sort the demes
       dp = dmap.begin();
        while(dp != dmap.end()){
            sort(dp->second.begin(), dp->second.end());
            
            dp++;
        }
        
        
        
        
        //leave successfully
        return 0; 
        
    
}
//---------------------------------------------
//Function Name:
//Function Definition: for reading in data with two infiles, an .arp and a summary file
//void world::readinfile(ifstream &in, ifstream &pin)
int world::readarpfile(ifstream &in, world* r)
{
	
	string buf, abuf;
	map<string, int> nameandsamp;
	vector<string> syslist, typelist;
	vector<vector<string> > ancstatelist;
	vector<double> mutratelist;
	vector<double> recratelist;
	int i;
	
	
	//Copy the loci info from r into the current world, as they must be the same
	vector<locinfo>::const_iterator lp = r->loci.begin();
	while(lp != r->loci.end()){
		
		locinfo newlocinfo;
		newlocinfo.type = lp->type;
		newlocinfo.mutrate = lp->mutrate;
		newlocinfo.ancstate = lp->ancstate;
		newlocinfo.recrate = lp->recrate;
		newlocinfo.phdist = lp->phdist;
		newlocinfo.length = lp->length;
		newlocinfo.offset = lp->offset;
		
		//BUT leave lmap until the data's been put in
		
		loci.push_back(newlocinfo);
		
		lp++;
	}
	
	
	
	//Base data size off of r data - Must be the same!
	datalength = r->datalength;
	datawidth = r->datawidth;
	
	
	data = new char *[datalength] ;
	for(int x = 0 ; x < datalength ; x++ ) data[x] = new char[datawidth];
	
	
	
	
	
	//---------------------
	//then read in arp
	in >> buf;
	if(buf != "#Arlequin" )  {
		cerr << "Rejstats:Either the infile is not an .arp file or there is a read error.";
		abort();
	}
	
	
	int numname = 0;
	int di = 0; //line counter
	
	map<string, vector<int> >::const_iterator rdemes = r->dmap.begin();
	
	
	while(in){
		while(buf != "{"){
			in >> buf;
			if(buf == "[[Structure]]"){//end of data in file
				
				//sort loci?
				
				for(unsigned int lp=0; lp < loci.size(); lp++) {
					
					
					set<string, int> sorted;
					for(i = 0; i < datalength; i++){
						
						string entry = getdata(i,lp);
						
						if(entry[0] != 'X'){
							loci[lp].lmap.push_back(i);
						}
						
					}
					
					
					sort(loci[lp].lmap.begin(), loci[lp].lmap.end());
					
				}
				
				
				
				//Sort the demes
				map<string, vector<int> >::iterator dp = dmap.begin();
				while(dp != dmap.end()){
					sort(dp->second.begin(), dp->second.end());
					
					dp++;
				}
				
				
				
				
				//leave successfully
				return 0; 
			}
		}
		
		
		
		//Make new deme entry
		//stringstream dn;
		//dn << numname;
		//string demename = dn.str();
		vector<int> newvec;
		if(rdemes ==r->dmap.end()) {
			cerr<< "Error! For Simcoal deme #: " << numname << " More demes that in real data." << endl;
			abort();
		}
		dmap[rdemes->first] = newvec; //assign deme names back to sim demes
		
		
		
		
		char achar;
		achar=in.peek();
		//read and discard all newlines at front
		while((achar == '\n') || (achar == '\r')){//read all newlines out
			achar=in.peek();
			if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
		}
		
		int entriesindeme = 0;
		
		while(true){//Loop won't end unless it hits }
			achar=in.peek();
			if(achar == '}') break; //will jump out as soon as it hits } since theree is no } until end
			buf.clear();
			
			
			//If there is genotypic data repeat with the same name 
			int hapgeno = 1; //assume this is haplotypic data unless specified
			if(genotypic) hapgeno = 2;
			for(int dh=0; dh<hapgeno; dh++){
				bool aread = true;
				
				while(aread){
					achar=in.peek();
					if((achar == '\n') || (achar == '\r')){ //either kind of pesky newline
						
						while((achar == '\n') || (achar == '\r')){//read all newlines out
							achar=in.peek();
							if((achar == '\n') || (achar == '\r')) in.get(achar); //if its not an extra carriage return, put back in buffer to be read on next line
						}
						
						aread = false;
						
						//end if
					}
					
					else{
						in.get(achar);
						buf += achar;
					}
					//end while	
				}
				//end for hapgeno
			}
			
			
			if(buf.size() > 0){
				// get name
				string::size_type thepos = buf.find_first_of("\t");
				if(thepos == string::npos){
					cerr << "Rejstats:Could not read name of segment from arp file";
					abort();
				}
				string name(buf,0,thepos);
				
				
				
				
				
				
				
				
				buf.erase(0, (buf.find_first_of("\t")+1));
				
				thepos = buf.find_first_of("\t");
				if(thepos == string::npos){
					cerr << "Rejstats:Could not read that second column from arp file";
					abort();
				}
				string secondcol(buf,0,thepos);
				buf.erase(0, (buf.find_first_of("\t")+1));
				
				
				//if genotypic, should have a doublesized in there which gets read into two separate lines of data
				for(int dh=0; dh<hapgeno; dh++){
					
					thepos = buf.find_first_not_of("\t");
					
					buf.erase(0, (thepos));
					
					nmap.push_back(name);
					
					
					//put into dmap
					dmap[rdemes->first].push_back(di);
					
					int lwidth = 0; //total size of a line
					//int dj =0;
					unsigned int dk = 0; //position in data row
					//	while(buf.size()>0){
					//now readin each piece of data
					
					//cout << name << " " << buf << endl;
					//cout << "datawidth: " << datawidth << endl;
					
					
					
					
					
					for(unsigned int dj=0; dj < loci.size(); dj++){
						//cout <<loci[dj].type << " ";
						
						
						
						if(loci[dj].type == "MICROSAT"){
							thepos = buf.find_first_not_of(" ");
							//cout << thepos << endl;
							buf.erase(0, (thepos));
							//cout << "BUF:" << buf << endl;
							
							thepos = buf.find_first_not_of("	");
							buf.erase(0, (thepos));
							
							
							thepos = buf.find_first_of(" ");
							string inseg(buf,0,thepos);
							buf.erase(0, (thepos+1));
							//cout << dj << ":  INSEG: " << inseg << " lwidth: " << lwidth << "  length: " << inseg.length()<< endl;
							
							int zeroes = loci[dj].length - inseg.length();
							
							for(int i=0; i<zeroes; i++) {
								data[di][dk] = '0';
								lwidth++; //TOTAL WIDTH ADDED TO HERE TOO
								dk++;
							}
							
							lwidth += inseg.length();
							if(lwidth > datawidth){
								cerr << "Error: For Tag " << name << " length of locus " << (dj+1) << " too large! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
								abort();
							}
							
							
							string::iterator sit;
							for ( sit=inseg.begin() ; sit < inseg.end(); sit++ ){
								//cout << "Entering a: " << *sit << endl;
								data[di][dk] = *sit;
								dk++;
							}
							
							
							
						}
						
						else if(loci[dj].type == "SNP"){
							thepos = buf.find_first_not_of(" ");
							//cout << thepos << endl;
							buf.erase(0, (thepos));
							//cout << "BUF:" << buf << endl;
							
							
							thepos = buf.find_first_of(" ");
							string inseg(buf,0,thepos);
							buf.erase(0, (thepos+1));
							//cout << "INSEG: " << inseg << endl;
							
							lwidth += inseg.length();
							if(lwidth > datawidth){
								cerr << "Error: For Tag " << name << " length too large at locus " << (dj+1) << "! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
								abort();
							}
							
							
							string::iterator sit;
							for ( sit=inseg.begin() ; sit < inseg.end(); sit++ ){
								//cout << "Entering a: " << *sit << endl;
								data[di][dk] = *sit;
								dk++;
							}
							
							
						}
						
						else if(loci[dj].type == "DNA"){
							
							thepos = buf.find_first_not_of(" ");
							//cout << thepos << endl;
							buf.erase(0, (thepos));
							//cout << "BUF:" << buf << endl;
							
							
							string inseg;
							unsigned int l;
							for(l=0; l<loci[dj].length; l++){
								inseg += buf[l];
								
							}
							buf.erase(0, (l));
							//cout << "INSEG: " << inseg << endl;
							
							lwidth += inseg.length();
							if(lwidth > datawidth){
								cerr << "Error: For Tag " << name << " length too large at locus " << (dj+1) << "! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
								abort();
							}
							
							
							string::iterator sit;
							for ( sit=inseg.begin() ; sit < inseg.end(); sit++ ){
								//cout << "Entering a: " << *sit << endl;
								data[di][dk] = *sit;
								dk++;
							}
							
							
						}
						
						
						
						
						
						
						//end read in each data line loop	
						
						
					}
					//cout << endl;
					di++; //only advance line count if this was a real line
					
					//end for dh hapgeno
				}
				
				
			}
			
			entriesindeme++;
			if(genotypic) entriesindeme++;
			
		}
		if(entriesindeme != rdemes->second.size()){
			cerr << "Error! For deme " << rdemes->first << " size should be " << rdemes->second.size() << " but sim deme size is " << entriesindeme << endl;
			abort();
		}
		
		rdemes++;
		numname++;
		//end whilein
	}
	
	
	
	
	//end fxn
	return 0;
}



//---------------------------------------------
//Function Name:
//Function Definition: Reads rej stats outfile and REPLACES the stats currently in the world
int world::readoutfile(istream &in)
{
	string buf, name, type;
	int i;
	
	vector<string> locallocilist;
	
	//FIRST, make sure the stats in the world have been cleared
	stats.clear();
	
	
	
	//if needed base this on the loci list in the world so we ensure it matches
	
	
	while(buf != "/---") in >> buf;
	
	//test outfile
	//ofstream out("rejectortest.out");
	
	//now take in each test one by one
	while(in){
		while(true){
			
			//New stat to be read, so make one
			statresult x;
			
			int cwidth = 0;
			in >> name; //test name
			
			in >> type; //test type
			
			
			
			
			//cout << name << endl;
			
			//ONEWAY
			
			
			if(type == "Oneway"){
				x.type = "Oneway";
				
				in >> buf; //System
				in >> buf; // ->
				
				sresult o;
				vector<string> onewaysyslist;
				while(true){
					in >> buf;
					if(buf == "Deme#") break;
					onewaysyslist.push_back(buf);
					cwidth++;
				}
				
				//test if thats the same as the locus size
				if(cwidth != loci.size()){
					cerr << "Error in readoutfile. Number of loci does not match for " << name << endl;
					abort();
				}
				
				
				in >> buf; //DemeName
				
				for(i = 0; i<cwidth; i++) in >> buf;
				map<string, map<string, vector<double> > > inDemelist;
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name. 
					
					
					vector<string>::const_iterator osl =  onewaysyslist.begin();
					for(i = 0; i<cwidth; i++){
						if(osl == onewaysyslist.end()){
							cerr << "Error: Oneway result readin. width of data row does not match width of system row!";
							abort();
						}
						
						in >> buf;
						
						
						vector<double> inVec =  inDeme[*osl]; //will create empty if one does not exist accoring to stl tutorial p 177
						//cout << "*osl: " << *osl << "  buf: " << buf << endl;
						
						if(buf == "NA" || buf == "nan") inVec.push_back(-999);
						else inVec.push_back(atof(buf.c_str()));
						
						
						
						
						inDeme[*osl] = inVec;
						osl++;
						
						/*
						 vector<double>::const_iterator ppp = inVec.begin();
						 cout << "invec\n";
						 while(ppp != inVec.end()){
						 cout << *ppp << endl;
						 ppp++;
						 }
						 */
						
						
					}
					inDemelist[demenum] = inDeme;
				}
				x.result["Oneway"] = inDemelist;
				
				
				
				//onewaytestprint(out, o, name);
				
			}
			
			
			else if(type == "Onewaydlist"){
				x.type = "Onewaydlist";
				
				in >> buf; //System
				in >> buf; // ->
				
				sresult o;
				vector<string> onewaysyslist;
				while(true){
					in >> buf;
					if(buf == "Deme#") break;
					onewaysyslist.push_back(buf);
					cwidth++;
				}
				
				in >> buf; //DemeName
				
				for(i = 0; i<cwidth; i++) in >> buf;
				map<string, map<string, vector<double> > > inDemelist;
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name. 
					
					
					vector<string>::const_iterator osl =  onewaysyslist.begin();
					for(i = 0; i<cwidth; i++){
						if(osl == onewaysyslist.end()){
							cerr << "Error: Onewaydlist result readin. width of data row does not match width of system row!";
							abort();
						}
						
						in >> buf;
						
						
						vector<double> inVec =  inDeme[*osl]; //will create empty if one does not exist accoring to stl tutorial p 177
						//cout << "*osl: " << *osl << "  buf: " << buf << endl;
						
						if(buf == "NA" || buf == "nan") inVec.push_back(-999);
						else inVec.push_back(atof(buf.c_str()));
						
						
						
						
						inDeme[*osl] = inVec;
						osl++;
						
						//vector<double>::const_iterator ppp = inVec.begin();
						//cout << "invec\n";
						//while(ppp != inVec.end()){
						//		cout << *ppp << endl;
						//	ppp++;
						//	}
						
						
					}
					inDemelist[demenum] = inDeme;
				}
				x.result["Onewaydlist"] = inDemelist;
				
				//onewaytestprint(cout, o, name);
				
			}
			
			
			//TWOWAY
			
			
			else if (type == "Twoway"){
				x.type = "Twoway";
				sresult o;
				vector<string> versussys;
				vector<string> versustype;
				vector<string> locuspos;
				
				in >> buf; //Versus#
				in >> buf; // ->
				
				while(true){
					in >> buf;
					if(buf == "->") break;
					cwidth++;
				}
				
				
				in >> buf; // DemeName
				
				for(i = 0; i<cwidth; i++){
					in >> buf; //versus names
					versussys.push_back(buf);
				}
				
				map<string, map<string, vector<double> > > inDemelist;
				
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name.  
					
					//cout << "----------" << endl << "For Deme: " << demenum << endl; 
					
					vector<string>::const_iterator vs =  versussys.begin();
					
					
					for(i = 0; i<cwidth; i++){
						if(vs == versussys.end()){
							cerr << "Error: Twoway result readin. width of data row does not match width of versus row!";
							abort();
						}
						
						in >> buf; //data
						
						if(buf != "."){
							vector<double> inVec =  inDeme[*vs]; //will create empty if one does not exist accoring to stl tutorial p 177
							// cout << "*vs: " << *vs << "  buf: " << buf << endl;
							
							if(anydigits(buf.c_str())) inVec.push_back(atof(buf.c_str()));								
							else   inVec.push_back(-999);
							
							/*
							 vector<double>::const_iterator ppp = inVec.begin();
							 cout << "invec\n";
							 while(ppp != inVec.end()){
							 cout << *ppp << endl;
							 ppp++;
							 }
							 */
							
							inDeme["Twoway"] = inVec;
						}
						
						
						
						
						
						inDemelist[*vs] = inDeme;
						
						vs++;
						
						
					}
					
					
					
					
					//end while true
					
					x.result[demenum] = inDemelist;
				}
				
				
				//twowaytestprint(cout, o, name);	
				
			}
			
			
			//ByTwoDeme
			
			
			else if (type == "ByTwoDeme"){
				x.type = "ByTwoDeme";
				sresult o;
				vector<string> versussys;
				vector<string> versustype;
				vector<string> locuspos;
				
				in >> buf; //Versus#
				in >> buf; // ->
				
				while(true){
					in >> buf;
					if(buf == "->") break;
					cwidth++;
				}
				
				
				in >> buf; // DemeName
				
				for(i = 0; i<cwidth; i++){
					in >> buf; //versus names
					versussys.push_back(buf);
				}
				
				map<string, map<string, vector<double> > > inDemelist;
				
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name.  
					
					//cout << "----------" << endl << "For Deme: " << demenum << endl; 
					
					vector<string>::const_iterator vs =  versussys.begin();
					
					
					for(i = 0; i<cwidth; i++){
						if(vs == versussys.end()){
							cerr << "Error: Twoway result readin. width of data row does not match width of versus row!";
							abort();
						}
						
						in >> buf; //data
						
						if(buf != "."){
							vector<double> inVec =  inDeme[*vs]; //will create empty if one does not exist accoring to stl tutorial p 177
							// cout << "*vs: " << *vs << "  buf: " << buf << endl;
							
							if(anydigits(buf.c_str())) inVec.push_back(atof(buf.c_str()));								
							else   inVec.push_back(-999);
							
							/*
							 vector<double>::const_iterator ppp = inVec.begin();
							 cout << "invec\n";
							 while(ppp != inVec.end()){
							 cout << *ppp << endl;
							 ppp++;
							 }
							 */
							
							inDeme["ByTwoDeme"] = inVec;
						}
						
						
						
						
						
						inDemelist[*vs] = inDeme;
						
						vs++;
						
						
					}
					
					
					
					
					//end while true
					
					x.result[demenum] = inDemelist;
				}
				
				
				//twowaytestprint(cout, o, name);	
				
			}
			
			//TWOWAYXLOCI
			
			
			else if (type == "TwowayXloci"){
				x.type = "TwowayXloci";
				sresult o;
				vector<string> versussys;
				vector<string> versustype;
				vector<string> locuspos;
				
				in >> buf; //Versus#
				in >> buf; // ->
				
				while(true){
					in >> buf;
					if(buf == "->") break;
					cwidth++;
				}
				
				
				in >> buf; // DemeName
				
				for(i = 0; i<cwidth; i++){
					in >> buf; //versus names
					versussys.push_back(buf);
				}
				
				map<string, map<string, vector<double> > > inDemelist;
				
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name.  
					
					//cout << "----------" << endl << "For Deme: " << demenum << endl; 
					
					vector<string>::const_iterator vs =  versussys.begin();
					
					
					for(i = 0; i<cwidth; i++){
						if(vs == versussys.end()){
							cerr << "Error: Twoway result readin. width of data row does not match width of versus row!";
							abort();
						}
						
						in >> buf; //data
						
						if(buf != "."){
							vector<double> inVec =  inDeme[*vs]; //will create empty if one does not exist accoring to stl tutorial p 177
							// cout << "*vs: " << *vs << "  buf: " << buf << endl;
							
							if(anydigits(buf.c_str())) inVec.push_back(atof(buf.c_str()));								
							else   inVec.push_back(-999);
							
							/*
							 vector<double>::const_iterator ppp = inVec.begin();
							 cout << "invec\n";
							 while(ppp != inVec.end()){
							 cout << *ppp << endl;
							 ppp++;
							 }
							 */
							
							inDeme["TwowayXloci"] = inVec;
						}
						
						
						
						
						
						inDemelist[*vs] = inDeme;
						
						vs++;
						
						
					}
					
					
					
					
					//end while true
					
					x.result[demenum] = inDemelist;
				}
				
				
				//twowaytestprint(cout, o, name);
				//twowaytestprint(cout, o, name);	
				
			}
			
			
			
			//TWOLOCI
			
			
			else if (type == "TwoLoci"){
				x.type = "TwoLoci";
				sresult o;
				vector<string> twolocdemenum;
				vector<string> twolocsys;
				vector<string> twolocpair;
				
				in >> buf; //Deme#
				
				while(true){
					in >> buf;
					if(buf == "DemeName") break;
					
					cwidth++;
				}
				
				
				int demenamenum = 0;
				while(true){
					in >> buf; //actual deme names. 
					twolocdemenum.push_back(buf);
					if(buf == "System") break;
					demenamenum++;
				}
				
				if(cwidth != demenamenum){
					cerr << "Error: Twoloci test " << name << "number of horizontal deme#'s does not number of deme names";
					abort();
				}
				
				
				int sysnum = 0;
				
				while(true){
					in >> buf;
					if(buf == "Locus") break; 
					twolocsys.push_back(buf);
					sysnum++;
				}
				
				
				if(cwidth != sysnum){
					cerr << "Error: Twoloci test " << name << "number of horizontal deme#'s does not number of systems";
					abort();
				}
				
				in >> buf; //Pair
				
				int pairnum = 0;
				while(true){
					in >> buf;
					if(buf == "Results") break;
					string firstlocus = buf;
					if(firstlocus.find_first_of("<")){
						firstlocus.erase(firstlocus.find_first_of("<"));
					}
					twolocpair.push_back(firstlocus);
					pairnum++;
				}
				
				if(cwidth != pairnum){
					cerr << "Error: Twoloci test " << name << "number of horizontal deme#'s does not number of locus pairs";
					abort();
				}
				
				
				
				vector<string>::const_iterator td = twolocdemenum.begin();
				vector<string>::const_iterator ts = twolocsys.begin();
				vector<string>::const_iterator tp = twolocpair.begin();
				for(i = 0; i<cwidth; i++){
					if(td == twolocdemenum.end()){
						cerr << "Error: Twoloci result readin. width of deme row incorrect!";
						abort();
					}
					if(ts == twolocsys.end()){
						cerr << "Error: Twoloci result readin. width of sys row incorrect!";
						abort();
					}
					if(tp == twolocpair.end()){
						cerr << "Error: Twoloci result readin. width of pair row incorrect!";
						abort();
					}
					
					
					
					
					
					
					in >> buf; //data
					
					
					map<string, map<string, vector<double> > > inDemelist = o[*td];
					map<string, vector<double> > inDeme = inDemelist[*ts];
					
					
					vector<double> inVec =  inDeme[*tp]; //will create empty if one does not exist accoring to stl tutorial p 177
					//cout << "*td: " << *td << "*ts: " << *ts << "*tp: " << *tp << "  buf: " << buf << endl;
					if(buf == "NA" || buf == "nan") inVec.push_back(-999);
					else inVec.push_back(atof(buf.c_str()));
					
					inDeme[*tp] = inVec;
					inDemelist[*ts] = inDeme;
					o[*td] = inDemelist;
					
					td++;
					ts++;
					tp++;
				}
				
				in >> buf; //this is the /--- since there is a single row in loci tests and it isn't read as a delimiter
				
				x.result = o;
				//	twolocitestprint(out, o, name);
				
				
			}
			
			
			//TWOLOCIadj
			
			
			else if (type == "TwoLociAdj"){
				x.type = "TwoLociAdj";
				sresult o;
				vector<string> twolocdemenum;
				vector<string> twolocsys;
				vector<string> twolocpair;
				
				in >> buf; //Deme#
				
				while(true){
					in >> buf;
					if(buf == "DemeName") break;
					
					cwidth++;
				}
				
				
				int demenamenum = 0;
				while(true){
					in >> buf; //actual deme names. 
					twolocdemenum.push_back(buf);
					if(buf == "System") break;
					demenamenum++;
				}
				
				if(cwidth != demenamenum){
					cerr << "Error: Twoloci test " << name << "number of horizontal deme#'s does not number of deme names";
					abort();
				}
				
				
				int sysnum = 0;
				
				while(true){
					in >> buf;
					if(buf == "Locus") break; 
					twolocsys.push_back(buf);
					sysnum++;
				}
				
				
				if(cwidth != sysnum){
					cerr << "Error: Twoloci test " << name << "number of horizontal deme#'s does not number of systems";
					abort();
				}
				
				in >> buf; //Pair
				
				int pairnum = 0;
				while(true){
					in >> buf;
					if(buf == "Results") break;
					string firstlocus = buf;
					if(firstlocus.find_first_of("<")){
						firstlocus.erase(firstlocus.find_first_of("<"));
					}
					twolocpair.push_back(firstlocus);
					pairnum++;
				}
				
				if(cwidth != pairnum){
					cerr << "Error: Twoloci test " << name << "number of horizontal deme#'s does not number of locus pairs";
					abort();
				}
				
				
				
				vector<string>::const_iterator td = twolocdemenum.begin();
				vector<string>::const_iterator ts = twolocsys.begin();
				vector<string>::const_iterator tp = twolocpair.begin();
				for(i = 0; i<cwidth; i++){
					if(td == twolocdemenum.end()){
						cerr << "Error: Twoloci result readin. width of deme row incorrect!";
						abort();
					}
					if(ts == twolocsys.end()){
						cerr << "Error: Twoloci result readin. width of sys row incorrect!";
						abort();
					}
					if(tp == twolocpair.end()){
						cerr << "Error: Twoloci result readin. width of pair row incorrect!";
						abort();
					}
					
					
					
					
					
					
					in >> buf; //data
					
					
					map<string, map<string, vector<double> > > inDemelist = o[*td];
					map<string, vector<double> > inDeme = inDemelist[*ts];
					
					
					vector<double> inVec =  inDeme[*tp]; //will create empty if one does not exist accoring to stl tutorial p 177
					//cout << "*td: " << *td << "*ts: " << *ts << "*tp: " << *tp << "  buf: " << buf << endl;
					if(buf == "NA" || buf == "nan") inVec.push_back(-999);
					else inVec.push_back(atof(buf.c_str()));
					
					inDeme[*tp] = inVec;
					inDemelist[*ts] = inDeme;
					o[*td] = inDemelist;
					
					td++;
					ts++;
					tp++;
				}
				
				in >> buf; //this is the /--- since there is a single row in loci tests and it isn't read as a delimiter
				
				x.result = o;
				//	twolocitestprint(out, o, name);
				
				
			}
			
			
			//PERSYSTEM
			
			
			else if (type == "PerSystem"){
				x.type = "PerSystem";
				
				in >> buf; //Deme#
				in >> buf; //DemeName
				
				sresult o;
				cwidth = 1; // column width equals ONE 7/9/09
				vector<string> persyssyslist;
				for(i = 0; i<cwidth; i++){
					in >> buf;
					persyssyslist.push_back(buf);
				}
				
				
				map<string, map<string, vector<double> > > inDemelist;
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name. 
					
					
					vector<string>::const_iterator psl =  persyssyslist.begin();
					for(i = 0; i<cwidth; i++){
						if(psl == persyssyslist.end()){
							cerr << "Error: Persys result readin. width of data column does not match width of system column!";
							abort();
						}
						
						
						in >> buf; //data
						
						
						vector<double> inVec =  inDeme[*psl]; //will create empty if one does not exist accoring to stl tutorial p 177
						// cout << "*psl: " << *psl << "  buf: " << buf << endl;
						
						if(buf == "NA" || buf == "nan") inVec.push_back(-999);
						else inVec.push_back(atof(buf.c_str()));
						
						inDeme[*psl] = inVec;
						psl++;
						
						//vector<string>::const_iterator ppp = inVec.begin();
						//cout << "invec\n";
						//while(ppp != inVec.end()){
						//	cout << *ppp << endl;
						//ppp++;
						//}
						
					}
					inDemelist[demenum] = inDeme;
				}
				x.result["Persystem"] = inDemelist;
				//persystestprint(out, o, name);
				
			}
			
			
			//BYDEME
			
			
			else if (type == "ByDeme"){
				x.type = "ByDeme";
				
				in >> buf; //Deme#
				in >> buf; //DemeName
				
				sresult o;
				cwidth = 1; // column width equals ONE 7/9/09
				vector<string> persyssyslist;
				for(i = 0; i<cwidth; i++){
					in >> buf;
					persyssyslist.push_back(buf);
				}
				
				
				map<string, map<string, vector<double> > > inDemelist;
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name. 
					
					
					vector<string>::const_iterator psl =  persyssyslist.begin();
					for(i = 0; i<cwidth; i++){
						if(psl == persyssyslist.end()){
							cerr << "Error: Persys result readin. width of data column does not match width of system column!";
							abort();
						}
						
						
						in >> buf; //data
						
						
						vector<double> inVec =  inDeme[*psl]; //will create empty if one does not exist accoring to stl tutorial p 177
						// cout << "*psl: " << *psl << "  buf: " << buf << endl;
						
						if(buf == "NA" || buf == "nan") inVec.push_back(-999);
						else inVec.push_back(atof(buf.c_str()));
						
						inDeme[*psl] = inVec;
						psl++;
						
						//vector<string>::const_iterator ppp = inVec.begin();
						//cout << "invec\n";
						//while(ppp != inVec.end()){
						//	cout << *ppp << endl;
						//ppp++;
						//}
						
					}
					inDemelist[demenum] = inDeme;
				}
				x.result["ByDeme"] = inDemelist;
				//persystestprint(out, o, name);
				
			}
			
			//Perlocusxpop
			
			
			else if(type == "PerLocusXPop"){
				x.type = "PerLocusXPop";
				
				in >> buf; //System
				in >> buf; // ->
				
				sresult o;
				
				
				cwidth = loci.size(); // column width equals number of systems
				vector<string> persyssyslist;
				for(i = 0; i<cwidth; i++){
					in >> buf;
					
				}
				
				
				for(i = 0; i<cwidth; i++){ 
					in >> buf;
					persyssyslist.push_back(buf);
				}
				map<string, map<string, vector<double> > > inDemelist;
				
				map<string, vector<double> > inDeme;
				
				
				vector<string>::const_iterator psl =  persyssyslist.begin();
				for(i = 0; i<cwidth; i++){
					if(psl == persyssyslist.end()){
						cerr << "Error: PerLocusXPop result readin. width of data row does not match width of system row!";
						abort();
					}
					
					in >> buf;
					
					
					vector<double> inVec =  inDeme[*psl]; //will create empty if one does not exist accoring to stl tutorial p 177
					//  cout << "*osl: " << *psl << "  buf: " << buf << endl;
					
					if(buf == "NA" || buf == "nan") inVec.push_back(-999);
					else inVec.push_back(atof(buf.c_str()));
					
					
					
					
					inDeme[*psl] = inVec;
					psl++;
					
					//vector<double>::const_iterator ppp = inVec.begin();
					//cout << "invec\n";
					//while(ppp != inVec.end()){
					//		cout << *ppp << endl;
					//	ppp++;
					//	}
					
					
				}
				inDemelist["PerLocusXPop"] = inDeme;
				
				x.result["PerLocusXPop"] = inDemelist;
				
				//onewaytestprint(out, o, name);
				
				
				in >> buf; // the /---
			}
			
			else if(type == "PerLocusXPopAv"){
				x.type = "PerLocusXPopAv";
				
				in >> buf; //Average
				in >> buf; // ->
				
				in >> buf; //OMG TEH DATA 
				
				sresult o;
				
				map<string, map<string, vector<double> > > inDemelist;
				
				map<string, vector<double> > inDeme;
				
				vector<double> inVec;
				
				if(buf == "NA" || buf == "nan") inVec.push_back(-999);
				else inVec.push_back(atof(buf.c_str()));
				
				
				inDeme["PerLocusXPopAv"] = inVec;
				
				inDemelist["PerLocusXPopAv"] = inDeme;
				
				x.result["PerLocusXPopAv"] = inDemelist;
				
				
				
				in >> buf; // the /---
				
				
			}
			
			
			else if (type == "OnewayAv"){
				x.type = "OnewayAv";
				in >> buf; //Type
				in >> buf; //->
				
				sresult o;
				vector<string> onewaysyslist;
				
				while(true){
					in >> buf;
					if(buf == "Deme#") break;
					onewaysyslist.push_back(buf);
					cwidth++;
				}
				
				in >> buf; //DemeName
				
				map<string, map<string, vector<double> > > inDemelist;
				while(true){
					map<string, vector<double> > inDeme;
					string demenum;
					in >> buf;//12/09/05 from now on, in rejeftor, the "nmae" is the first column, the deme number, so that real and simcoal results match up
					if(buf == "/---") break;
					in >> demenum; //actual deme name. 
					
					vector<string>::const_iterator osl =  onewaysyslist.begin();
					for(i = 0; i<cwidth; i++){
						if(osl ==  onewaysyslist.end()){
							cerr << "Error: Onewayav result readin. width of data column does not match width of system column!";
							abort();
						}
						
						in >> buf; //data
						vector<double> inVec =  inDeme[*osl]; //will create empty if one does not exist accoring to stl tutorial p 177
						// cout << "*osl: " << *osl << "  buf: " << buf << endl;
						
						if(buf == "NA" || buf == "nan") inVec.push_back(-999);
						else inVec.push_back(atof(buf.c_str()));
						
						inDeme[*osl] = inVec;
						osl++;
					}
					
					//cout << "demenum: " << demenum << endl; 
					inDemelist[demenum] = inDeme;
					
				}
				x.result["OnewayAv"] = inDemelist;
				//onewayavtestprint(cout, o, name);
				
				
			}
			
			else if (type == "TESTS"){
				
				
				
				
				
				
				return 0;
			}
			else{
				cerr << "Error,  test " << type << "  has no match!";
				return -1;
			}
			
			//while true
			
			//put the statreult in the map
			stats[name] = x;
			
		}
		
	}
	
	return -1; //shouldnt get to here
}





//---------------------------------------------
//Function Name: 
//Function Definition: constructor creates paramters from prior distribution

int world::makescinfile(world *r, priordist *pt, ofstream &pout, bool poplock) 
{
	string buf;
	
	
	
	double g;
	
	pout << "//Parameters for the coalescence simulation program : simcoal.exe" << endl;
	
	unsigned int i;
	
	//number of samples to simulate = demes = size of popsizes
	pout << pt->popsizes.size() << " samples to simulate" << endl;
	
	map<string, vector<int> >::const_iterator dp = r->dmap.begin();
	
	if(r->dmap.size() != pt->popsizes.size()){
		cerr << "Error in makscinfile! Number of populaitons in input file: " << r->dmap.size() << " but number referred to in priors: " << pt->popsizes.size() << endl;
		abort();
	}
	
	//pop sizes
	pout << "//Population effective sizes  - number of genes 2*diploids" << endl;
	vector<distn>::const_iterator popsizesp = pt->popsizes.begin();
	if(poplock){
		double onepopforall = ((long)((dv8(popsizesp->v1, popsizesp->v2, popsizesp->name, true, false))+0.5));
		while(popsizesp != pt->popsizes.end()){
			tparams.demesize[dp->first] = onepopforall;
			pout << onepopforall << endl;
			popsizesp++;
			dp++;
		}
		
	}
	
	else{
		while(popsizesp != pt->popsizes.end()){
			g = ((long)((dv8(popsizesp->v1, popsizesp->v2, popsizesp->name, true, false))+0.5));
			tparams.demesize[dp->first] =g;
			pout << g << endl;
			popsizesp++;
			dp++;
		}
	}
	
	dp = r->dmap.begin();
	if(r->dmap.size() != pt->demesampsizes.size()){
		cerr << "Error in makscinfile! Number of samples in input file: " << r->dmap.size() << " but number referred to in priors: " << pt->demesampsizes.size() << endl;
		abort();
	}
	
	//Sample Sizes set as equal to sizes of real demes
	pout << "//Samples sizes - number of genes 2*diploids" << endl;
	map<string, int>::const_iterator dss = pt->demesampsizes.begin();
	while(dss != pt->demesampsizes.end()){
		tparams.samplesize[dss->first] = dss->second;
		pout << dss->second << endl;
		dss++;
	}
	
	
	dp = r->dmap.begin();
	if(r->dmap.size() != pt->growthrates.size()){
		cerr << "Error in makscinfile! Number of growth rates in input file: " << r->dmap.size() << " but number referred to in priors: " << pt->growthrates.size() << endl;
		abort();
	}
	
	pout << "//Growth rates	: negative growth implies population expansion" << endl;
	vector<distn>::const_iterator growthratesp = pt->growthrates.begin();
	while(growthratesp != pt->growthrates.end()){
		g = dv8(growthratesp->v1, growthratesp->v2, growthratesp->name, false, false);
		tparams.growthrate[dp->first] =g;
		pout << g << endl;
		growthratesp++;
		dp++;
	}
	
	pout << "//Number of migration matrices : 0 implies no migration between demes" << endl;
	
	if(pt->migrmat.empty()){
		pout << "0" << endl;
	}
	else{
		pout << pt->migrmat.size() << endl;
		
		i = 0;
		
		vector <distnmm>::const_iterator mmptr = pt->migrmat.begin();
		while(mmptr != pt->migrmat.end()){
			pout << "//Migration rates matrix " << i << endl;
			migmatrix inMM;
			vector<vector<double> >::const_iterator v1p = mmptr->v1.begin();
			vector<vector<double> >::const_iterator v2p = mmptr->v2.begin();
			while(v1p != mmptr->v1.end()){
				vector<double> indVec;
				vector<double>::const_iterator v1q = v1p->begin();
				vector<double>::const_iterator v2q = v2p->begin();
				while(v1q !=  v1p->end()){
					double tset = dv8(*v1q, *v2q, mmptr->name, false, true);
					pout << tset << " ";
					indVec.push_back(tset);
					v1q++;
					v2q++;
				}
				inMM.push_back(indVec);
				v1p++;
				v2p++;
				pout << endl;
			}
			mmptr++;
			i++;
			tparams.migmat.push_back(inMM);
		}
		
	}
	
	pout << "//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index" << endl;
	pout << pt->histev.size() << " historical events" << endl;
	//HISTORICAL EVENTS - NOT A LOOP AS SOURCE, SINK, AND MIGRATION MATRIX DO NOT VARY
	vector <distnhe>::const_iterator heptr =pt-> histev.begin();
	
	
	while(heptr != pt->histev.end()){
		vector<double> thishe;
		vector<double>::const_iterator hv1 = heptr->v1.begin();
		vector<double>::const_iterator hv2 = heptr->v2.begin();
		
		g =(long)(dv8(hv1[0], hv2[0], heptr->name, true, false)+0.5);
		thishe.push_back(g);
		pout << g << " "; // timing
		
		g = (int)(heptr->v1[1]+0.5);
		thishe.push_back(g);
		pout << g << " ";
		
		
		if((int) heptr->v1[2] == -1) { //Added to pass -1's for forcing snps 11/20/08
			g =  (int)(heptr->v1[2]);
		}
		else {
			
			g= (int)(heptr->v1[2]+0.5);
		}
		thishe.push_back(g);
		pout << g << " ";
		
		g=  dv8(hv1[3], hv2[3], heptr->name, false, true);
		thishe.push_back(g);
		pout << g << " ";
		
		g = dv8(hv1[4], hv2[4], heptr->name, false, true);
		thishe.push_back(g);
		pout << g << " ";
		
		g =dv8(hv1[5], hv2[5], heptr->name, false, false);
		thishe.push_back(g);
		pout << g << " ";
		
		g = (int)(heptr->v1[6]+0.5);
		thishe.push_back(g);
		pout << g << " ";
		
		pout << endl;
		heptr++;
		
		tparams.histev.push_back(thishe);
	}
	
	
	
	string lasttype = "Nothing";
	double lastmutrate = -1;
	double lastrecrate = -1;
	int lastlength = -1;
	
	
	bool firsttime = true;
	
	
	vector<sc2indchrom> sc2indchroms;
	int numinblock =0;
	vector<sc2block> sc2blocks;
	
	//map<string, vector<block> >::const_iterator sl = blocks.begin();
	

	
	for(unsigned int i=0; i<r->loci.size(); i++){
		
		
		if(r->loci[i].recrate < 0.5){		
			if(r->loci[i].type == lasttype && r->loci[i].mutrate == lastmutrate && r->loci[i].recrate == lastrecrate && r->loci[i].length == lastlength){
				
				
				
				
				numinblock++;
				
				
				//end if same
			}
			
			else {
				if(firsttime){ 
					numinblock++;
					firsttime = false;
				}
				else{
					sc2block newblock;
					newblock.type = lasttype;
					newblock.mutrate = lastmutrate;
					newblock.recrate = lastrecrate;
					newblock.length = lastlength;
					newblock.numinblock = numinblock;
					sc2blocks.push_back(newblock);
					numinblock = 1;
				}
				
				
				lasttype = r->loci[i].type;
				lastmutrate = r->loci[i].mutrate;
				lastrecrate = r->loci[i].recrate;
				lastlength = r->loci[i].length;			
				
			}
		}
		else{ //new ind chrom
			
			
			if(firsttime){
				numinblock++;
				firsttime = false;
				
			}
			/*
			if(numinblock > 1) {//put exisiting block in first, then new one
				sc2block newblock;
				newblock.type = lasttype;
				newblock.mutrate = lastmutrate;
				newblock.recrate = lastrecrate;
				newblock.length = lastlength;
				newblock.numinblock = numinblock;
				sc2blocks.push_back(newblock);
			}
			*/
			
			lasttype = r->loci[i].type;
			lastmutrate = r->loci[i].mutrate;
			lastrecrate = r->loci[i].recrate;
			lastlength = r->loci[i].length;	
			
			sc2block newblock;
			newblock.type = lasttype;
			newblock.mutrate = lastmutrate;
			newblock.recrate = lastrecrate;
			newblock.length = lastlength;
			newblock.numinblock = numinblock;
			sc2blocks.push_back(newblock);
			sc2indchrom newchrom;
			newchrom.blocks = sc2blocks;
			sc2indchroms.push_back(newchrom);
			
			

			
			
			
			sc2blocks.clear();
			numinblock = 1;
			

			
		}
		
		


		

		


	}
	
	if(lastrecrate < 0.5){
		//and the last one
		sc2block newblock;
		newblock.type = lasttype;
		newblock.mutrate = lastmutrate;
		newblock.recrate = lastrecrate;
		newblock.length = lastlength;
		newblock.numinblock = numinblock;
		sc2blocks.push_back(newblock);
		sc2indchrom newchrom;
		newchrom.blocks = sc2blocks;
		sc2indchroms.push_back(newchrom);
		

		
	}
	 
	

	
	//now print blocks
	

	
	
	pout << "//Number of independent loci [chromosome]" << endl;
	pout << sc2indchroms.size() << " 1" << endl;
	
	int chromi = 1;
	
	vector<sc2indchrom>::const_iterator sci = sc2indchroms.begin();
	while(sci != sc2indchroms.end()){
		
		pout << "//Number of contiguous linkage blocks in chromosome " << chromi << " :" << endl;
		pout << sci->blocks.size()<< endl; 
		pout << "//Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters" << endl;
		
		vector<sc2block>::const_iterator scp = sci->blocks.begin();
		while(scp != sci->blocks.end()){
			
			
			if (scp->type == "UEP") pout << "SNP " << scp->numinblock <<  " " << scp->recrate << " "; //lets UEP's mimic SNPs in simcoal 6/1/06 mjj
			else if (scp->type == "DNA") pout << scp->type << " " << scp->length <<  " " << scp->recrate << " ";
			
			else pout << scp->type << " " << scp->numinblock << " " << scp->recrate << " "; //number of sites now put in 10/13/06
			
			if (scp->type == "MICROSAT"){
				pout << scp->mutrate << " 0.00 " << (pt->msrc.v2-pt->msrc.v1);
			}
			else if (scp->type == "SNP" || scp->type == "UEP"){
				
				if(scp->mutrate == -1){ //Added 11/20/08 to pass -1 to shut off mutation through the minor allele freq
					pout << " -1";
				}
				else {							
					pout << " 0.00"; //minor allele can have a 0 frequency
				}
			}
			else if (scp->type == "DNA"){
				pout << scp->mutrate << " 0.33"; // 0.33 means equal chan ce of all three other bases being the one mutated to
			}
			
			
			
			pout << endl;

			scp++;
		
		}
		
		sci++;
		chromi++;
		//end indchroms loop
		
	}
	
	return 0;
}




//---------------------------------------------
//Function Name: 
//Function Definition: constructor creates paramters for ms from prior distribution

string world::makemsline(world *r, priordist *pt, bool poplock, bool macs) 
{
	
	
	
	double g;
	stringstream msl;
	if(macs) msl << "./macs "; 
	else msl << "./msHOT "; //only one "sample" or replication per run
	
	
	
	
	int totnumsamps = 0;
	map<string, int>::const_iterator dss = pt->demesampsizes.begin();
	while(dss != pt->demesampsizes.end()){
		totnumsamps += dss->second;
		dss++;
	}
	msl << totnumsamps << " ";
	
	long totsize = 0;
	if(macs) {
		
		for(unsigned int i = 0; i<r->loci.size(); i++){
				totsize += r->loci[i].phdist;

		}
		msl << totsize << " ";
		
	}
	
	if(!macs) msl << " 1 ";
	
	//map<string, vector<block> >::const_iterator sl = blocks.begin();
	//for(unsigned int i=0; i<r->loci.size(); i++){
	
	double firstmutrate = 0.0;
	
	
	
	//vector<block>::const_iterator ssl = sl->second.begin();
	unsigned int i=0;
	
	//Mutation rate for ms is for all sites
	/*
	 Also suppose that the neutral mutation rate is 10−8 per site per generation 
	 and that we are considering a segment 8,000 base pairs long. In this case, if 
	 we take Nto be 20,000, we have θ = 4  ∗ 20, 000 ∗ 10− 
	 8 
	 ∗ 8, 000 = 6.40. 
	 */
	//So sum all mutation rates
	firstmutrate = r->loci[i].mutrate;
	
	
	
	
	//check if it matches what msHOT needs
	for(i=0; i<r->loci.size(); i++){
		//while(ssl != sl->second.end()){
		if (r->loci[i].type != "SNP" && r->loci[i].type != "UEP"){
			cerr << "Error: SNPs only in msHOT and macs!" << endl; //all loci must be SNPs
			return "error";
		}
		if(r->loci[i].mutrate != firstmutrate){
			cerr << "Error: All mutation rates must match in msHOT and macs!" << endl; //all mutation rates must be the same
			return "error";
			
		}
	}
	
	
	map<string, vector<int> >::const_iterator dp = r->dmap.begin();
	
	//n0 is the DIPLOID size so divide by 2 from rejector
	
	
	
	//Set N0 through mutation parameter and size of first population
	vector<distn>::const_iterator popsizesp = pt->popsizes.begin();
	double n0 =  (long)((((dv8(popsizesp->v1, popsizesp->v2, popsizesp->name, true, false))+0.5))/2.0);
	
	//cout << "n0: " << n0 << endl;


	if(macs){
		msl << "-t " << 4.0*n0*firstmutrate << " ";
		
		
		msl << "-r " << 4.0*n0*ADJ_RECOMB_PROB << " ";
		
		//msl << "-r " << 0 << " ";
	}
	else {
		//limit number of seg sites
		msl << "-s " << r->loci.size() << " ";
		
		
		//Recombination - Calculate rho
		
		string recline = r->recrateandhotspots(n0);
		
		msl << recline << " ";
	}
	
	
	//sample sizes
	dp = r->dmap.begin();
	if(r->dmap.size() != pt->demesampsizes.size()){
		cerr << "Error in makemsinfile! Number of samples in input file: " << r->dmap.size() << " but number referred to in priors: " << pt->demesampsizes.size() << endl;
		abort();
	}
	msl << "-I " << pt->demesampsizes.size() << " ";
	dss =  pt->demesampsizes.begin();
	while(dss !=  pt->demesampsizes.end()){
		
		msl << dss->second << " ";
		tparams.samplesize[dss->first] = dss->second;
		
		dp++;
		dss++;
	}
	
	//The size of the first population will be taken as "N0" and the rest will be modified accordingly	
	
	dp = r->dmap.begin();
	if(r->dmap.size() != pt->popsizes.size()){
		cerr << "Error in makemsinfile! Number of populaitons in input file: " << r->dmap.size() << " but number referred to in priors: " << pt->popsizes.size() << endl;
		abort();
	}
	i=0;
	
	vector<double> cursizes; //retained through historical events so that we can calculate what happens to a population through multiple relative size changes rleative to n0
	
	while(popsizesp !=  pt->popsizes.end()){
		
		if(popsizesp !=  pt->popsizes.begin()){
			g = (long)(((dv8(popsizesp->v1, popsizesp->v2, popsizesp->name, true, false))+0.5)/2.0);
			tparams.demesize[dp->first] =(g*2.0);//restore to chromosomes for consistency, only use diploid when talking to msHOT
			cursizes.push_back(g);
			
			
			msl << "-n " << (i+1) << " " << (g/n0) << " ";
			
			//cout << "pop " << i+1 << ": " << g/n0 << endl;
			
		}
		else{
			tparams.demesize[dp->first] =(n0*2.0); //restore to chromosomes for consistency, only use diploid when talking to msHOT
			cursizes.push_back(n0);
			
		}
		
		dp++;
		popsizesp++;
		i++;
	}
	
	
	//GROWTH RATES
	dp = r->dmap.begin();
	if(r->dmap.size() != pt->growthrates.size()){
		cerr << "Error in makemsline! Number of growth rates in input file: " << r->dmap.size() << " but number referred to in priors: " << pt->growthrates.size() << endl;
		abort();
	}
	i=0;
	vector<distn>::const_iterator gp = pt->growthrates.begin();
	while(gp !=  pt->growthrates.end()){
		g= (dv8(gp->v1, gp->v2, gp->name, false, false));
		
		if(r->dmap.size() == 1){ //Don't know why, but with only one pop msHOT won't accept -g, but can use -G anyway
			msl << "-G " << g << " ";
			
		}
		else{
			
			msl << "-g " << (i+1) << " " << g << " ";
			
		}
		tparams.growthrate[dp->first] =g;
		
		dp++;
		gp++;
		i++;
	}
	
	
	//INITIAL MIGRATION MATRIX
	if(!pt->migrmat.empty()){
		
		msl << "-ma ";
		
		
		
		vector <distnmm>::const_iterator mmptr = pt->migrmat.begin(); //JUST THE FIRST ONE
		migmatrix inMM;
		
		
		vector<vector<double> >::const_iterator v1p = mmptr->v1.begin();
		vector<vector<double> >::const_iterator v2p = mmptr->v2.begin();
		while(v1p != mmptr->v1.end()){
			vector<double> indVec;
			vector<double>::const_iterator v1q = v1p->begin();
			vector<double>::const_iterator v2q = v2p->begin();
			while(v1q !=  v1p->end()){
				double tset = dv8(*v1q, *v2q, mmptr->name, false, true);
				msl << tset << " ";
				indVec.push_back(tset);
				v1q++;
				v2q++;
			}
			inMM.push_back(indVec);
			v1p++;
			v2p++;
			msl << " ";
		}
		tparams.migmat.push_back(inMM);
		
	}
	
	
	
	
	//Hisotical events
	vector <distnhe>::const_iterator heptr =  pt->histev.begin();
	
	int lastmigmat = 0;
	
	double lastrelsize = -1;
	
	while(heptr !=  pt->histev.end()){
		vector<double>::const_iterator hv1 = heptr->v1.begin();
		vector<double>::const_iterator hv2 = heptr->v2.begin();
		
		vector<double> thishe;
		
		
		
		double eventtime = (long)(dv8(hv1[0], hv2[0], heptr->name, true, false)+0.5); // timing
		int source =  (int)(heptr->v1[1]+0.5);
		int sink = (int)(heptr->v1[2]+0.5);
		double prop = dv8(hv1[3], hv2[3], heptr->name, false, true);
		
		//BOTTLENECK  - if it is a -1 it tries to set it to the inverse of the previous hist event
		double relsize;
		if(hv1[4] == -1){ 
			if(heptr == pt->histev.begin()){
				cerr << "Error in makemsline. You cannot set a bottleneck using -1 in relatove size on a historical event on the FIRST historical event." << endl;
				abort();
			}
			relsize = 1/lastrelsize; 
		}
		else {
			
			relsize = dv8(hv1[4], hv2[4], heptr->name, false, true); //Rejector records the size change relative to whatever that population was to that point
			

		}
		
		double sizeratio = (cursizes[sink]/n0);  // For msHOT, first figure out what the ratio of the current size is to n0
		double en = sizeratio * relsize;
		//double actchk = cursizes[sink] * relsize; 
		cursizes[sink] = cursizes[sink] * relsize; //but store relative to its original size
		
		
		
		double newgrowth = dv8(hv1[5], hv2[5], heptr->name, false, false);
		int migid = (int)(heptr->v1[6]+0.5);
		
		thishe.push_back(eventtime);
		thishe.push_back(source);
		thishe.push_back(sink);
		thishe.push_back(prop);
		thishe.push_back(relsize);
		thishe.push_back(newgrowth);
		thishe.push_back(migid);
		
		
		//convert time into units of 4N0 generations
		double fnotime = eventtime / (4*n0);
        double incfnotime = fnotime;
        double fnotimemininc = fnotime / 1000;
        if(fnotimemininc < 0.00001){
            fnotimemininc = 0.00001;
        }
		
		//check for movement
		if(source != sink){ //some movement
			if(prop == 1) { //-ej
				msl << "-ej " << incfnotime << " " << (source+1) <<  " " << (sink+1) << " ";
                incfnotime += fnotimemininc;
			}
			else if (prop == 0){
			}
			else {
				cout << "Error: msHOT does not recognize historical events with partial migration. Reconfigure as a migraiton matrix." << endl; // cannot have halfway mig events FIX INTO MIG RATES LATER
				return "error";
			}
			//end if source != sink
		}
		
		
		
		//NEW SIZE
		if(relsize != 1) { // if there is some change in SIZE // SETS GRWOTH TO ZERO
			msl << "-en " << incfnotime <<  " " <<  (sink+1) << " " << en << " ";
            incfnotime += fnotimemininc;
		}
		
		//NEW GROWTH RATE
		if(true) { // if there is some change in growth // JUST PUT IN ALL??
			msl << "-eg " << incfnotime <<  " " <<  (sink+1) << " " << newgrowth << " ";
            incfnotime += fnotimemininc;
		}
		
		//NEW MIGMAT
		if(migid != lastmigmat) { // if there is some change in growth // JUST PUT IN ALL??
			lastmigmat = migid;
			int numpops = r->dmap.size();
			msl << "-ema " << " " << incfnotime << " " << numpops << " ";
            incfnotime += fnotimemininc;
			
			if(migid > pt->migrmat.size()){
				cout << "Error: You appear to have specified a change to a migration matrix " << migid << " that does not exist." << endl; // cannot have halfway mig events FIX INTO MIG RATES LATER
				return "error";
				
			}
			vector <distnmm>::const_iterator mmptr = pt->migrmat.begin(); //JUST THE FIRST ONE
			for(i=0; i< migid; i++) mmptr++;
			
			vector<vector<double> >::const_iterator v1p = mmptr->v1.begin();
			vector<vector<double> >::const_iterator v2p = mmptr->v2.begin();
			while(v1p != mmptr->v1.end()){
				vector<double>::const_iterator v1q = v1p->begin();
				vector<double>::const_iterator v2q = v2p->begin();
				while(v1q !=  v1p->end()){
					double tset = dv8(*v1q, *v2q, mmptr->name, false, true);
					msl << tset << " ";
					
					v1q++;
					v2q++;
				}
				
				v1p++;
				v2p++;
				msl << " ";
			}
		}
		
		
		
		lastrelsize = relsize;
		tparams.histev.push_back(thishe);	
		heptr++;
		
		
	}
	
	
	
	//if(macs) msl << " 2>/dev/null | ./msformatter"; 
	if(macs) msl << " 2>/dev/null"; 
	
	
	//end outputms
	return msl.str();
}


//---------------------------------------------
//Function Name:
//Function Definition: Read ms output from stdin (piped)
int world::readmacsunformatted(world* r, string msline){
	double mutrate = 0.0;
	double recrate = 0.0;
	int segsites = 0;
	double initmig =0.0;
	vector<int> nbypop;
	map<int, double> popsizes;
	map<int, double> grrates;
	vector <double> inmutRates;
	vector <vector <string> > inancStates;
	string buf;
	
//    cout << "Using macs command: " << endl;
//    cout << msline << endl;
	//RUN MACS
	char c;
	string msout;
	FILE *msp;
	msp = popen(msline.c_str(), "r");
	do {
		c = fgetc (msp);
		if (c != EOF) msout += c;
	} while (c != EOF);
	pclose(msp);
	
	//cout << "Macs complete." << endl;
	
	istringstream ciss(msout, istringstream::in);
	string mscmdline;
	getline(ciss, mscmdline);
	//cout << "CMDLINE:" << mscmdline << endl;
	
	
	istringstream mcmd(mscmdline, istringstream::in);
	deque<string> mstokens;
	string msclpart;
	while( mcmd >> msclpart )     
	{
		
		//cout << msclpart << endl;
		mstokens.push_back(msclpart);
		
	}
	
	
	if(mstokens.empty()){ //some mistake made in ms, no output
		cerr << "ERROR in runsim: No output from MaCS" << endl;
		return -1;
		
	}
	
	string mspart;
	
	
	//the command
	mspart = mstokens.front();
	//cout << mspart << endl;
	if(mspart != "COMMAND:"){
		cerr << "ERROR in runsim: First word of read output is not COMMAND:" << endl;
		return -1;
	}
	mstokens.pop_front();
	
	//the command
	mspart = mstokens.front();
	//cout << mspart << endl;
	if(mspart != "./macs"){
		cerr << "ERROR in runsim: Second word of read output is not ./macs" << endl;
		return -1;
	}
	mstokens.pop_front();
	
	//Next is total number of samples
	//i++;
	int tot = atoi(mstokens.front().c_str());
	//cout << "Tot samples from macs: " << tot << endl;
	mstokens.pop_front();
	
	//Next is total number of bp
	int basepairs = atoi(mstokens.front().c_str());
	//cout << "Tot basepairs from macs: " << basepairs << endl;
	
	mstokens.pop_front();
	
	
	
	string msw;
	
	//Teh branchz0rs ms line reader starts here
	while (!mstokens.empty()){
		//Grab the first part. That will tell you what to do with the rest
		msw = mstokens.front();
		mstokens.pop_front();
		
		
		if(msw == "-t"){
			mutrate  = atof(mstokens.front().c_str());
			//cout << "Theta now: " << theta << endl;
			mstokens.pop_front();
		}
		else if (msw == "-s"){
			segsites = atoi(mstokens.front().c_str());
			mstokens.pop_front();
		}
		else if (msw == "-I"){
			int npop = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			for(int i =0; i<npop; i++){
				nbypop.push_back(atoi(mstokens.front().c_str()));
				
				mstokens.pop_front();
			}
			//check if migration rate added
			string migchk = mstokens.front();
			if(migchk.find_first_of("-") != string::npos){
				
			}
			else{
				initmig = atof(mstokens.front().c_str());
				mstokens.pop_front();
			}
			
		}
		else if (msw == "-n"){
			int npop = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			double popsz = atof(mstokens.front().c_str());
			mstokens.pop_front();
			popsizes[npop] = popsz;
		}
		
		else if (msw == "-g"){
			int npop = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			double grrate = atof(mstokens.front().c_str());
			mstokens.pop_front();
			grrates[npop] = grrate;
		}
		
		else if (msw == "-G"){
			double grrate = atof(mstokens.front().c_str());
			mstokens.pop_front();
			int gpr = 0;
			map<int, double>::iterator gp = grrates.begin();
			while(gp != grrates.end()){
				grrates[gpr] = grrate;
				gp++;
				gpr++;
			}
		}
		
		else if (msw == "-ej"){
			
			mstokens.pop_front();
			
			mstokens.pop_front();
			mstokens.pop_front();
		}
		
		else if (msw == "-eg"){
			
			mstokens.pop_front();
			
			mstokens.pop_front();
			mstokens.pop_front();
		}
		
		
		else if (msw == "-en"){ //pop size change
			mstokens.pop_front();
			
			mstokens.pop_front();
			mstokens.pop_front();
			
		}
		
		
		else if (msw == "-r"){ //recrate
			recrate = atof(mstokens.front().c_str());
			mstokens.pop_front();
			
			
		}
		
		else if (msw == "-v"){ //hotspots
			int nhot = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			for(int i =0; i<nhot; i++){
				mstokens.pop_front();
				mstokens.pop_front();
				mstokens.pop_front();
			}
			
			
		}
		
		else if (msw == "-ma"){ //migration
			
			//double newmig = 0.0;
			string newmigstr = "nothing";
			int migwidth = 0;
			
			mstokens.pop_front(); //x
			mstokens.pop_front(); //x
			
			
			while(newmigstr != "x"){
				newmigstr = mstokens.front();
				mstokens.pop_front();
				if(newmigstr == "x") break;
				//double newmig = atof(newmigstr.c_str());
				migwidth++;
				
			}
			

			
			for(int nlen=1; nlen < migwidth; nlen++){
				mstokens.pop_front(); //x	
				
				cout << endl << endl << mstokens.front() << endl;
				
				for(int nwid=0; nwid < migwidth; nwid++) mstokens.pop_front(); //not actually reading them in yet
				
				mstokens.pop_front(); //x	
				
			}
			
			cout << endl << endl << mstokens.front() << endl;
		}
			
			
		
		else if (msw == "0"){ //seems to attach itself when redirect is added
			
		}
		
		else {
			
			cerr << "ERROR in readmsline: Do not recognize command-line switrch " << msw << endl;
			return -1;
			break;
		}
		
		//end while mstoken not empty
	}
	
	
	//cout << "Command line tokenized and read" << endl;
	
	
	string macsline;
	
	bool firstindiv = true;
	datawidth = 0;
	
	vector<string> thedata;
	
	int totalsamples = -1;
	int totalsites = -1;
	
	
	while (ciss) { //run to end
		
		getline(ciss, macsline);
//        cout << "mline: " << macsline << endl;
		istringstream mline(macsline, istringstream::in);
		deque<string> macstokens;
		string macspart;
		while( mline >> macspart )     
		{
			macstokens.push_back(macspart);
		}
		
		string mf;
		while (!macstokens.empty()){
			
			//Get first token, to tell you what sort of line it is
			
			mf = macstokens.front();
			macstokens.pop_front();
			
			
		if(mf == "SEED:"){
			macstokens.pop_front();
		}
		else if(mf == "SITE:"){
			macstokens.pop_front();
			macstokens.pop_front();
            macstokens.pop_front();
			

			
			if(firstindiv){//First person through check 
				
				mf = macstokens.front();
				
				datalength = mf.size();

				if(datalength != tot){
					cerr << "ERROR in readmacsunformatted. Number of read samples " << datalength << " does not match the number of samples " << tot << " specified on the macs command line." << endl;
					abort();
				}
				
				
				
			}
			

			thedata.push_back(mf);
			macstokens.pop_front();

			datawidth++;
			
			
			
		}
		else if(mf == "TOTAL_SAMPLES:"){
			totalsamples = atoi(macstokens.front().c_str());
			macstokens.pop_front();
				
			}
			
		else if(mf == "TOTAL_SITES:"){
			totalsites = atoi(macstokens.front().c_str());
			macstokens.pop_front();
			
		}
			
			

			
			
		}
		
		
	}
	
	
	//Assign all to one population if -I not used
	if(nbypop.empty()) nbypop.push_back(tot);
	
	//Check that the sum of the pops adds to total
	int totchk = 0;
	vector<int>::const_iterator nbp = nbypop.begin();
	while(nbp != nbypop.end()){
		totchk += *nbp;
		
		nbp++;
	}
	if(totchk != tot){
		cerr << "ERROR in rejstats :Total calculated from populations in -I switch not equal to nsam "  << endl;
		return -1;
	}
	
	
	map<string, int>::const_iterator ssp = tparams.samplesize.begin();
	
	if(tparams.samplesize.size() != nbypop.size()){
		cerr << "Error in readmslime! Number of samples in input file: " << tparams.samplesize.size() << " but number referred to in priors: " << nbypop.size() << endl;
		abort();
	}
	
	
	if (datalength != totalsamples) {
		cerr << "ERROR in readmacsunformatted. Number of read samples " << datalength << " does not match the number of samples " << totalsamples << " specified in macs output TOTAL_SAMPLES:." << endl;
		abort();
	}
	
	if (datawidth != totalsites) {
		cerr << "ERROR in readmacsunformatted. Number of read sites " << datawidth << " does not match the number of samples " << totalsites << " specified in macs output TOTAL_SITES:." << endl;
		abort();
	}
	
	
	
	//Create data
	data = new char *[datalength] ;
	for(int x = 0 ; x < datalength ; x++ ) data[x] = new char[datawidth];
	
	
	//cout << "Sim data created with datalength " << datalength << " and datawidth " << datawidth << endl;
	long arraysize = 0;
	for(int d=0; d<datalength; d++){
		for(int e=0; e<datawidth; e++){
			arraysize+=sizeof(data[d][e]);
		}
	}
	//cout << "Calculated zie: " << arraysize << endl;
		
	
	segsites = datawidth;
	
	//add loci
	//************************
	vector<string> newas;
	newas.push_back("A");
	
	int curlength = 0;
	
	for(int i=0;i<segsites;i++){
		

		
		locinfo newlocinfo;
		newlocinfo.type = "SNP";
		newlocinfo.mutrate = mutrate; //FIXME div by 4N?
		newlocinfo.ancstate = newas;
		newlocinfo.recrate = recrate; //FIXME div by 4N?
		newlocinfo.phdist = basepairs/segsites; //FIXME - be more exact site to site?
		newlocinfo.length = 1;
		newlocinfo.offset = curlength;
		curlength += newlocinfo.length;
		
		
		loci.push_back(newlocinfo);
		
		
	}
	
	
	
	//names
	
	//************************
	
	//cout << "Begin reading data" << endl;
	
	map<string, vector<int> >::const_iterator rdemes = r->dmap.begin();
	int di = 0;
	int i = 0;
	nbp = nbypop.begin();
	ssp = tparams.samplesize.begin();
	while(nbp != nbypop.end()){ //iterator through the populations
		
		
		
		int dpop = *nbp; //THIS IS THE SAMPLE SIZE
		//cout << "Reading Pop: " << i << " which has " << dpop << " individuals." << endl;
		
		
		
		vector<int> newvec;
		if(rdemes ==r->dmap.end()) {
			cerr<< "Error! For macs deme #: " << i << " More demes that in real data." << endl;
			abort();
		}
		dmap[ssp->first] = newvec; //assign deme names back to sim demes
		
		
		
		
		
		
		for(int j = 0; j < dpop ; j++){//for each individual in the population
			
			
			//construct a name
			stringstream inn;
			if(genotypic){
				if(j%2==0) inn << i << "_" << j;
				else inn << i << "_" << (j-1);	
			}
			else inn << i << "_" << j;
			string inname = inn.str();
			//cout << "d00d#: " << inname << endl;
			
			
			

			//cout << endl;
			
			nmap.push_back(inname);
			
			//put into dmap
			dmap[ssp->first].push_back(di);
			
			
			di++;
			//end j< dpop	
		}
		
		
		
		//end while nbp !=nbypop
		rdemes++;
		nbp++;
		ssp++;
		i++;
	}
	
	
	//Add in the data
	
	di =0;
	
//    cout << "thedata size: " << thedata.size() << endl;
	
	vector<string>::iterator td = thedata.begin();
	while(td != thedata.end()){
//        cout << "Locus: " << di << " " << td->size() << endl;


        for(int k=0; k<td->size(); k++){
            string inlocus = (td->substr(k, 1));
//            cout << "inlocus: " << inlocus << " ";
            
            char thesnp;
            
//            cout << k << "| ";
            
            if(inlocus == "1")thesnp = 'T';
            else if(inlocus == "0")thesnp = 'A';
            else {
                cerr << "ERROR in readmacsunfornmatted: Do not recognize data: " << inlocus  << endl;
                return -1;
            }
            
//            //check if length matches expected length
//            if(loci[k].length != 1){
//                cerr << "Error: length of locus! Expected " << loci[k].length << endl;
//                abort();
//            }
            
//            cout << di << " ";
            
//            cout << "Entering: " << thesnp << " at " << k << "-" << di << " \r" << std::flush;
            data[k][di] = thesnp;
		
		
	}
		
		di++;
		td++;
	}
	
	
	
	//sort loci?
	
	//cout << "Done data input" << endl;
	
	
	
	for(unsigned int lp=0; lp < loci.size(); lp++) {
		
		//cout << "Locus: " << lp << " | ";
		
		for(i = 0; i < datalength; i++){
			
			string entry = getdata(i,lp);
			//cout << entry << " ";
			
			if(entry[0] != 'X'){
				loci[lp].lmap.push_back(i);
			}
			
		}
		
		
		sort(loci[lp].lmap.begin(), loci[lp].lmap.end());
		//cout << endl;
		
		
	}
	
	
	//cout << "Done making lmaps" << endl;
	
	//Sort the demes
	map<string, vector<int> >::iterator dp = dmap.begin();
	while(dp != dmap.end()){
		sort(dp->second.begin(), dp->second.end());
		
		dp++;
	}
	
	//cout << "Done done sorting demes" << endl;
	
	
	dp = dmap.begin();
	while(dp != dmap.end()){
		//cout << endl <<"______________" << endl << "DEME: " << dp->first <<endl;
		
		vector<int>::const_iterator ddp = dp->second.begin();
		while(ddp != dp->second.end()){
			//cout << *ddp << "| ";
			for(i=0; i<loci.size(); i++){
				string entry = getdata(*ddp,i);
				//cout << entry << " ";
			}
			//cout << endl;
			ddp++;
		}
		
		dp++;
	}
	
	
	return 0;
}

//---------------------------------------------
//Function Name:
//Function Definition: Read ms output from stdin (piped)
int world::readmacsinfile(world* r, string msline){
	double mutrate = 0.0;
	double recrate = 0.0;
	int segsites = 0;
	double initmig =0.0;
	vector<int> nbypop;
	map<int, double> popsizes;
	map<int, double> grrates;
	vector <double> inmutRates;
	vector <vector <string> > inancStates;
	string buf;
	
	
	
	//cout << "Begin readmacsinfile" << endl;
	
	

	
	

	
	
	
	
	
	//RUN MACS
	char c;
	string msout;
	FILE *msp;
	msp = popen(msline.c_str(), "r");
	do {
		c = fgetc (msp);
		if (c != EOF) msout += c;
	} while (c != EOF);
	pclose(msp);
	
	//cout << "Macs complete." << endl;
	
	//TOKENIZE THE STRINGSTREAM

	
	istringstream ciss(msout, istringstream::in);
	string mscmdline;
	getline(ciss, mscmdline);
	//cout << "CMDLINE:" << mscmdline << endl;
	
	istringstream mcmd(mscmdline, istringstream::in);
	deque<string> mstokens;
	string msclpart;
	while( mcmd >> msclpart )     
	{
		
		//cout << msclpart << endl;
		mstokens.push_back(msclpart);
		
	}
	
	
	if(mstokens.empty()){ //some mistake made in ms, no output
		cerr << "ERROR in runsim: No output from MaCS" << endl;
		return -1;
		
	}
	
	
	
	

	
	
	
	
	

	string mspart;
	

		
		
		//the command
		mspart = mstokens.front();
		//cout << mspart << endl;
		if(mspart != "./macs"){
			cerr << "ERROR in runsim: First word of read output is not ./macs" << endl;
			return -1;
		}
		mstokens.pop_front();
		
		//Next is total number of samples
		//i++;
		int tot = atoi(mstokens.front().c_str());
		//cout << "Tot samples from macs: " << tot << endl;
		mstokens.pop_front();
		
		//Next is total number of bp
		int basepairs = atoi(mstokens.front().c_str());
		//cout << "Tot basepairs from macs: " << basepairs << endl;

		mstokens.pop_front();
		
		
		
		string msw;

		//Teh branchz0rs ms line reader starts here
		while (!mstokens.empty()){
			//Grab the first part. That will tell you what to do with the rest
			msw = mstokens.front();
			mstokens.pop_front();
			
			
			if(msw == "-t"){
				mutrate  = atof(mstokens.front().c_str());
				//cout << "Theta now: " << theta << endl;
				mstokens.pop_front();
			}
			else if (msw == "-s"){
				segsites = atoi(mstokens.front().c_str());
				mstokens.pop_front();
			}
			else if (msw == "-I"){
				int npop = atoi(mstokens.front().c_str());
				mstokens.pop_front();
				for(int i =0; i<npop; i++){
					nbypop.push_back(atoi(mstokens.front().c_str()));
					
					mstokens.pop_front();
				}
				//check if migration rate added
				string migchk = mstokens.front();
				if(migchk.find_first_of("-") != string::npos){
					
				}
				else{
					initmig = atof(mstokens.front().c_str());
					mstokens.pop_front();
				}
				
			}
			else if (msw == "-n"){
				int npop = atoi(mstokens.front().c_str());
				mstokens.pop_front();
				double popsz = atof(mstokens.front().c_str());
				mstokens.pop_front();
				popsizes[npop] = popsz;
			}
			
			else if (msw == "-g"){
				int npop = atoi(mstokens.front().c_str());
				mstokens.pop_front();
				double grrate = atof(mstokens.front().c_str());
				mstokens.pop_front();
				grrates[npop] = grrate;
			}
			
			else if (msw == "-G"){
				double grrate = atof(mstokens.front().c_str());
				mstokens.pop_front();
				int gpr = 0;
				map<int, double>::iterator gp = grrates.begin();
				while(gp != grrates.end()){
					grrates[gpr] = grrate;
					gp++;
					gpr++;
				}
			}
			
			else if (msw == "-ej"){
				
				mstokens.pop_front();
				
				mstokens.pop_front();
				mstokens.pop_front();
			}
			
			else if (msw == "-eg"){
				
				mstokens.pop_front();
				
				mstokens.pop_front();
				mstokens.pop_front();
			}
			
			
			else if (msw == "-en"){ //pop size change
				mstokens.pop_front();
				
				mstokens.pop_front();
				mstokens.pop_front();
				
			}
			
			
			else if (msw == "-r"){ //recrate
				recrate = atof(mstokens.front().c_str());
				mstokens.pop_front();
							  
				
			}
			
			else if (msw == "-v"){ //hotspots
				int nhot = atoi(mstokens.front().c_str());
				mstokens.pop_front();
				for(int i =0; i<nhot; i++){
					mstokens.pop_front();
					mstokens.pop_front();
					mstokens.pop_front();
				}
				
				
			}
			else if (msw == "0"){ //seems to attach itself when redirect is added
				
			}
				
			else {
				
				cerr << "ERROR in readmsline: Do not recognize command-line switrch " << msw << endl;
				return -1;
				break;
			}
			
			//end while mstoken not empty
		}
		
		
		//cout << "Command line tokenized and read" << endl;

	//Assign all to one population if -I not used
	if(nbypop.empty()) nbypop.push_back(tot);
	
	//Check that the sum of the pops adds to total
	int totchk = 0;
	vector<int>::const_iterator nbp = nbypop.begin();
	while(nbp != nbypop.end()){
		totchk += *nbp;
		
		nbp++;
	}
	if(totchk != tot){
		cerr << "ERROR in rejstats :Total calculated from populations in -I switch not equal to nsam "  << endl;
		return -1;
	}
	
	
	map<string, int>::const_iterator ssp = tparams.samplesize.begin();
	
	if(tparams.samplesize.size() != nbypop.size()){
		cerr << "Error in readmslime! Number of samples in input file: " << tparams.samplesize.size() << " but number referred to in priors: " << nbypop.size() << endl;
		abort();
	}
	

	

	
buf = "notsegsites";

while (buf != "segsites:") {
	ciss >> buf;
	//cout << buf << endl;
}

	//segsites
	ciss >> buf;
//Number of segregating sites - use to determine datawidth
segsites = atoi(buf.c_str());
	
	int i;

	//cout << "Segsites: " << segsites << endl;
	
	
	


	while (buf != "positions:") {
		ciss >> buf;
		//cout << buf << endl;
	}
	
	
	
	string posline;
	getline(ciss, posline);
	
	
	istringstream pos(posline, istringstream::in);
	deque<string> postokens;
	string pospart;
	vector<double> positions;
	for (i=0; i<segsites; i++) {
		
		pos >> pospart;
		//cout << pospart << endl;
		positions.push_back(atof(pospart.c_str()));
		
	}
	
	//if this does not exhaust the line something is wrong
	if(!postokens.empty()){ //some mistake made in ms, no output
		cerr << "ERROR in runsim: Number of entires in positions line does not match the number of segregating sites" << endl;
		return -1;
		
	}
	
	

	
	
	//************************
	//number of sites may be different that world r!! base data needed off of received number of samples and sites
	
	datalength = tot;
	datawidth = 0;
	
	
	
	
	//************************
	vector<string> newas;
	newas.push_back("A");
	
	for(i=0;i<segsites;i++){
		
		//cout << i << endl;
		
		locinfo newlocinfo;
		newlocinfo.type = "SNP";
		newlocinfo.mutrate = mutrate; //FIXME div by 4N?
		newlocinfo.ancstate = newas;
		newlocinfo.recrate = recrate; //FIXME div by 4N?
		newlocinfo.phdist = basepairs/segsites; //FIXME - be more exact site to site?
		newlocinfo.length = 1;
		newlocinfo.offset = datawidth;
		datawidth += newlocinfo.length;
		
		
		loci.push_back(newlocinfo);
		
		
	}
	
	
	//cout << "datalength: " << datalength << endl;
	//cout << "datawidth: " << datawidth << endl;
	
	
	data = new char *[datalength] ;
	for(int x = 0 ; x < datalength ; x++ ) data[x] = new char[datawidth];
	
	
	//************************
	
	//cout << "Begin reading data" << endl;
	
	map<string, vector<int> >::const_iterator rdemes = r->dmap.begin();
	int di = 0;
	i = 0;
	nbp = nbypop.begin();
	ssp = tparams.samplesize.begin();
	while(nbp != nbypop.end()){ //iterator through the populations
		
		
		
		int dpop = *nbp; //THIS IS THE SAMPLE SIZE
		//cout << "Reading Pop: " << i << " which has " << dpop << " individuals." << endl;
		
		
		
		vector<int> newvec;
		if(rdemes ==r->dmap.end()) {
			cerr<< "Error! For macs deme #: " << i << " More demes that in real data." << endl;
			abort();
		}
		dmap[ssp->first] = newvec; //assign deme names back to sim demes
		
		
		
		
		
		
		for(int j = 0; j < dpop ; j++){//for each individual in the population
			string chrm;
			getline (ciss, chrm);
			//cout << j << ": ********" << endl << chrm << endl << endl;
			if(chrm.size() != segsites){
				cerr << "ERROR in rejector2 reading macs: Size of a chromosome's sites: " << chrm.size() << " unequal to declared segsites: " << segsites << endl;
				return -1;
			}
			
			//construct a name
			stringstream inn;
			if(genotypic){
				if(j%2==0) inn << i << "_" << j;
				else inn << i << "_" << (j-1);	
			}
			else inn << i << "_" << j;
			string inname = inn.str();
			//cout << "d00d#: " << inname << endl;
			
			
			
			
			int lwidth = 0;
			
			
			for(int k=0; k<chrm.size(); k++){
				string inlocus = (chrm.substr(k, 1));
				//cout << inlocus << " ";
				
				char thesnp;
				
				if(inlocus == "1")thesnp = 'T';
				else if(inlocus == "0")thesnp = 'A';
				else {
					cerr << "ERROR in rejstats -m: Do not recognize data: " << inlocus  << endl;
					return -1;
				}
				
				//check if length matches expected length
				if(loci[k].length != 1){
					cerr << "Error: For Tag " << inname << " length of locus! Expected " << loci[k].length << endl;
					
					abort();
				}
				
				lwidth++;
				if(lwidth > datawidth){
					cerr << "Error: For Tag " << inname << " length of locus " << (k+1) << " too large! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
					abort();
				}
				
				//cout << "Entering: " << thesnp << " at " << di << "-" << k << endl;
				data[di][k] = thesnp;
				
				
			}
			
			
			
			
			
			
			//cout << endl;
			
			nmap.push_back(inname);
			
			//put into dmap
			dmap[ssp->first].push_back(di);
			
			
			di++;
			//end j< dpop	
		}
		
		
		
		//end while nbp !=nbypop
		rdemes++;
		nbp++;
		ssp++;
		i++;
	}
	
	
	//cout << "End reading data" << endl;
	
	//sort loci?
	
	
	
	
	for(unsigned int lp=0; lp < loci.size(); lp++) {
		
		
		
		for(i = 0; i < datalength; i++){
			
			string entry = getdata(i,lp);
			//cout << entry << " ";
			
			if(entry[0] != 'X'){
				loci[lp].lmap.push_back(i);
			}
			
		}
		
		
		sort(loci[lp].lmap.begin(), loci[lp].lmap.end());
		//cout << endl;
		
		
	}
	
	
	
	
	//Sort the demes
	map<string, vector<int> >::iterator dp = dmap.begin();
	while(dp != dmap.end()){
		sort(dp->second.begin(), dp->second.end());
		
		dp++;
	}
	
	
return 0;	
}

//---------------------------------------------
//Function Name:
//Function Definition: Read ms output from stdin (piped)
int world::readmsinfile(world* r, string msline){
	double theta = 0.0;
	int segsites = 0;
	double initmig =0.0;
	vector<int> nbypop;
	map<int, double> popsizes;
	map<int, double> grrates;
	vector <double> inmutRates;
	vector <vector <string> > inancStates;
	string buf;
	
	
	
	//cout << "Begin readmsinfile" << endl;
	
	//Copy the loci info from r into the current world, as they must be the same
	
	
	vector<locinfo>::const_iterator lp = r->loci.begin();
	while(lp != r->loci.end()){
		
		locinfo newlocinfo;
		newlocinfo.type = lp->type;
		newlocinfo.mutrate = lp->mutrate;
		newlocinfo.ancstate = lp->ancstate;
		newlocinfo.recrate = lp->recrate;
		newlocinfo.phdist = lp->phdist;
		newlocinfo.length = lp->length;
		newlocinfo.offset = lp->offset;
		
		//BUT leave lmap until the data's been put in
		
		loci.push_back(newlocinfo);
		
		lp++;
	}
	
	
	//Base data size off of r data - Must be the same!
	datalength = r->datalength;
	datawidth = r->datawidth;
	
	
	data = new char *[datalength] ;
	for(int x = 0 ; x < datalength ; x++ ) data[x] = new char[datawidth];
	
	//ofstream wtf("msouts.txt", ios_base::app);
	
	
	//cout << "popen starting" << endl;
	//RUN MS
	char c;
	string msout;
	FILE *msp;
	msp = popen(msline.c_str(), "r");
	do {
		c = fgetc (msp);
		if (c != EOF) msout += c;
	} while (c != EOF);
	pclose(msp);
	
	
	//wtf << msout << endl << endl;
	
	//wtf.close();
	
	//cout << "popen complete" << endl;
	
	
	
	/*
	 string mscmdline;
	 getline (cin, mscmdline);
	 //cout << "CMDLINE:" << mscmdline << endl;
	 
	 istringstream ciss(mscmdline, istringstream::in);
	 
	 deque<string> mstokens;
	 string msclpart;
	 while( ciss >> msclpart )     
	 {
	 
	 //cout << msclpart << endl;
	 mstokens.push_back(msclpart);
	 
	 }
	 */
	
	
	
	//TOKENIZE THE STRINGSTREAM
	
	istringstream ciss(msout, istringstream::in);
	string mscmdline;
	getline(ciss, mscmdline);
	//cout << "CMDLINE:" << mscmdline << endl;
	
	istringstream mcmd(mscmdline, istringstream::in);
	deque<string> mstokens;
	string msclpart;
	while( mcmd >> msclpart )     
	{
		
		//cout << msclpart << endl;
		mstokens.push_back(msclpart);
		
	}
	
	
	if(mstokens.empty()){ //some mistake made in ms, no output
		cerr << "ERROR in runsim: No output from msHOT" << endl;
		return -1;
		
	}
	
	string mspart;
	
	//the command
	mspart = mstokens.front();
	//cout << mspart << endl;
	if(mspart != "./msHOT"){
		cerr << "ERROR in runsim: First word of read output from msHOt is not ./msHOT" << endl;
		return -1;
	}
	mstokens.pop_front();
	
	//Next is total number of samples
	//i++;
	int tot = atoi(mstokens.front().c_str());
	//cout << tot << endl;
	mstokens.pop_front();
	
	//Next is total number of samples. FOR NOW It MUST be 1
	int smpls = atoi(mstokens.front().c_str());
	//cout << smpls << endl;
	if(smpls != 1){
		cerr << "ERROR in runsim: Only accepts data from 1 sample from msHOt (i.e. ./msHOT 4 1 etc....)" << endl;
		return -1;
	}
	mstokens.pop_front();
	
	
	
	string msw;
	//Teh branchz0rs ms line reader starts here
	while (!mstokens.empty()){
		//Grab the first part. That will tell you what to do with the rest
		msw = mstokens.front();
		mstokens.pop_front();
		
		
		if(msw == "-t"){
			theta  = atof(mstokens.front().c_str());
			//cout << "Theta now: " << theta << endl;
			mstokens.pop_front();
		}
		else if (msw == "-s"){
			segsites = atoi(mstokens.front().c_str());
			mstokens.pop_front();
		}
		else if (msw == "-I"){
			int npop = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			for(int i =0; i<npop; i++){
				nbypop.push_back(atoi(mstokens.front().c_str()));
				
				mstokens.pop_front();
			}
			//check if migration rate added
			string migchk = mstokens.front();
			if(migchk.find_first_of("-") != string::npos){
				
			}
			else{
				initmig = atof(mstokens.front().c_str());
				mstokens.pop_front();
			}
			
		}
		else if (msw == "-n"){
			int npop = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			double popsz = atof(mstokens.front().c_str());
			mstokens.pop_front();
			popsizes[npop] = popsz;
		}
		
		else if (msw == "-g"){
			int npop = atoi(mstokens.front().c_str());
			mstokens.pop_front();
			double grrate = atof(mstokens.front().c_str());
			mstokens.pop_front();
			grrates[npop] = grrate;
		}
		
		else if (msw == "-G"){
			double grrate = atof(mstokens.front().c_str());
			mstokens.pop_front();
			int gpr = 0;
			map<int, double>::iterator gp = grrates.begin();
			while(gp != grrates.end()){
				grrates[gpr] = grrate;
				gp++;
				gpr++;
			}
				  }
				  
				  else if (msw == "-ej"){
					  
					  mstokens.pop_front();
					  
					  mstokens.pop_front();
					  mstokens.pop_front();
				  }
				  
				  else if (msw == "-eg"){
					  
					  mstokens.pop_front();
					  
					  mstokens.pop_front();
					  mstokens.pop_front();
				  }
				  
				  else if (msw == "-seeds"){
					  
					  mstokens.pop_front();
					  
					  mstokens.pop_front();
					  mstokens.pop_front();
				  }
				  
				  else if (msw == "-en"){ //pop size change
					  mstokens.pop_front();
					  
					  mstokens.pop_front();
					  mstokens.pop_front();
					  
				  }
		
		
				  else if (msw == "-r"){ //recrate
					  
					  mstokens.pop_front();
					  
					  mstokens.pop_front();			  
					  
				  }
		
				  else if (msw == "-v"){ //hotspots
					  int nhot = atoi(mstokens.front().c_str());
					  mstokens.pop_front();
					  for(int i =0; i<nhot; i++){
						  mstokens.pop_front();
						  mstokens.pop_front();
						  mstokens.pop_front();
					  }
					  
					  
				  }
				  
				  else {
					  
					  cerr << "ERROR in readmsline: Do not recognize command-line switrch " << msw << endl;
					  return -1;
					  break;
				  }
				  
				  //end while mstoken not empty
				  }
				  
				  
				  //cout << "Command line tokenized and read" << endl;
				  
				  
				  
				  //Assign all to one population if -I not used
				  if(nbypop.empty()) nbypop.push_back(tot);
				  
				  //Check that the sum of the pops adds to total
				  int totchk = 0;
				  vector<int>::const_iterator nbp = nbypop.begin();
				  while(nbp != nbypop.end()){
					  totchk += *nbp;
					  
					  nbp++;
				  }
				  if(totchk != tot){
					  cerr << "ERROR in rejstats :Total calculated from populations in -I switch not equal to nsam "  << endl;
					  return -1;
				  }
				  
				  //cout << "Totchk: " << totchk << endl;
				  
				  //CHECK THAT SIZES MATCH
				  //sample sizes
				  
				  map<string, int>::const_iterator ssp = tparams.samplesize.begin();
				  
				  if(tparams.samplesize.size() != nbypop.size()){
					  cerr << "Error in readmslime! Number of samples in input file: " << tparams.samplesize.size() << " but number referred to in priors: " << nbypop.size() << endl;
					  abort();
				  }
				  
				  
				  
				  string seeds;
				  getline (ciss, seeds);
				  //cout << "SEEDS:" << seeds << endl;
				  
				  
				  //doubleslash
				  ciss >> buf;
				  //cout << buf << endl;
				  if(buf != "//") return -1;
				  
				  
				  ciss >> buf; //segsites
				  if(buf != "segsites:") return -1; // Failed
				  
				  
				  //	cout << buf << endl;
				  ciss >> buf; 
				  //	cout << buf << endl;
				  segsites = atoi(buf.c_str());
				  
				  //cout << "segsites: " << segsites << endl;
				  
				  while(buf != "positions:")	ciss >> buf; 
				  string pos;
				  getline (ciss, pos);
				  //cout << "POS:" << pos << endl;
				  
				  
				  
				  //cout << "Begin reading data" << endl;
				  
				  map<string, vector<int> >::const_iterator rdemes = r->dmap.begin();
				  int di = 0;
				  int i = 0;
				  nbp = nbypop.begin();
				  ssp = tparams.samplesize.begin();
				  while(nbp != nbypop.end()){ //iterator through the populations
					  
					  
					  
					  int dpop = *nbp; //THIS IS THE SAMPLE SIZE
					  //cout << "Reading Pop: " << i << " which has " << dpop << " individuals." << endl;
					  
					  
					  
					  vector<int> newvec;
					  if(rdemes ==r->dmap.end()) {
						  cerr<< "Error! For ms deme #: " << i << " More demes that in real data." << endl;
						  abort();
					  }
					  dmap[ssp->first] = newvec; //assign deme names back to sim demes
					  
					  
					  
					  
					  
					  
					  for(int j = 0; j < dpop ; j++){//for each individual in the population
						  string chrm;
						  getline (ciss, chrm);
						  if(chrm.size() != segsites){
							  cerr << "ERROR in rejstats -m: Size of a chromosome's sites: " << chrm.size() << " unequal to declared segsites: " << segsites << endl;
							  return -1;
						  }
						  
						  //construct a name
						  stringstream inn;
						  if(genotypic){
							  if(j%2==0) inn << i << "_" << j;
							  else inn << i << "_" << (j-1);	
						  }
						  else inn << i << "_" << j;
						  string inname = inn.str();
						  //cout << "d00d#: " << inname << endl;
						  
						  
						  
						  
						  int lwidth = 0;
						  
						  
						  for(int k=0; k<chrm.size(); k++){
							  string inlocus = (chrm.substr(k, 1));
							  //cout << inlocus << " ";
							  
							  char thesnp;
							  
							  if(inlocus == "1")thesnp = 'T';
							  else if(inlocus == "0")thesnp = 'A';
							  else {
								  cerr << "ERROR in rejstats -m: Do not recognize data: " << inlocus  << endl;
								  return -1;
							  }
							  
							  //check if length matches expected length
							  if(loci[k].length != 1){
								  cerr << "Error: For Tag " << inname << " length of locus! Expected " << loci[k].length << endl;
								  
								  abort();
							  }
							  
							  lwidth++;
							  if(lwidth > datawidth){
								  cerr << "Error: For Tag " << inname << " length of locus " << (k+1) << " too large! Total length so far in line: " << lwidth << " Datawidth is:  " << datawidth << endl;
								  abort();
							  }
							  
							  //cout << "Entering: " << thesnp << " at " << di << "-" << k << endl;
							  data[di][k] = thesnp;
							  
							  
						  }
						  
						  
						  
						  
						  
						  
						  //cout << endl;
						  
						  nmap.push_back(inname);
						  
						  //put into dmap
						  dmap[ssp->first].push_back(di);
						  
						  
						  di++;
						  //end j< dpop	
					  }
					  
					  
					  
					  //end while nbp !=nbypop
					  rdemes++;
					  nbp++;
					  ssp++;
					  i++;
				  }
				  
				  
				  //cout << "End reading data" << endl;
				  
				  //sort loci?
				  
				  
				  
				  
				  for(unsigned int lp=0; lp < loci.size(); lp++) {
					  
					  
					  
					  for(i = 0; i < datalength; i++){
						  
						  string entry = getdata(i,lp);
						  //cout << entry << " ";
						  
						  if(entry[0] != 'X'){
							  loci[lp].lmap.push_back(i);
						  }
						  
					  }
					  
					  
					  sort(loci[lp].lmap.begin(), loci[lp].lmap.end());
					  //cout << endl;
					  
					  
				  }
				  
				  
				  
				  
				  //Sort the demes
				  map<string, vector<int> >::iterator dp = dmap.begin();
				  while(dp != dmap.end()){
					  sort(dp->second.begin(), dp->second.end());
					  
					  dp++;
				  }
				  
				  
				  
				  
				  
			/*
			 cin >> buf;
			 for(i = 0;i < datalength ; i++){
			 
			 cout << i << " ";
			 //need offset here
			 
			 for(int j=0; j < datawidth; j++){
			 
			 
			 cout << data[i][j];
			 
			 
			 }
			 cout << endl;
			 }
			 */
				  
				  //cout << "End readmsinfile" << endl;
				  
				  return 0;
}
				  
				  
world runsim(world *r, priordist *pt, bool poplock, string basefile)//this will be where prior does its thing
			{
				int syschk = 0;
				world s;
				
				
				if(simtype == 0) { //0 means simcoal2
					//CREATE SC2 .par FILE
					string poutfile = basefile + ".par";
					ofstream pout(poutfile.c_str());
					if(!pout) {
						cerr << "Error! Cannot create file " << poutfile << " for simcoal2. Loop failed." << endl;
						s.broken = true;
						return s;
					}
					
					//s.tparams = makescissfile(r, pt, pout, poplock);
					syschk = s.makescinfile(r, pt, pout, poplock);
					
					pout.close();
					
					//RUN SIMCOAL
					string simcoalcmd = "./simcoalrej2_1_2 ";
					simcoalcmd += poutfile;
					if(genotypic) simcoalcmd += " 1 1";
					else simcoalcmd += " 1 0"; //one run, haplotypic data
					
					//cout << simcoalcmd << endl;
					
					syschk = system (simcoalcmd.c_str());
					if (syschk == -1){
						cerr << "Error executing simcoal!";
						s.broken = true;
						return s;
					}
					
					
					string soutfile =  basefile + "_0.arp";
					
					//load simcoal output
					ifstream sin(soutfile.c_str());
					if(!sin) {
						cerr << "Cannot find " << soutfile << " generated by simcoal!" << endl;
						s.broken = true;
						return s;
					}
					
					s.readarpfile(sin, r);
					
					
					
					
					//kill simcoal files
					killscfiles(basefile);
					
					
					s.broken = false;
					return s;
				}
				
				else if(simtype == 1){ // 1 means msHOT
					//string msline =  makemsline(r, pt, poplock);
					string msline =  s.makemsline(r, pt, poplock, false);
					
//                    cout << "msline: " << msline << endl;
					
					
					
					
					if(msline == "error") {
						cerr << "Error! Bad string for ms. Loop failed." << endl;
						s.broken = true;
						return s;
					}
					
					syschk = s.readmsinfile(r, msline);
					if (syschk == -1){
						cerr << "Error executing ms!";
						s.broken = true;
						return s;
					}
					
					
					
					s.broken = false;
					return s;
				}
				
				else if(simtype == 2){ // 2 means MaCS
					string msline =  s.makemsline(r, pt, poplock, true);
					
					if(msline == "error") {
						cerr << "Error! Bad string for macs. Loop failed." << endl;
						s.broken = true;
						return s;
					}
					
					//cout << "macsline: " << msline << endl;
					
					//syschk = s.readmacsinfile(r, msline);
					syschk = s.readmacsunformatted(r, msline);
					if (syschk == -1){
						cerr << "Error executing macs!";
						s.broken = true;
						return s;
					}
					
					
					
					s.broken = false;
					return s;
					
					
				}
				
				
				//If it gets here there's troubel
				s.broken = true;
				return s;
			}
				  
				  
	  //---------------------------------------------
	  //Function Name:
	  //Function Definition: Get rid of the  prior file in case of a problem with simcoal
	  
	  void killscfiles(string basefile){
		  int syschk;
		  
		  string parname, arpname, arpoutname, arbname, simparamname, snpname;
		  
		  parname = basefile + ".par";
		  arpname = basefile + "_0.arp";
		  arbname = basefile + ".arb";
		  simparamname = basefile + ".simparam";
		  arpoutname = basefile + "_0.out";
		  snpname = basefile + "-snp.txt";
		  
		  
		  
		  
		  string killfile = "rm -f ";
		  
		  killfile += parname;
		  
		  
		  
		  syschk = system (killfile.c_str());
		  if (syschk==-1){
			  cerr << "Error executing " << killfile << " !";
		  }
		  
		  
		  
		  killfile = "rm -f ";
		  
		  killfile += arpname;
		  
		  
		  
		  syschk = system (killfile.c_str());
		  if (syschk==-1){
			  cerr << "Error executing " << killfile << " !";
		  }
		  
		  
		  
		  
		  
		  killfile = "rm -f ";
		  
		  killfile += arpoutname;
		  
		  
		  
		  syschk = system (killfile.c_str());
		  if (syschk==-1){
			  cerr << "Error executing " << killfile << " !";
		  }
		  
		  killfile = "rm -f ";
		  
		  killfile += arbname;
		  
		  
		  
		  syschk = system (killfile.c_str());
		  if (syschk==-1){
			  cerr << "Error executing " << killfile << " !";
		  }
		  
		  
		  killfile = "rm -f ";
		  
		  killfile += simparamname;
		  
		  
		  
		  syschk = system (killfile.c_str());
		  if (syschk==-1){
			  cerr << "Error executing " << killfile << " !";
		  }
		  
		  
		  killfile = "rm -f ";
		  
		  killfile += snpname;
		  
		  
		  
		  syschk = system (killfile.c_str());
		  if (syschk==-1){
			  cerr << "Error executing " << killfile << " !";
		  }
		  
		  
	  }
				  
				  

//---------------------------------------------
//Function Name:
//Function Definition: outputs world in rejstats input file format
void world::outputinfile(ostream &out, priordist pt){
	
	
	out << "REJECTOR" << endl << endl;
	
	out << "Number of Populations" << endl;
	out << dmap.size() << endl;
	
	
	
	out << endl;
    
    
    
    out << "/--Priors" << endl;
    out << "Population Sizes" << endl;
    vector<distn>::const_iterator ps = pt.popsizes.begin();
    while(ps != pt.popsizes.end()){
        out << ps->name << endl;
        out << ps->v1 << endl;
        out << ps->v2 << endl;
        ps++;
    }
    out << endl;
    
    out << "Growth Rates" << endl;
    ps = pt.growthrates.begin();
    while(ps != pt.growthrates.end()){
        out << ps->name << endl;
        out << ps->v1 << endl;
        out << ps->v2 << endl;
        ps++;
    }
    out << endl;
    
    
    out << "Number of Migration Matrices" << endl;
    out << "0" << endl;
    out << endl;
    
    
    out << "Number of Historical Events" << endl;
    out << "0" << endl;
    out << endl;
    
    out << "Microsat Range Constraint" << endl;
    out << pt.msrc.name << endl;
    out << pt.msrc.v1 << endl;
    out << pt.msrc.v2 << endl;
    out << endl;
    
    
    
    
	
	out << "/--Statlist" << endl;
	
	map<string, double>::const_iterator tl = pt.alphalist.begin();
	
	while(tl != pt.alphalist.end()){
		out << tl->first << "\t" << tl->second << endl;
		tl++;
	}
	
	out << endl;
	
	out << "/--Data" << endl;
	out << "Vnaught\t" << Vnaught << endl << endl;
	
	vector<locinfo>::const_iterator lp = loci.begin();
	out << "Loci" << "\t";
	while(lp != loci.end()){
		out << lp->type << "\t";
		lp++;
	}
	out << endl;
	
	lp = loci.begin();
	out << "MutRates" << "\t";
	while(lp != loci.end()){
		out << lp->mutrate << "\t";
		lp++;
	}
	out << endl;
	
	lp = loci.begin();
	out << "Ancestral" << "\t";
	while(lp != loci.end()){
		vector<string>::const_iterator llp = lp->ancstate.begin();
		if(lp->ancstate.size() >1 ){
			out << "\"";
			while(llp != lp->ancstate.end()){
				out << *llp;
				if(llp != lp->ancstate.end()) out << ",";
				llp++;
			}
			out << "\"";
		}
		else out << *llp;		  
		
		
		out << "\t";
		lp++;
	}
	out << endl;
	
	lp = loci.begin();
	out << "RecombRt" << "\t";
	while(lp != loci.end()){
		out << lp->recrate << "\t";
		lp++;
	}
	out << endl;
	
	
	lp = loci.begin();
	out << "NumLoci" << "\t";
	while(lp != loci.end()){
		out << "1\t";
		lp++;
	}
	out << endl;
	
	lp = loci.begin();
	out << "Length" << "\t";
	while(lp != loci.end()){
		out << lp->length << "\t";
		lp++;
	}
	out << endl << endl;
	
	
	
	out << "Tag\tPopulation" << endl;
	
	
	
    map<string, vector<int> >::const_iterator dp  = dmap.begin();
	while(dp != dmap.end()){
        
        vector<int>::const_iterator dpp = dp->second.begin();
        while(dpp != dp->second.end()){
            out << nmap[*dpp] << " " << dp->first << " ";
           // cout << *dpp << " " << nmap[*dpp] << dp->first << " ";
            
            for(unsigned int lp=0; lp < loci.size(); lp++) {
                
                for(unsigned int lpp=loci[lp].offset; lpp < (loci[lp].offset+loci[lp].length); lpp++){ 
                  //  cout << data[*dpp][lpp];
                    out << data[*dpp][lpp];
                    
                }
              //  cout << " ";
                out << " ";
                
            }
           // cout << endl;
            out << endl;
            
            dpp++;
        }
        
        dp++;
        
    }
	
	
	
}
				  
				  
				  

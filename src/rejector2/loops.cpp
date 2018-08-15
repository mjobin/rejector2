//------------------------------------------------
//
// REJSTATS
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropology
// Stanford University, 2007
//
//------------------------------------------------



#include <cmath>
#include <climits>
#include <sstream>
#include "rdeme.h"

using namespace std;

void world::printstatresult(statresult* sr, ostream &out)
{

	//cout << sr->type << endl;
	if(sr->type == "PerSystem") persystemprint(sr,out);
	if(sr->type == "PerLocusXPop") perlocusxpopprint(sr, out);
	if(sr->type == "PerLocusXPopAv") perlocusxpopavprint(sr, out);
	if(sr->type == "Oneway") onewayprint(sr, out);
	if(sr->type == "OnewayAv") onewayavprint(sr, out);
	if(sr->type == "ModOnewayAv") modonewayavprint(sr, out);
	if(sr->type == "Twoway") twowayprint(sr, out);
	if(sr->type == "TwowayXloci") twowayxlociprint(sr, out);
	if(sr->type == "TwoLoci") twolociprint(sr, out);
	if(sr->type == "TwoLociAdj") twolociadjprint(sr, out);
	if(sr->type == "ByDeme") bydemeprint(sr,out);
	if(sr->type == "ByTwoDeme") bytwodemeprint(sr,out);
	
	

}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs persystem tests
statresult world::persystemloop (string locustype, string ancdec, double (world::*tc) (vector<int>*)){
	statresult x;
	x.type = "PerSystem";
	
	map<string, vector<int> >::iterator di = dmap.begin();
	map<string, map<string, vector<double> > > in2dememap;
	while(di != dmap.end()){ //iterate through demes
	


        map<string, vector<double> > inSys;



				
				
				vector<double> inVec; // will always have precisely one entry in persystem loops
				
				if(di->second.size() >= DEME_CUTOFF){
					
					vector<int> y;
					
					//Get the deme vector
					
					//	do a set—intersection on the deme vector and the locus iterator to get the y 
					//they are now both sorted
					
					
					if(ancdec == "ANC"){
						y = getancdecvec(&di->second, true);
						
					}
					else if(ancdec == "DEC"){
						y = getancdecvec(&di->second, false);
						

						
					}

					else{
						y = di->second;
											
					}
					

					if(y.empty()) inVec.push_back(0); //for persystems a empty deme means a zero, not an NA
					else inVec.push_back((*this.*tc)(&y));
					
                }
				else{
					inVec.push_back(-999);
				}
				
				
				
				inSys["PerSystem"] = inVec;


        in2dememap[di->first]=inSys;
        di++;
    }
	x.result["Persystem"]=in2dememap;
	
	return x;
}

//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs tests by passing name of deme
statresult world::bydemeloop (string locustype, string ancdec, double (world::*tc) (string)){
	statresult x;
	x.type = "ByDeme";
	
	map<string, vector<int> >::iterator di = dmap.begin();
	map<string, map<string, vector<double> > > in2dememap;
	while(di != dmap.end()){ //iterate through demes
		
		
		
        map<string, vector<double> > inSys;
		
		
		
		
		
		vector<double> inVec; // will always have precisely one entry in persystem loops
		
		if(di->second.size() >= DEME_CUTOFF){
			
			inVec.push_back((*this.*tc)(di->first));
			
		}
		else{
			inVec.push_back(-999);
		}
		
		
		
		inSys["ByDeme"] = inVec;
		
		
        in2dememap[di->first]=inSys;
        di++;
    }
	x.result["ByDeme"]=in2dememap;
	
	return x;
}


/*
 
 while(di != (--dmap.end())){   //WATCH THIS BIT
 //cout << "Deme:" << di->second.name << endl;
 map<string, vector<int> >::const_iterator dj = di;
 dj++;
 map<string, map<string, vector<double> > > in2dememap;
 while(dj != dmap.end()){
 //cout << "Versus:" << dj->second.name << endl;
 
 map<string, vector<double> > in2deme;
*/ 

//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs tests by passing name of deme
statresult world::bytwodemeloop (string locustype, string ancdec, double (world::*tc) (string, string)){
	statresult x;
	x.type = "ByTwoDeme";
	
	map<string, vector<int> >::iterator di = dmap.begin();
	
	if(dmap.size() < 2){//
		map<string, map<string, vector<double> > > in2dememap;
		map<string, vector<double> > in2deme;
		vector<double> inVec;
		inVec.push_back(-999);
		in2deme["ByTwoDeme"] = inVec;
		in2dememap["ByTwoDeme"]=in2deme;
		x.result[di->first]=in2dememap;
		
		return x;
	}
	
	
	while(di != (--dmap.end())){   //WATCH THIS BIT
		map<string, vector<int> >::const_iterator dj = di;
		dj++;
		map<string, map<string, vector<double> > > in2dememap;
		while(dj != dmap.end()){ //iterate through demes
		
		
		
        map<string, vector<double> > inSys;
		
		
		
		
		
		vector<double> inVec; // will always have precisely one entry in persystem loops
		
		if(di->second.size() >= DEME_CUTOFF){
			
			inVec.push_back((*this.*tc)(di->first, dj->first));
			
		}
		else{
			inVec.push_back(-999);
		}
		
		
		
		inSys["ByTwoDeme"] = inVec;
		
		
        in2dememap[dj->first]=inSys;
			dj++;
		}
		
		//end while di
		x.result[di->first]=in2dememap;
		di++;
	}

	
	return x;
}




//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs perlocusxpop tests
statresult world::perlocusxpoploop (string locustype, string ancdec, double (world::*tc) (int)){
	statresult x;
	x.type = "PerLocusXPop";
	
	
	map<string, map<string, vector<double> > > in2dememap;
	map<string, vector<double> > indeme;
	vector<double> inVec;
	
	for(unsigned int j =0; j<loci.size(); j++){
		
		if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
			inVec.push_back((*this.*tc)(j));
		}
		else inVec.push_back(-999);
	}
	indeme["PerLocusXPop"] = inVec;
	in2dememap["PerLocusXPop"] = indeme;
	
	
	
	x.result["PerLocusXPop"]=in2dememap;
	
	
	return x;
}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs perlocusxpop tests
statresult world::perlocusxpopavloop (string locustype, string ancdec, double (world::*tc) (int)){
	statresult x;
	x.type = "PerLocusXPopAv";
	
	
	map<string, map<string, vector<double> > > in2dememap;
	map<string, vector<double> > indeme;
	vector<double> inVec;
	
	double tot, count = 0.0;
	
	for(unsigned int j =0; j<loci.size(); j++){
		
		if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
			//inVec.push_back((*this.*tc)(j));
			tot += ((*this.*tc)(j));
			count++;
		}
			//	else inVec.push_back(-999);
	}
	
	if(count > 0.0){ //i.e. at least one result was not -999 or nan
		inVec.push_back(tot/count);
	}
	else {
		inVec.push_back(-999);
	}
	
	
	/*
	double tot = 0.0;
	
	vector<double>::const_iterator av = inVec.begin();
	while(av != inVec.end()){
		tot += *av;
		
		av++;
	}
	
	vector<double> avfromvec;
	avfromvec.push_back(tot/(double)inVec.size()); //only one entry
	*/
	indeme["PerLocusXPopAv"] = inVec;
	//indeme["PerLocusXPop"] = avfromvec;
	in2dememap["PerLocusXPopAv"] = indeme;
	


	
	
	/*
	//map<string, map<string, vector<double> > > testsbysys = ((*this.*tc)(&dmap));
	
	double syscount = 0.0;
	double systot = 0.0;
	
	
	map<string, map<string, vector<double> > >::const_iterator tp = testsbysys.begin();
	while(tp != testsbysys.end()){
		map<string, vector<double> >::const_iterator tq = tp->second.begin();
		while(tq != tp->second.end()){
			vector<double>::const_iterator tr = tq->second.begin();
			while(tr != tq->second.end()){
				double tnum = *tr;
				if(isnan(tnum) || tnum == -999){
					tr++;
				}
				else{
					systot += *tr;
					tr++;
					syscount++;
				}
			}
			
			tq++;
		}
		tp++;
	}
	
	double avfst = systot/syscount;
	
	vector<double> insys;
	insys.push_back(avfst);
	
	map<string, vector<double> > indeme;
	
	indeme["ALL"] = insys;
	
	map<string, map<string, vector<double> > > in2dememap;
	
	in2dememap["ALL"] = indeme;
	*/
	
	x.result["PerLocusXPopAv"]=in2dememap;
	
	
	return x;
	
	
	
}




//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs oneway tests
statresult world::onewayloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, int)){
	statresult x;
	x.type = "Oneway";
	map<string, vector<int> >::const_iterator di = dmap.begin();
	map<string, map<string, vector<double> > > in2dememap;
	while(di != dmap.end()){ //iterate through demes
		
		//int j=0;
		map<string, vector<double> > inSys;
		
		vector<double> inVec;
		for(unsigned int j =0; j<loci.size(); j++){
		//vector<locinfo*>::const_iterator li = loci.begin();
		//while(li != loci.end()){
			
			if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
				//pair<multimap<string, int>::const_iterator, multimap<string, int>::const_iterator> bounds;
				//bounds = dmap.equal_range(di->first);

				//cout << "Deme: " << di->first << " Locus: " << j << " Type: " << loci[j].type << endl;
				vector<int> y;
								
				vector<int>::iterator lp, lq;
				unsigned int vecsz = 0;
				
				 if(ancdec == "ANC"){
					 vector<int> adlist = getancdecvec(&loci[j].lmap, true);
					 lp = adlist.begin();
					 lq = adlist.end();
					 vecsz = adlist.size();
				 
				 }
				 else if(ancdec == "DEC"){
				 vector<int> adlist = getancdecvec(&loci[j].lmap, true);
					 lp = adlist.begin();
					 lq = adlist.end();
					  vecsz = adlist.size();
				 
				 }
				 //get a vec of all the ints that match
				 else{
				
					 lp = loci[j].lmap.begin();
					 lq = loci[j].lmap.end();
					 vecsz = loci[j].lmap.size();
				 
				 }
				
				if(vecsz >=  DEME_CUTOFF){
				
					set_intersection(di->second.begin(), di->second.end(), lp, lq, std::back_inserter(y));
					
					if(y.size() >=DEME_CUTOFF) {
						inVec.push_back((*this.*tc)(&y, j));
					}
					else inVec.push_back(-999); //if fails cutoff
					
					
				}
				
				else inVec.push_back(-999); 


					

				
				//NOW RUN THE TEST ON THE llist if there are still enough to work on

			
			}
			
			else inVec.push_back(-999); //if fails cutoff or wrong type
			

		}
		inSys["Oneway"] = inVec;
		
		in2dememap[di->first]=inSys;
		di = dmap.upper_bound(di->first); //advance to next allele
	}
	x.result["Oneway"] = in2dememap;
	return x;
	
	/*
	 if(ancdec == "ANC"){
	 vector<int> anclist = getancdecvec(&loci[j].lmap, true);
	 set_intersection(di->second.begin(), di->second.end(), anclist.begin(), anclist.end(), std::back_inserter(y));
	 
	 }
	 else if(ancdec == "DEC"){
	 vector<int> declist = getancdecvec(&loci[j].lmap, true);
	 set_intersection(di->second.begin(), di->second.end(), declist.begin(), declist.end(), std::back_inserter(y));
	 
	 }
	 //get a vec of all the ints that match
	 else{
	 set_intersection(di->second.begin(), di->second.end(), loci[j].lmap.begin(), loci[j].lmap.end(), std::back_inserter(y));
	 
	 
	 }
	 */
	
	/*
	 cout << "----------" << endl;
	 cout << "Oneway Loop for Deme: " << di->first << " and Locus: " << j << endl;
	 vector<int>::const_iterator py = y.begin();
	 while(py != y.end()){
	 cout << *py << " ";
	 py++;
	 }
	 cout << endl << endl;
	 */
	
}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs oneway tests
statresult world::onewayavloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, int)){
	
	statresult x;
	x.type = "OnewayAv";
	map<string, vector<int> >::const_iterator di = dmap.begin();
	map<string, map<string, vector<double> > > in2dememap;
	while(di != dmap.end()){ //iterate through demes
		
		//int j=0;
		//map<string, vector<double> > inSys; //Defined by type here, SNP, STR et5c.
		map<string, loopav > inSys; //Defined by type here, SNP, STR et5c.
		
		
		
		
		for(unsigned int j =0; j<loci.size(); j++){
			//vector<locinfo*>::const_iterator li = loci.begin();
			//while(li != loci.end()){
			
			if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 

				
				//cout << "Deme: " << di->first << " Locus: " << j << " Type: " << loci[j].type << endl;
				vector<int> y;
				
				vector<int>::iterator lp, lq;
				unsigned int vecsz = 0;
				
				if(ancdec == "ANC"){
					vector<int> adlist = getancdecvec(&loci[j].lmap, true);
					lp = adlist.begin();
					lq = adlist.end();
					vecsz = adlist.size();
					
				}
				else if(ancdec == "DEC"){
					vector<int> adlist = getancdecvec(&loci[j].lmap, true);
					lp = adlist.begin();
					lq = adlist.end();
					vecsz = adlist.size();
					
				}
				//get a vec of all the ints that match
				else{
					
					lp = loci[j].lmap.begin();
					lq = loci[j].lmap.end();
					vecsz = loci[j].lmap.size();
					
				}
				
				if(vecsz >=  DEME_CUTOFF){
					
					set_intersection(di->second.begin(), di->second.end(), lp, lq, std::back_inserter(y));
					
					if(y.size() >=DEME_CUTOFF) {
						
						if(inSys.find(loci[j].type) == inSys.end()){ //put in a new vector is this type not yet used
							loopav newloopav;
							newloopav.tot = 0.0;
							newloopav.count = 0.0;
							inSys[loci[j].type] = newloopav;
						}
						//Add to type vector
						double theresult = ((*this.*tc)(&y, j));
						//inSys[loci[j].type].second.push_back((*this.*tc)(&y, j));
						if(theresult != -999){
							inSys[loci[j].type].tot += theresult;
							inSys[loci[j].type].count++;
							
						}
					
						
						//inVec.push_back((*this.*tc)(&y, j));
						//tot += ((*this.*tc)(&y, j));
						//count++;
					}
					//else inVec.push_back(-999); //if fails cutoff
					
					
				}
				
				//else inVec.push_back(-999); 
				
				
				
				
				
				//NOW RUN THE TEST ON THE llist if there are still enough to work on
				
				
			}
			
			//else inVec.push_back(-999); //if fails cutoff or wrong type
			
			
		}
		
		//check if whole types empty, and insert a dummy if so
		for(unsigned int j =0; j<loci.size(); j++){
			if(inSys.find(loci[j].type) == inSys.end()){
				loopav newloopav;
				newloopav.tot = -999;
				newloopav.count = 0.0;
				inSys[loci[j].type] = newloopav;
				
			}
			
		}
		
		
		//Average through
		map<string, vector<double> > avinSys;
		map<string, loopav>::const_iterator isp = inSys.begin();
		while(isp != inSys.end()){
			if(isp->second.count > 0.0){
				double vecav = isp->second.tot/isp->second.count;
				vector<double> avinVec;
				avinVec.push_back(vecav);
				avinSys[isp->first] = avinVec;
				
			}
			
			else {

				vector<double> avinVec;
				avinVec.push_back(-999);
				avinSys[isp->first] = avinVec;
				
			}
			
			
			
			isp++;
		}
		
		/*
		if(count > 0.0){ //i.e. at least one result was not -999 or nan
			inVec.push_back(tot/count);
		}
		else {
			inVec.push_back(-999);
		}
		 */
		
		
		
		//inSys["OnewayAv"] = inVec;
		
		in2dememap[di->first]=avinSys;
		di = dmap.upper_bound(di->first); //advance to next allele
	}
	x.result["OnewayAv"] = in2dememap;
	return x;
	
	/*
	
	statresult x;
	x.type = "OnewayAv";
	map<string, vector<int> >::const_iterator di = dmap.begin();
	map<string, map<string, vector<double> > > in2dememap;
	

    while(di != dmap.end()){ //iterate through each deme

        map<string, vector<double> > inSys; //here defines firsts as "SNP, STR..."
		

		
		vector<double> inVec;
		
		for(unsigned int j =0; j<loci.size(); j++){
		
		

					if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
						
						vector<int> y;
						
						vector<int>::iterator lp, lq;
						unsigned int vecsz = 0;
						
						if(ancdec == "ANC"){
							vector<int> adlist = getancdecvec(&loci[j].lmap, true);
							lp = adlist.begin();
							lq = adlist.end();
							vecsz = adlist.size();
							
						}
						else if(ancdec == "DEC"){
							vector<int> adlist = getancdecvec(&loci[j].lmap, true);
							lp = adlist.begin();
							lq = adlist.end();
							vecsz = adlist.size();
							
						}
						//get a vec of all the ints that match
						else{
							
							lp = loci[j].lmap.begin();
							lq = loci[j].lmap.end();
							vecsz = loci[j].lmap.size();
							
						}
						
						if(vecsz >=  DEME_CUTOFF){
							
							set_intersection(di->second.begin(), di->second.end(), lp, lq, std::back_inserter(y));
							
							if(y.size() >=DEME_CUTOFF) {
								inVec.push_back((*this.*tc)(&y, j));
							}
							else inVec.push_back(-999); //if fails cutoff
							
							
						}
						
						else inVec.push_back(-999); 
						//inSys["OnewayAv"] = inVec;
						

						
						
					}
					else{

						inVec.push_back(-999);
						//inSys["OnewayAv"] = inVec;
					}
				}
				
		inSys["OnewayAv"] = inVec;

		//AVERAGING

		
		//average through, but ignore the -999's !!
		map<string, vector<double> > avinSys;
		
		map<string, vector<double> >::const_iterator isp = inSys.begin();
		while(isp != inSys.end()){
			//cout << "sys: " << isp->first << endl;
			unsigned int ignore = 0;
			double vectot = 0.0;
			vector<double>::const_iterator vp = isp->second.begin();
			while (vp != isp->second.end()){
				//cout << *vp << endl;
				if(*vp == -999 || isnan(*vp)) ignore++;
				//if(*vp == -999) ignore++;
				else vectot += *vp;
				vp++;
			}
			
			unsigned int denom = isp->second.size() - ignore;
			double vecav;
			if (denom > 0) vecav = vectot / ((double) denom);
			else vecav = -999;
			
			vector<double> avinVec;
			avinVec.push_back(vecav);
			
			avinSys[isp->first] = avinVec;
			
			isp++;
		}
		
		
		
		
		in2dememap[di->first]=avinSys;
        di++;
    }
	
	x.result["OnewayAv"] = in2dememap;
	return x;
	*/
	
}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs oneway tests
statresult world::modonewayavloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, int)){
	
	statresult x;
	x.type = "ModOnewayAv";
	map<string, vector<int> >::const_iterator di = dmap.begin();
	
	map<string, map<string, map<string, vector<double> > > > preresults; 
	
	map<string, map<string, vector<double> > > prein2dememap;
    while(di != dmap.end()){ //iterate through each deme
		// cout << "Deme: " << di->first << endl;

        map<string, vector<double> > inSys; //here defines firsts as "SNP, STR..."
		
		vector<double> inVec;
		
		for(unsigned int j =0; j<loci.size(); j++){
		


				


					
					if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
						
						vector<int> y;
						
						vector<int>::iterator lp, lq;
						unsigned int vecsz = 0;
						
						if(ancdec == "ANC"){
							vector<int> adlist = getancdecvec(&loci[j].lmap, true);
							lp = adlist.begin();
							lq = adlist.end();
							vecsz = adlist.size();
							
						}
						else if(ancdec == "DEC"){
							vector<int> adlist = getancdecvec(&loci[j].lmap, true);
							lp = adlist.begin();
							lq = adlist.end();
							vecsz = adlist.size();
							
						}
						//get a vec of all the ints that match
						else{
							
							lp = loci[j].lmap.begin();
							lq = loci[j].lmap.end();
							vecsz = loci[j].lmap.size();
							
						}
						
						if(vecsz >=  DEME_CUTOFF){
							
							set_intersection(di->second.begin(), di->second.end(), lp, lq, std::back_inserter(y));
							
							if(y.size() >=DEME_CUTOFF) {
								inVec.push_back((*this.*tc)(&y, j));
							}
							else inVec.push_back(-999); //if fails cutoff
							
							
						}
						
						else inVec.push_back(-999); 

						
						/*
						vector<int> y;
						
						if(ancdec == "ANC"){
							vector<int> anclist = getancdecvec(&loci[j].lmap, true);
							set_intersection(di->second.begin(), di->second.end(), anclist.begin(), anclist.end(), std::back_inserter(y));
							
						}
						else if(ancdec == "DEC"){
							vector<int> declist = getancdecvec(&loci[j].lmap, true);
							set_intersection(di->second.begin(), di->second.end(), declist.begin(), declist.end(), std::back_inserter(y));
							
						}
						else{
							set_intersection(di->second.begin(), di->second.end(), loci[j].lmap.begin(), loci[j].lmap.end(), std::back_inserter(y));
							
							
						}
						
						
						//NOW RUN THE TEST ON THE llist if there are still enough to work on
						if(y.size() >=DEME_CUTOFF) {
							inVec.push_back((*this.*tc)(&y, j));
						}
						else inVec.push_back(-999); //if fails cutoff
						*/
						
						
						
						inSys["ModOnewayAv"] = inVec;
						
					}
					else{
						//vector<double> inVec = inSys[sl->second[j]];
						inVec.push_back(-999);
						inSys["ModOnewayAv"] = inVec;
					}
					
					
					
					
					
					//end for j loop
				}
				
				

			

        
		
		
		
		
		prein2dememap[di->first]=inSys;
        di++;
    }
	preresults["ModOnewayAv"]=prein2dememap;
	
	//cout << "preresults" << endl;
	//Go through and check if for any deme as a locus there is an NA, make all corresponding loci in other demes also NA
	
	map<string, map<string, map<string, vector<double> > > >::iterator rp = preresults.begin(); //oneway
	map<string, map<string, vector<double> > >::iterator rq = rp->second.begin();
	//cout << "Deme\t" << rq->first << endl; //first deme only here
	map<string, vector<double> >::iterator rr = rq->second.begin();
	while(rr != rq->second.end()){
		//cout << "System\t" << rr->first << endl;
		int i = 0;
		vector<double>::iterator rs = rr->second.begin();
		while(rs != rr->second.end()){
			//cout << "Locus\t"<< i << endl << endl;
			bool anynas = false;
			map<string, map<string, vector<double> > >::iterator qq = rp->second.begin();
			while(qq != rp->second.end()){
				//cout << "Compare Deme\t" << qq->first << endl; 
				map<string, vector<double> >::iterator qr = qq->second.find(rr->first);
				if(qr != qq->second.end()){
					
					//cout << "Compare System\t" << qr->first << endl; 
					vector<double>::iterator qs = qr->second.begin();
					//cout << "Compare Locus\t" << i << endl;
					for(int j = 0; j<i; j++) qs++;
					//cout << "Compare\t" << *qs << endl;
					if(isnan(*qs) || *qs == -999){
						anynas = true;
						
						break;
					}
					
					
				}
				qq++;
			}
			
			//now if anynas is true -999 to all!
			if(anynas){
				qq = rp->second.begin();
				while(qq != rp->second.end()){
					//cout << "Change Deme\t" << qq->first << endl; 
					map<string, vector<double> >::iterator qr = qq->second.find(rr->first);
					if(qr != qq->second.end()){
						//cout << "Change System\t" << qr->first << endl; 
						//cout << "Change Locus\t" << i << endl;
						vector<double>::iterator qs = qr->second.begin();
						for(int j = 0; j<i; j++) qs++;
						//cout << endl << "before \t" << *qs << endl;
						*qs = -999;
						//cout << "AND NOW " << *qs << endl;
						
					}
					qq++;
				}
				//end if anynas
			}
			
			rs++;
			i++;
			//end first deme segment
		}
		//cout << endl;
		rr++;
		//end first deme system loop
	}
	///end first deme deme "loop"
	
	// now average across loci
	rp = preresults.begin(); //oneway
	map<string, map<string, vector<double> > > in2dememap;
	rq = rp->second.begin();
	while(rq != rp->second.end()){
		//cout << "Deme: " << rq->first << endl;
		
		//average through, but ignore the -999's !!
		map<string, vector<double> > avinSys;
		
		map<string, vector<double> >::iterator isp = rq->second.begin();
		while(isp != rq->second.end()){
			//cout << "Sytem: " << isp->first << endl;
			
			
			//cout << "wtf: " << isp->first << endl;
			unsigned int ignore = 0;
			double vectot = 0.0;
			vector<double>::const_iterator vp = isp->second.begin();
			while (vp != isp->second.end()){
				//	cout << *vp << endl;
				if(*vp == -999 || isnan(*vp)) ignore++;
				//if(*vp == -999) ignore++;
				else vectot += *vp;
				vp++;
			}
			unsigned int denom = isp->second.size() - ignore;
			double vecav;
			if (denom > 0) vecav = vectot / ((double) denom);
			else vecav = -999;
			
			vector<double> avinVec;
			avinVec.push_back(vecav);
			
			avinSys[isp->first] = avinVec;
			
			isp++;
			
		}
		in2dememap[rq->first]=avinSys;
		rq++;
	}
	
	
		x.result["ModOnewayAv"] = in2dememap;
		return x;
	
	//end fxn
}



//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs twoway tests
statresult world::twowayloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, vector<int>*, int)){
	statresult x;
	x.type = "Twoway";
	map<string, vector<int> >::const_iterator di = dmap.begin();
	
	if(dmap.size() < 2){//
		map<string, map<string, vector<double> > > in2dememap;
		map<string, vector<double> > in2deme;
		vector<double> inVec;
		inVec.push_back(-999);
		in2deme["Twoway"] = inVec;
		in2dememap["Twoway"]=in2deme;
		x.result[di->first]=in2dememap;
		
		return x;
	}

	while(di != (--dmap.end())){   //WATCH THIS BIT
		//cout << "Deme:" << di->second.name << endl;
		map<string, vector<int> >::const_iterator dj = di;
		dj++;
		map<string, map<string, vector<double> > > in2dememap;
		while(dj != dmap.end()){
			//cout << "Versus:" << dj->second.name << endl;

			map<string, vector<double> > in2deme;
			
			double vecav = 0.0;
			double avtot = 0.0;
			double vecavsz = 0.0;
			


					vector<double> inVec;

						for(unsigned int j =0; j<loci.size(); j++){


							if((loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
							

								
								vector<int> y,z;
								
								vector<int>::iterator lp, lq;
								unsigned int vecsz = 0;
								
								
								if(ancdec == "ANC"){
									vector<int> adlist = getancdecvec(&loci[j].lmap, true);
									lp = adlist.begin();
									lq = adlist.end();
									vecsz = adlist.size();
									

									
									
								}
								else if(ancdec == "DEC"){
									vector<int> adlist = getancdecvec(&loci[j].lmap, false);
									lp = adlist.begin();
									lq = adlist.end();
									vecsz = adlist.size();

									
								}
								//get a vec of all the ints that match
								else{
									
									lp = loci[j].lmap.begin();
									lq = loci[j].lmap.end();
									vecsz = loci[j].lmap.size();
									

									
								}
								
								if(vecsz >=  DEME_CUTOFF){
									
									set_intersection(di->second.begin(), di->second.end(), lp, lq, std::back_inserter(y));
									set_intersection(dj->second.begin(), dj->second.end(), lp, lq, std::back_inserter(z));
									
									//NOW RUN THE TEST ON THE llist if there are still enough to work on
									if(y.size() >=DEME_CUTOFF && z.size() >=DEME_CUTOFF) {
										inVec.push_back((*this.*tc)(&y, &z, j));
									}
									else inVec.push_back(-999); //if fails cutoff
									
									
								}
								
								else inVec.push_back(-999); 
								
								
								
								/*
								vector<int> y,z;
								
								if(ancdec == "ANC"){
									vector<int> ancxlist = getancdecvec(&loci[i].lmap, true);
									set_intersection(di->second.begin(), di->second.end(), ancxlist.begin(), ancxlist.end(), std::back_inserter(y));
									
									vector<int> ancylist = getancdecvec(&loci[j].lmap, true);
									set_intersection(di->second.begin(), di->second.end(), ancylist.begin(), ancylist.end(), std::back_inserter(z));
									
								}
								else if(ancdec == "DEC"){
									vector<int> decxlist = getancdecvec(&loci[i].lmap, true);
									set_intersection(di->second.begin(), di->second.end(), decxlist.begin(), decxlist.end(), std::back_inserter(y));
									
									vector<int> decylist = getancdecvec(&loci[j].lmap, true);
									set_intersection(di->second.begin(), di->second.end(), decylist.begin(), decylist.end(), std::back_inserter(z));
									
								}
								else{
									set_intersection(di->second.begin(), di->second.end(), loci[i].lmap.begin(), loci[i].lmap.end(), std::back_inserter(y));
									set_intersection(di->second.begin(), di->second.end(), loci[j].lmap.begin(), loci[j].lmap.end(), std::back_inserter(z));
									
								}


								//NOW RUN THE TEST ON THE llist if there are still enough to work on
								if(y.size() >=DEME_CUTOFF && z.size() >=DEME_CUTOFF) {
									inVec.push_back((*this.*tc)(&y, &z, j));
								}
								else inVec.push_back(-999); //if fails cutoff
								*/
								
								
								// end if((sl->second[k] == locustype || locustype == "ALL") && (sysi->second.lmap[k].size() >= DEME_CUTOFF) && (sysj->second.lmap[k].size() >= DEME_CUTOFF))
							}
							
							

							else{
								inVec.push_back(-999);
							}
							
						
						
					}
					
					//CHANGE HERE TO AVERAGE ACROSS ALL SYSTEMS IN DEME 1/16/06. VECTOR ENTERED WILL HAVE ONLY ONE ENTRY FROM NOW ON!
					vector<double>::const_iterator vecp = inVec.begin();
					while(vecp != inVec.end()){
						if(*vecp != -999 && !isnan(*vecp)){
							avtot += *vecp;
							vecavsz += 1.0;
						}
						vecp++;
					}
					


				//end i loop
			
			if(vecavsz > 0.0){ //i.e. at least one result was not -999 or nan
				vecav = avtot / vecavsz;
			}
			else {
				vecav = -999;
			}
			//add this average to the vector holding the averages for all the systems for this deme matching
			vector<double> avinVec;
			avinVec.push_back(vecav);
			
			in2deme["Twoway"] = avinVec;
			
			in2dememap[dj->first]=in2deme;
			dj++;
		}
		x.result[di->first]=in2dememap;
		di++;
	}
	
	
	
	return x;
}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs twoway tests
statresult world::twowayxlociloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, vector<int>*, string)){
	statresult x;
	x.type = "TwowayXloci";
	map<string, vector<int> >::iterator di = dmap.begin();
	
	if(dmap.size() < 2){//
		map<string, map<string, vector<double> > > in2dememap;
		map<string, vector<double> > in2deme;
		vector<double> inVec;
		inVec.push_back(-999);
		in2deme["TwowayXloci"] = inVec;
		in2dememap["TwowayXloci"]=in2deme;
		x.result[di->first]=in2dememap;
		
		return x;
	}

	while(di != (--dmap.end())){   //WATCH THIS BIT

		map<string, vector<int> >::iterator dj = di;
		dj++;
		map<string, map<string, vector<double> > > in2dememap;
		while(dj != dmap.end()){
			map<string, vector<double> > in2deme;
			vector<double> inVec;

			if(di->second.empty() || dj->second.empty()) inVec.push_back(-999);
			else inVec.push_back((*this.*tc)(&di->second, &dj->second, locustype));
			
			in2deme["TwowayXloci"] = inVec;
			
			in2dememap[dj->first]=in2deme;
			dj++;
			//end while dj
		}
		x.result[di->first]=in2dememap;
		di++;
	}
	
	return x;
}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs tests on two-locus comparisons
statresult world::twolociloop (string locustype, string ancdec, double (world::*tc) (vector<int>* , int , int )){
	statresult x;
	x.type = "TwoLoci";
	map<string, vector<int> >::iterator di = dmap.begin();
	

	map<string, map<string, vector<double> > > indeme;
	while(di != dmap.end()){


				map<string, vector<double> > inmsystem;
		
		for(unsigned int i =0; i<loci.size(); i++){
			
			vector<double> oneloci;
			
			for(unsigned int j=(i+1); j<loci.size(); j++){
		


						if((loci[i].type == locustype || locustype == "ALL") && (loci[i].lmap.size() >= DEME_CUTOFF) && (loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 

							vector<int> y,z;
							
							vector<int>::iterator ip, iq, jp, jq;
							unsigned int vecszi, vecszj = 0;
							
							
							if(ancdec == "ANC"){
								vector<int> iadlist = getancdecvec(&loci[i].lmap, true);
								ip = iadlist.begin();
								iq = iadlist.end();
								vecszi = iadlist.size();
								
								vector<int> jadlist = getancdecvec(&loci[j].lmap, true);
								jp = jadlist.begin();
								jq = jadlist.end();
								vecszj = jadlist.size();
								
								
							}
							else if(ancdec == "DEC"){
								vector<int> iadlist = getancdecvec(&loci[i].lmap, false);
								ip = iadlist.begin();
								iq = iadlist.end();
								vecszi = iadlist.size();
								
								vector<int> jadlist = getancdecvec(&loci[j].lmap, false);
								jp = jadlist.begin();
								jq = jadlist.end();
								vecszj = jadlist.size();
								
								
							}
							//get a vec of all the ints that match
							else{
								
								ip = loci[i].lmap.begin();
								iq = loci[i].lmap.end();
								vecszi = loci[i].lmap.size();
								
								jp = loci[j].lmap.begin();
								jq = loci[j].lmap.end();
								vecszj = loci[j].lmap.size();
								
								
								
							}
							
							if(vecszi >=  DEME_CUTOFF && vecszj >=  DEME_CUTOFF){
								
								set_intersection(di->second.begin(), di->second.end(), ip, iq, std::back_inserter(y));
								set_intersection(di->second.begin(), di->second.end(), jp, jq, std::back_inserter(z));
								
								//NOW RUN THE TEST ON THE llist if there are still enough to work on
								if(y.size() >=DEME_CUTOFF && z.size() >=DEME_CUTOFF) {
								//	cout << ((*this.*tc)(&y,i,j)) << endl;
									oneloci.push_back((*this.*tc)(&y,i,j));
								}
								else oneloci.push_back(-999); //if fails cutoff
								
								
							}
							
							else oneloci.push_back(-999);
							

							
							
						}
						else{
							oneloci.push_back(-999);
						}
					}
					stringstream ss;
					ss << i;
					
			/*
					cout << "Put in: " << ss.str() << ": ";
			vector<double>::const_iterator wth = oneloci.begin();
			while(wth != oneloci.end()){
					cout << *wth << " ";
				
				wth++;
			}
			cout << endl;
			*/
			
					inmsystem[ss.str()] = oneloci; //loccus number becomes name in the map
				}
				if(!inmsystem.empty()){//only enter systems taht had more than one locus
					indeme["TwoLoci"]=inmsystem;
				}

			

		
		
		x.result[di->first]=indeme;
		di++;
	}
	return x;
}


//---------------------------------------------
//Function Name: foreach
//Function Definition: loop to iterate through demes and performs tests on two-locus comparisons
statresult world::twolociadjloop (string locustype, string ancdec, double (world::*tc) (vector<int>* , int , int )){
	statresult x;
	x.type = "TwoLociAdj";
	map<string, vector<int> >::iterator di = dmap.begin();
	
	
	map<string, map<string, vector<double> > > indeme;
	while(di != dmap.end()){
		
		
		map<string, vector<double> > inmsystem;
		
		for(unsigned int i =0; i<(loci.size()-1); i++){
			
			vector<double> oneloci;
			
			unsigned int j=(i+1);
				
				
				
				if((loci[i].type == locustype || locustype == "ALL") && (loci[i].lmap.size() >= DEME_CUTOFF) && (loci[j].type == locustype || locustype == "ALL") && (loci[j].lmap.size() >= DEME_CUTOFF)){ 
					
					vector<int> y,z;
					
					vector<int>::iterator ip, iq, jp, jq;
					unsigned int vecszi, vecszj = 0;
					
					
					if(ancdec == "ANC"){
						vector<int> iadlist = getancdecvec(&loci[i].lmap, true);
						ip = iadlist.begin();
						iq = iadlist.end();
						vecszi = iadlist.size();
						
						vector<int> jadlist = getancdecvec(&loci[j].lmap, true);
						jp = jadlist.begin();
						jq = jadlist.end();
						vecszj = jadlist.size();
						
						
					}
					else if(ancdec == "DEC"){
						vector<int> iadlist = getancdecvec(&loci[i].lmap, false);
						ip = iadlist.begin();
						iq = iadlist.end();
						vecszi = iadlist.size();
						
						vector<int> jadlist = getancdecvec(&loci[j].lmap, false);
						jp = jadlist.begin();
						jq = jadlist.end();
						vecszj = jadlist.size();
						
						
					}
					//get a vec of all the ints that match
					else{
						
						ip = loci[i].lmap.begin();
						iq = loci[i].lmap.end();
						vecszi = loci[i].lmap.size();
						
						jp = loci[j].lmap.begin();
						jq = loci[j].lmap.end();
						vecszj = loci[j].lmap.size();
						
						
						
					}
					
					if(vecszi >=  DEME_CUTOFF && vecszj >=  DEME_CUTOFF){
						
						set_intersection(di->second.begin(), di->second.end(), ip, iq, std::back_inserter(y));
						set_intersection(di->second.begin(), di->second.end(), jp, jq, std::back_inserter(z));
						
						//NOW RUN THE TEST ON THE llist if there are still enough to work on
						if(y.size() >=DEME_CUTOFF && z.size() >=DEME_CUTOFF) {
							//	cout << ((*this.*tc)(&y,i,j)) << endl;
							oneloci.push_back((*this.*tc)(&y,i,j));
						}
						else oneloci.push_back(-999); //if fails cutoff
						
						
					}
					
					else oneloci.push_back(-999);
					
					
					
					
				}
				else{
					oneloci.push_back(-999);
				}
			
			stringstream ss;
			ss << i;
			

			
			inmsystem[ss.str()] = oneloci; //loccus number becomes name in the map
		}
		if(!inmsystem.empty()){//only enter systems taht had more than one locus
			indeme["TwoLociAdj"]=inmsystem;
		}
		
		
		
		
		
		x.result[di->first]=indeme;
		di++;
	}
	return x;
}



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------
//
// PRINTING LOOPS
//
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------



//---------------------------------------------
//Function Name: printresult
//Function Definition: output for oneway
void world::persystemprint(statresult* sr, ostream &out){

	map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
	s = sr->result.begin();
    out << "PerSystem" << endl;

	out << "Deme#\tDemeName\t";
    p = s->second.begin();  
	
	q = p->second.begin();
	while(q != p->second.end()){
		
		
		out << q->first << "\t";
		
		
		q++;
	}
    out << endl;
    
	
    
    p = s->second.begin(); 
	int demenum = 0;
    while(p != s->second.end()){
        out << demenum << "\t" << p->first << "\t";
        q = p->second.begin();
        while(q != p->second.end()){
            r = q->second.begin();
			
			if(*r == -999) out << "NA" << "\t";
			else out << *r << "\t";
            
            q++;
        }
        p++;
		demenum++;
        out << endl;
        
    }
    out << endl;
}

//---------------------------------------------
//Function Name: printresult
//Function Definition: output for oneway
void world::bydemeprint(statresult* sr, ostream &out){
	
	map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
	s = sr->result.begin();
    out << "ByDeme" << endl;
	
	out << "Deme#\tDemeName\t";
    p = s->second.begin();  
	
	q = p->second.begin();
	while(q != p->second.end()){
		
		
		out << q->first << "\t";
		
		
		q++;
	}
    out << endl;
    
	
    
    p = s->second.begin(); 
	int demenum = 0;
    while(p != s->second.end()){
        out << demenum << "\t" << p->first << "\t";
        q = p->second.begin();
        while(q != p->second.end()){
            r = q->second.begin();
			
			if(*r == -999) out << "NA" << "\t";
			else out << *r << "\t";
            
            q++;
        }
        p++;
		demenum++;
        out << endl;
        
    }
    out << endl;
}

//---------------------------------------------
//Function Name: printresult
//Function Definition: output for 
void world::perlocusxpopprint(statresult* sr, ostream &out){
	map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
	s = sr->result.begin();
    out << "PerLocusXPop" << endl; //per
	out << "System\t->\t";
    p = s->second.begin();  
	
	q = p->second.begin();
	while(q != p->second.end()){
		r = q->second.begin();
		int locusnum = 0;
		while(r != q->second.end()){
			
			out << q->first << "\t";
			r++;
			locusnum++;
		}
		
		q++;
	}
    out << endl;
    
	
    //first go through for headers
    out << "\t\t";
    p = s->second.begin();  
	q = p->second.begin();
	while(q != p->second.end()){
		r = q->second.begin();
		int locusnum = 0;
		while(r != q->second.end()){
			
			out << locusnum << "-" <<loci[locusnum].type << "\t";
			r++;
			locusnum++;
		}
		
		q++;
	}
    out << endl;
	
    
    p = s->second.begin();  
    while(p != s->second.end()){
		
		
        out  << "\t"  << "\t";
		
        q = p->second.begin();
        while(q != p->second.end()){
            r = q->second.begin();
            int locusnum = 0;
            while(r != q->second.end()){
				if(*r == -999) out << "NA" << "\t";
				else out << *r << "\t";
                r++;
                locusnum++;
            }
            
            q++;
        }
        p++;
        out << endl;
        
    }
    out << endl;
	
}


//---------------------------------------------
//Function Name: printresult
//Function Definition:
void world::perlocusxpopavprint(statresult* sr, ostream &out){
	map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
	s = sr->result.begin();
    out << "PerLocusXPopAv" << endl; //per
	p = s->second.begin(); 
	q = p->second.begin();
	r = q->second.begin();
	out << "Average\t->\t" << *r << endl << endl;
	
}

//---------------------------------------------
//Function Name: printresult
//Function Definition: output for oneway
void world::onewayprint(statresult* sr, ostream &out)
//void oneway::printresult(ostream &out, string locustype)
{
	map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
	s = sr->result.begin();
    out << "Oneway" << endl;
	out << "System\t->\t";
    p = s->second.begin();  
	
	q = p->second.begin();
	while(q != p->second.end()){
		r = q->second.begin();
		int locusnum = 0;
		while(r != q->second.end()){
			
			out << q->first << "\t";
			r++;
			locusnum++;
		}
		
		q++;
	}
    out << endl;
    
	
    //first go through for headers
    out << "Deme#\tDemeName\t";
    p = s->second.begin();  
	q = p->second.begin();
	while(q != p->second.end()){
		r = q->second.begin();
		int locusnum = 0;
		while(r != q->second.end()){
			
			out << locusnum << "-" << loci[locusnum].type << "\t";
			
			r++;
			locusnum++;
		}
		
		q++;
	}
    out << endl;
	
    
    p = s->second.begin();  
	int demenum = 0;
    while(p != s->second.end()){

		
		out << demenum << "\t" << p->first << "\t";
		
        q = p->second.begin();
        while(q != p->second.end()){
            r = q->second.begin();
            int locusnum = 0;
            while(r != q->second.end()){
				if(*r == -999) out << "NA" << "\t";
				else out << *r << "\t";
                r++;
                locusnum++;
            }
            
            q++;
        }
        p++;
		demenum++;
        out << endl;
        
    }
    out << endl;
}


//---------------------------------------------
//Function Name: printresult
//Function Definition: output for oneway
void world::onewayavprint(statresult* sr, ostream &out)

{
	map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
	s = sr->result.begin();
    out << "OnewayAv" << endl;
	out << "Type\t->\t";
    p = s->second.begin();  
	
	q = p->second.begin();
	while(q != p->second.end()){
		r = q->second.begin();
		int locusnum = 0;
		while(r != q->second.end()){
			
			out << q->first << "\t";
			r++;
			locusnum++;
		}
		
		q++;
	}
    out << endl;
    
	
    //first go through for headers
    out << "Deme#\tDemeName\t";
    out << endl;
	
    
    p = s->second.begin();  
	int demenum = 0;
    while(p != s->second.end()){
		
		out << demenum << "\t" << p->first << "\t";
		
        q = p->second.begin();
        while(q != p->second.end()){
            r = q->second.begin();
            int locusnum = 0;
            while(r != q->second.end()){
				if(*r == -999) out << "NA" << "\t";
				else out << *r << "\t";
                r++;
                locusnum++;
            }
            
            q++;
        }
        p++;
		demenum++;
        out << endl;
        
    }
    out << endl;
}

	
	//---------------------------------------------
	//Function Name: printresult
	//Function Definition: output for oneway
void world::modonewayavprint(statresult* sr, ostream &out)
	
	{
		map<string, map<string, map<string, vector<double> > > >::const_iterator s;
		map<string, map<string, vector<double> > >::const_iterator p;
		map<string, vector<double> >::const_iterator q;
		vector<double>::const_iterator r;
		
		s = sr->result.begin();
		out << "OnewayAv" << endl;
		out << "Type\t->\t";
		p = s->second.begin();  
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			int locusnum = 0;
			while(r != q->second.end()){
				
				out << q->first << "\t";
				r++;
				locusnum++;
			}
			
			q++;
		}
		out << endl;
		
		
		//first go through for headers
		out << "Deme#\tDemeName\t";
		out << endl;
		
		
		p = s->second.begin();  
		int demenum = 0;
		while(p != s->second.end()){
			
			out << demenum << "\t" << p->first << "\t";
			
			q = p->second.begin();
			while(q != p->second.end()){
				r = q->second.begin();
				int locusnum = 0;
				while(r != q->second.end()){
					if(*r == -999) out << "NA" << "\t";
					else out << *r << "\t";
					r++;
					locusnum++;
				}
				
				q++;
			}
			p++;
			demenum++;
			out << endl;
			
		}
		out << endl;
	}
	

//---------------------------------------------
//Function Name: printresult
//Function Definition: out putfor twoway
void world::twowayprint(statresult* sr, ostream &out)	
{
    map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
    
    //----------
	s = sr->result.begin();
    out << "Twoway" << endl;
	
	out << "Deme#\t->\t";
	
	p = s->second.begin();
	int versusnum = 1;
	while(p != s->second.end()){
		//versus

		
		
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			int locusnum = 0;
			while(r != q->second.end()){

				out << versusnum << "\t";
				r++;
				locusnum++;
				
			}
			q++;
		}
		p++;
		versusnum++;
		
	}
	
	out << endl;
	
	out << "->\tDemeName\t";
	
	p = s->second.begin();
	
	while(p != s->second.end()){
		//versus
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			
			while(r != q->second.end()){
				out << p->first << "\t";
				r++;
				
			}
			q++;
		}
		p++;
		
	}
	
    out << endl;
	
	
	
	
	//spacer needed.. but multiplied by total number of systems/loci
	int spaces = 0;
    
	int demenum = 0;
    //----------
	s = sr->result.begin();
    while(s != sr->result.end()){
		p = s->second.begin();


        out << demenum << "\t" << s->first << "\t";
		
		for(int i=0; i<spaces; i++) out << ".\t";  
		
        while(p != s->second.end()){
			//versus   
            q = p->second.begin();
            while(q != p->second.end()){
                
                r = q->second.begin();
                
                while(r != q->second.end()){
					
					
					if(*r == -999) out << "NA" << "\t";
					else out << *r << "\t";
					
                    r++;
					
                    
                }
                q++;
            }
            p++;
        }
        out << endl; //to make a two-way table
		spaces++;
		demenum++;
        s++;
    }
    out << endl;
}


//---------------------------------------------
//Function Name: printresult
//Function Definition: out putfor twoway
void world::twowayxlociprint(statresult* sr, ostream &out)	
{
    map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
    
    //----------
	s = sr->result.begin();
    out << "TwowayXloci" << endl;
	
	out << "Deme#\t->\t";
	
	p = s->second.begin();
	int versusnum = 1;
	while(p != s->second.end()){
		//versus
		
		
		
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			int locusnum = 0;
			while(r != q->second.end()){
				
				out << versusnum << "\t";
				r++;
				locusnum++;
				
			}
			q++;
		}
		p++;
		versusnum++;
		
	}
	
	out << endl;
	
	out << "->\tDemeName\t";
	
	p = s->second.begin();
	
	while(p != s->second.end()){
		//versus
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			
			while(r != q->second.end()){
				out << p->first << "\t";
				r++;
				
			}
			q++;
		}
		p++;
		
	}
	
    out << endl;
	
	
	
	
	//spacer needed.. but multiplied by total number of systems/loci
	int spaces = 0;
    
	int demenum = 0;
    //----------
	s = sr->result.begin();
    while(s != sr->result.end()){
		p = s->second.begin();
		
		
        out << demenum << "\t" << s->first << "\t";
		
		for(int i=0; i<spaces; i++) out << ".\t";  
		
        while(p != s->second.end()){
			//versus   
            q = p->second.begin();
            while(q != p->second.end()){
                
                r = q->second.begin();
                
                while(r != q->second.end()){
					
					
					if(*r == -999) out << "NA" << "\t";
					else out << *r << "\t";
					
                    r++;
					
                    
                }
                q++;
            }
            p++;
        }
        out << endl; //to make a two-way table
		spaces++;
		demenum++;
        s++;
    }
    out << endl;
}

//---------------------------------------------
//Function Name: printresult
//Function Definition: out putfor twoway
void world::bytwodemeprint(statresult* sr, ostream &out)	
{
    map<string, map<string, map<string, vector<double> > > >::const_iterator s;
    map<string, map<string, vector<double> > >::const_iterator p;
    map<string, vector<double> >::const_iterator q;
    vector<double>::const_iterator r;
	
    
    //----------
	s = sr->result.begin();
    out << "ByTwoDeme" << endl;
	
	out << "Deme#\t->\t";
	
	p = s->second.begin();
	int versusnum = 1;
	while(p != s->second.end()){
		//versus
		
		
		
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			int locusnum = 0;
			while(r != q->second.end()){
				
				out << versusnum << "\t";
				r++;
				locusnum++;
				
			}
			q++;
		}
		p++;
		versusnum++;
		
	}
	
	out << endl;
	
	out << "->\tDemeName\t";
	
	p = s->second.begin();
	
	while(p != s->second.end()){
		//versus
		
		q = p->second.begin();
		while(q != p->second.end()){
			r = q->second.begin();
			
			while(r != q->second.end()){
				out << p->first << "\t";
				r++;
				
			}
			q++;
		}
		p++;
		
	}
	
    out << endl;
	
	
	
	
	//spacer needed.. but multiplied by total number of systems/loci
	int spaces = 0;
    
	int demenum = 0;
    //----------
	s = sr->result.begin();
    while(s != sr->result.end()){
		p = s->second.begin();
		
		
        out << demenum << "\t" << s->first << "\t";
		
		for(int i=0; i<spaces; i++) out << ".\t";  
		
        while(p != s->second.end()){
			//versus   
            q = p->second.begin();
            while(q != p->second.end()){
                
                r = q->second.begin();
                
                while(r != q->second.end()){
					
					
					if(*r == -999) out << "NA" << "\t";
					else out << *r << "\t";
					
                    r++;
					
                    
                }
                q++;
            }
            p++;
        }
        out << endl; //to make a two-way table
		spaces++;
		demenum++;
        s++;
    }
    out << endl;
}


//---------------------------------------------
//Function Name: printresult
//Function Definition: out putfor twoway
void world::twolociprint(statresult* sr, ostream &out){
	unsigned int i;
	map<string, map<string, map<string, vector<double> > > >::const_iterator s; //demes
	map<string, map<string, vector<double> > >::const_iterator p; //msystems
	map<string, vector<double> >::const_iterator q; //oneloci at a time
	vector<double>::const_iterator r; // compared loci
	
	
	
	out << "TwoLoci\n";
	
	
	//deme numbers
	out << "Deme#\t";
	s = sr->result.begin();
	int demenum = 0;
	while(s != sr->result.end()){
		
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << demenum << "\t";
					
					r++;
				}
				q++;
			}
			p++;
			
		}
		s++;
		demenum++;
	}
	out << endl;
	
	
	//deme names
	out << "DemeName\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << s->first << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
	
	
	// system names
	out << "System\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << p->first << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
    
	
	
	// loci titles
	out << "Locus Pair\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << q->first << "<->" << i << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
	
	
	
	// results
	out << "Results\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					if(*r == -999) out << "NA" << "\t";
					else out << *r << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
	
}



//---------------------------------------------
//Function Name: printresult
//Function Definition: out putfor twoway
void world::twolociadjprint(statresult* sr, ostream &out){
	unsigned int i;
	map<string, map<string, map<string, vector<double> > > >::const_iterator s; //demes
	map<string, map<string, vector<double> > >::const_iterator p; //msystems
	map<string, vector<double> >::const_iterator q; //oneloci at a time
	vector<double>::const_iterator r; // compared loci
	
	
	
	out << "TwoLociAdj\n";
	
	
	//deme numbers
	out << "Deme#\t";
	s = sr->result.begin();
	int demenum = 0;
	while(s != sr->result.end()){
		
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << demenum << "\t";
					
					r++;
				}
				q++;
			}
			p++;
			
		}
		s++;
		demenum++;
	}
	out << endl;
	
	
	//deme names
	out << "DemeName\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << s->first << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
	
	
	// system names
	out << "System\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << p->first << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
    
	
	
	// loci titles
	out << "Locus Pair\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					out << q->first << "<->" << i << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
	
	
	
	// results
	out << "Results\t";
	s = sr->result.begin();
	while(s != sr->result.end()){
		p = s->second.begin();
		while(p != s->second.end()){
			q = p->second.begin();
			while(q != p->second.end()){
				i = atoi(q->first.c_str());
				r = q->second.begin();
				while(r != q->second.end()){
					i++;
					if(*r == -999) out << "NA" << "\t";
					else out << *r << "\t";
					r++;
				}
				q++;
			}
			p++;
		}
		s++;
	}
	out << endl;
	
}



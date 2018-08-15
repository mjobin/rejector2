//------------------------------------------------
//
// REJECTOR
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropological Sciences
// Stanford University, 2006
//
//------------------------------------------------

#include <sstream>
#include <cmath>
#include "rejector.h"

using namespace std;


//---------------------------------------------
//Function Name:
//Function Definition: ensures all parts of a string are integers
bool rtestr::stralphaTest(const char* Str) const
{
    const int len = strlen(Str);
    for(int i = 0;i < len;i++) if (!isalpha(Str[i])) return false;
    return true;
}


//---------------------------------------------
//Function Name:
//Function Definition: tests if any of a string is a digit
bool rtestr::anydigits(const char* Str) const
{
    const int len = strlen(Str);
    for(int i = 0;i < len;i++) if (isdigit(Str[i])) return true;
    return false;
}

//---------------------------------------------
//Function Name:
//Function Definition:
bool istresultempty(const tresult *r){
	tresult::const_iterator rp = r->begin();
	
	
	while(rp != r->end()){
		map<string, map<string, vector<double> > >::const_iterator rq = rp->second.begin();
		
		
		while(rq != rp->second.end()){
			map<string, vector<double> >::const_iterator rr = rq->second.begin();
			
			
			while(rr != rq->second.end()){
				if(!rr->second.empty()) return false; // if any vecotr is not empty return false
				
				rr++;
			}
			rq++;
		}
		rp++;
	}
	
	return true;
	
}

//---------------------------------------------
//Function Name:
//Function Definition: takes two worlds (usually r is real data nd s is a sim) and compares their  results
rejresult rejectcompare(world *r, world *s, priordist* pt){
	rejresult o;
	list<string> acc;
	double alpha;
	bool anyaccepted = false;
	bool anyrejected = false;
	
	map<string, statresult>::const_iterator rtests = r->stats.begin();
	while(rtests != r->stats.end()){
		map<string, statresult>::const_iterator stests = s->stats.find(rtests->first);
		//cout << "Testing " << rtests->first << endl;
		
		if(stests != s->stats.end()){ //make sure its not returning an empty test.. there shoud be a match for each test
									//	cout << " against " << stests->first << endl;
			
				if (rtests->first.find("SampleSize") != string::npos){ //just put X's for smaple size tests
				o.resultacc[rtests->first] = "X";
				o.numacc[rtests->first] = -999;
				o.numtot[rtests->first] = 1;
				}
			
			else{
				
				if(pt->alphalist.find(rtests->first) != pt->alphalist.end()) alpha = pt->alphalist[rtests->first];
				else if (pt->alphalist.find("ALL") != pt->alphalist.end()) alpha = pt->alphalist["ALL"];
				else{
					cerr << "Error in rejector! Cannot find test name " << rtests->first << " alpha value nor a value for all tests!" << endl;
					abort();
				}
				
				
				
				

					if(rejectcomparestats(r, s, alpha, rtests->first, &o)){
						o.resultacc[rtests->first] = "T";
						
						//cout << "*****ACCEPTED****  " << "NUMBERACC: " << o.numacc << " NUMBERTOT: " << o.numtot << endl << endl;

						
						acc.push_back(rtests->first);
						anyaccepted = true;
					}
					else{//case where comparison run but all data is missing
						if (o.numtot[rtests->first] == 0) {
							o.resultacc[rtests->first] = "X";
						}
						else{
						o.resultacc[rtests->first] = "F";
						anyrejected = true;
						}
					}
				
				
				
				

				
	
				
			}
			//end if not empty test
		}
		
		else{
			cout << endl << "No match in simulated test list for " << rtests->first << endl;
			o.resultacc[rtests->first] = "X";
			o.numacc[rtests->first] = -999;
			o.numtot[rtests->first] = 1;
		}
		rtests++;
	}
	
	return o;
	

}

//---------------------------------------------
//Function Name:
//Function Definition:
bool rejresult::rejectordecision(){
if(allornothing){
	if(anyaccepted && !anyrejected) return true;
	else return false;
}
else return anyaccepted;
}

//---------------------------------------------
//Function Name:
//Function Definition: Called by compare worlds to compare the results of two tests results of the same type
bool rejectcomparestats(world* r, world* s, double alpha, string testname, rejresult* o){
	map<string, statresult>::const_iterator rtests = r->stats.find(testname);
	map<string, statresult>::const_iterator stests = s->stats.find(testname);
	
	if(rtests == r->stats.end() || stests == s->stats.end()){
		cerr << "Error in rejector: compareresults. Comparertestrs found " << testname << " in real and simulated results, but compareresults did not! This indicates a programming error.";
		abort();
	}
	
	//cout << "rejectcomparestats on " << testname << endl;
	
	sresult::const_iterator rp = rtests->second.result.begin();

	int numberacc = 0;
	int numbertot = 0;
	bool alltrue = true;
	
	
	bool calcguard = false;  //turned true once there has been at least one accepted calculation, so that you can get a false positive with all skippecd comparisons
	
	while(rp != rtests->second.result.end()){
		sresult::const_iterator sp = stests->second.result.find(rp->first);
		if(sp != stests->second.result.end()){
			//cout << "rp: " << rp->first << "  sp: " << sp->first << endl;
			map<string, map<string, vector<double> > >::const_iterator rq = rp->second.begin();
			
			while(rq != rp->second.end()){
				map<string, map<string, vector<double> > >::const_iterator sq = sp->second.find(rq->first);
				//cout << "rq: " << rq->first << endl;
				if(sq != sp->second.end()){
					//cout << "rq: " << rq->first << "  sq: " << sq->first << endl;
					map<string, vector<double> >::const_iterator rr = rq->second.begin();
					
					while(rr != rq->second.end()){
						map<string, vector<double> >::const_iterator sr = sq->second.find(rr->first);
						if(sr != sq->second.end()){
							//cout << "rr: " << rr->first << "  sr: " << sr->first << endl;
							vector<double>::const_iterator rs = rr->second.begin(); 
							vector<double>::const_iterator ss = sr->second.begin();//no if statemtn here as its assumed vecotrs are of same length
								while(rs != rr->second.end()){
									//cout << "rs: " << *rs << "  ss: " << *ss << endl;
									//if(*rs != -999 && *ss != -999 && *rs != 0.0){//don't test if either datum is missing or if there's a division by zero
									if(*rs != -999 && *ss != -999){//don't test if either datum is missing 
										numbertot++;
										double delta;
										calcguard = true;
										
										delta = abs((*rs-*ss)/(*rs)); 
										
										if (delta > alpha) alltrue = false;
										
										else numberacc++;
										
										
										
									}
									
									
									rs++;
									ss++;
								}
								
								
								//if sr	
						}
						
						rr++;
						//cout << endl;
					}
					
					//if sq	
				}
				
				rq++;
			}
			
			//if sp	
		}
		
		rp++;
	}
	

	
	o->numacc[testname] = numberacc;
	o->numtot[testname] = numbertot;


	
	if(calcguard && alltrue) return true;
	else return false;
}


//---------------------------------------------
//Function Name:
//Function Definition: takes mean at each result
statmap averagestatresults(list<statmap> avlist){



	
	statmap avout;
	statmap testnums;
	

	int i =0;
	list<statmap>::const_iterator avp = avlist.begin();
	bool firsttime = true;
	while(avp != avlist.end()){
		//cout << "**********" << endl << "num: " << i << endl;
		statmap::const_iterator avt = avp->begin();
		while(avt != avp->end()){
			//cout << endl << "test name: " << avt->first << " Type: " << avt->second.type << endl;

			if(firsttime){
				//create empty result
				//vector<double> xs;
				//map<string, vector<double> > xr;
				//map<string, map<string, vector<double> > > xq;
				//sresult xp;
				
				
				
				statresult newstatresult;
				newstatresult.type = avt->second.type;
				sresult newsresult;
				newstatresult.result = newsresult;
				avout[avt->first] = newstatresult;
				

				
				statresult newtestnum;
				newtestnum.type = avt->second.type;
				sresult othernewsresult;
				newtestnum.result = othernewsresult;
				testnums[avt->first] = newtestnum;
				
			}
			
			
			//avtresult thisavtest = fived[avt->first];
			sresult thissresult = avout[avt->first].result;
			sresult thissnum = testnums[avt->first].result;
			
			sresult::const_iterator tp = avt->second.result.begin();
			//avtresult::iterator tvp = thisavtest.begin();
			
			while(tp != avt->second.result.end()){
				

				
				map<string, map<string, vector<double> > >::const_iterator tq = tp->second.begin();
				//map<string, map<string, vector<vector<double> > > >::iterator tvq = tvp->second.begin();
				//map<string, map<string, vector<vector<double> > > > thisq = thisavtest[tp->first];
				map<string, map<string, vector<double> > > thisq = thissresult[tp->first];
				map<string, map<string, vector<double> > > tnumq = thissnum[tp->first];
				
				while(tq != tp->second.end()){
					

					
					map<string, vector<double> >::const_iterator tr = tq->second.begin();
					//map<string, vector<vector<double> > >::iterator tvr = tvq->second.begin();
					//map<string, vector<vector<double> > > thisr = thisq[tq->first];
					map<string, vector<double> > thisr = thisq[tq->first];
					map<string, vector<double> > tnumr = tnumq[tq->first];
					
					while(tr != tq->second.end()){
						

						
						vector<double>::const_iterator ts = tr->second.begin();
						//vector<vector<double> >::iterator tvs = tvr->second.begin();
						//vector<vector<double> > thiss = thisr[tr->first];
						vector<double> thiss = thisr[tr->first];
						vector<double> tnums = tnumr[tr->first];
						
						vector<double>::iterator tvs = thiss.begin();
						vector<double>::iterator tns = tnums.begin();
						

						
						while(ts != tr->second.end()){
							
							//cout << *ts << endl;
	
							if(firsttime){
								
								if(*ts == -999 || isnan(*ts)){
									thiss.push_back(-999);
									tnums.push_back(0.0);
								}
								else{


									thiss.push_back(*ts);
									tnums.push_back(1.0);
								}
			


									
								ts++;
							}
							else{
								
								if(*ts != -999 && !isnan(*ts)){
									//now increase 

									*tns = *tns + 1.0;

									if(*tvs != -999){
										double thissfrac;
										if(*tvs == 0.0) thissfrac = 0.0;
										else thissfrac = *tvs * ((*tns-1) / *tns);
										double tsfrac;
										if(*ts == 0.0) tsfrac = 0.0;
										else tsfrac = *ts * (1/ *tns);
										double aved = thissfrac + tsfrac;
										*tvs = aved;
									}
									else {
										*tvs = *ts; //if the current average is -99 all NA's so far so just replace with that value
										
									}
									
								}
								ts++;
								
								tvs++;
								tns++;
								
							}
								
						

							
						

							
							
						}

						thisr[tr->first] = thiss;
						tnumr[tr->first] = tnums;
						tr++;
						
					}
					thisq[tq->first] = thisr;
					tnumq[tq->first] = tnumr;
					tq++;
					
				}
				thissresult[tp->first] = thisq;
				thissnum[tp->first] = tnumq;
				tp++;
				
			}
			
			
			avout[avt->first].result = thissresult;
			testnums[avt->first].result = thissnum;
			avt++;
		}
		firsttime = false;
		avp++;
		i++;
	}
	
	
	
	
	
	return avout;
}




//---------------------------------------------
//Function Name:
//Function Definition: spews out the averages. So very quick, and so very dirty
void avdistout(list<statmap> avlist, ostream &dout){
	
	list<statmap>::const_iterator avp = avlist.begin();
	bool firsttime = true;
	
	
	
	
	
	
	while(avp != avlist.end()){
		statmap::const_iterator avt = avp->begin();

		
		if(firsttime == true){
			

			
			while(avt != avp->end()){
				

				
				tresult::const_iterator tp = avt->second.result.begin();
				
				while(tp != avt->second.result.end()){
					
					map<string, map<string, vector<double> > >::const_iterator tq = tp->second.begin();
					int outside = 1;
					while(tq != tp->second.end()){
						
						map<string, vector<double> >::const_iterator tr = tq->second.begin();
						int inside = outside+1;
						while(tr != tq->second.end()){
							
							vector<double>::const_iterator ts = tr->second.begin();
							
							while(ts != tr->second.end()){
								
								if (avt->second.type == "Twoway"){
								
									dout << avt->first << "-" << outside << "-" << inside << "\t";
									
								}
								else {
									
									dout << avt->first << "-"<< outside  << "\t";
									
								}
									
								
								ts++;
							}
							inside++;
							tr++;
						}
						outside++;
						tq++;
					}
					
					tp++;
				}
				avt++;
			}
			
			dout << endl;
			
			//end if firsttime
			firsttime = false;
			avt = avp->begin(); //reset so first results actually print!
		}
		
		//lets do the numbers
		
		while(avt != avp->end()){
			
			sresult::const_iterator tp = avt->second.result.begin();
	
				while(tp != avt->second.result.end()){
					
					map<string, map<string, vector<double> > >::const_iterator tq = tp->second.begin();
					
					while(tq != tp->second.end()){
						
						map<string, vector<double> >::const_iterator tr = tq->second.begin();
						
						while(tr != tq->second.end()){
							
							vector<double>::const_iterator ts = tr->second.begin();
							
							while(ts != tr->second.end()){
								
								dout << *ts << "\t";
								
								ts++;
							}
							tr++;
						}
						tq++;
					}
					tp++;
				}
				avt++;
		}
		dout << endl;
		avp++;
	}
	//end fxn
}



//---------------------------------------------
//Function Name:
//Function Definition:
int compareworlds(world* r, world* s, priordist* pt, ostream &rout, bool hdr, bool keepstats){

	if(comparetype == 0){
		rejresult o = rejectcompare(r, s, pt);
		if(hdr) s->outfileprint(rout, true, pt, &o, keepstats); //print header the first time
		
		if(o.rejectordecision() || alwaysprintoutputline) {  // Print an output line if the test accepts OR if instructed to print a line for each iteration
			successcount++;
			s->outfileprint(rout, false, pt, &o, keepstats);
			
		}
	
		
		return 0;
	}
	
	//Another system may allow for the alteration of the priordist here
	
	
	cerr << "Unknown comparison type " << comparetype << endl;
	return -1;
}

//---------------------------------------------
//Function Name: outfileprint
//Function Definition: prints outfile
void world::outfileprint(ostream &rout, bool hdr, priordist* pt, rejresult* o, bool keepstats){
	//read in from .par for spacing
	string buf;
	int i, j;
	int mm;
	
	//map<string, vector<int> >::const_iterator dp = dmap.begin();
	

	
	map<string, double>::const_iterator dsp = tparams.demesize.begin();
	

	
	while(dsp != tparams.demesize.end()){
		if(hdr) rout << "DemeSize-" << dsp->first << "\t";
		else rout << dsp->second << "\t";
		dsp++;
	}
	
	

	

	map<string, int>::const_iterator ssp = tparams.samplesize.begin();
	while(ssp != tparams.samplesize.end()){
		if(hdr) rout << "SampSize-" << ssp->first << "\t";
		else rout << ssp->second << "\t";
		ssp++;
	}
	
	

	

	map<string, double>::const_iterator gsp = tparams.growthrate.begin();
	while(gsp != tparams.growthrate.end()){
		if(hdr) rout << "Gr.Rate-" << gsp->first << "\t";
		else rout << gsp->second << "\t";
		gsp++;

	}
	
	list<vector<vector <double> > >::const_iterator mmp = tparams.migmat.begin();
	mm = 0;
	while(mmp != tparams.migmat.end()){
		i = 0;
		vector<vector<double> >::const_iterator mmpi = mmp->begin();
		while(mmpi != mmp->end()){
			vector<double>::const_iterator mmpj = mmpi->begin();
			j=0;
			while(mmpj != mmpi->end()){
				if(hdr) rout << "MM" << mm << "-" << i << "," << j << "\t";
				else rout << *mmpj << "\t";
				mmpj++;
				j++;
			}
			mmpi++;
			i++;
		}
		
		
		mmp++;
		mm++;
	}
	
	
	
	
	
	i = 0;
	list<vector<double> >::const_iterator hsp = tparams.histev.begin();
	while(hsp != tparams.histev.end()){
		vector<double>::const_iterator hhsp = hsp->begin();
		j = 0;
		while(hhsp != hsp->end()){
			if(hdr) {
				rout << "HEvent" << i << "-";
				switch (j) {
					case 0:
						rout << "Time\t";
						break;
					case 1:
						rout << "Source\t";
						break;
					case 2:
						rout << "Sink\t";
						break;
					case 3:
						rout << "Prop\t";
						break;
					case 4:
						rout << "RelSize\t";
						break;
					case 5:
						rout << "NewGrRate\t";
						break;
					case 6:
						rout << "MMid\t";
						break;
					default:
						cerr << "Error in rejector! Historical event header has illegal entry." << endl;
						break;
				}
				
				
				
				j++;
			}
			else rout << *hhsp << "\t";
			
			
			hhsp++;
		}
		hsp++;
		i++;
	}
	
	
	
	
	
	rout << "\t";
	
	if(hdr){
		
		rout << "Test(Alpha)->\t";
		map<string, double>::const_iterator alphp = pt->alphalist.begin();
		if(pt->alphalist.size() == 1 && alphp->first == "ALL"){ //special case for ALL print same alpha under each test
			map<string,string>::const_iterator racc = o->resultacc.begin(); 
			while(racc != o->resultacc.end()){
				rout << racc->first << "(" << alphp->second << ")\t";
				racc++;
			}
		}
		
		else {
			if(pt->alphalist.size() != o->resultacc.size()){
				cerr << "Error in rejector. List of alphas and list of tests not the same size! Aborting!" << endl;
				
				map<string, double>::const_iterator alphp = pt->alphalist.begin();
				while(alphp != pt->alphalist.end()){
					cout << alphp->first << endl;
					alphp++;
				}
				
				map<string,string>::const_iterator racc = o->resultacc.begin(); 
				while(racc != o->resultacc.end()){
					cout << racc->first << endl;
					racc++;
				}
				
				abort();
			}
			
			
			map<string,string>::const_iterator racc = o->resultacc.begin();
			
			
			while(racc != o->resultacc.end()){
				rout << racc->first << "(" << alphp->second << ")\t";
				alphp++;
				racc++;
			}
		}
		
	}
	
	
	else{
		rout << "\t";
		map<string,string>::const_iterator racc = o->resultacc.begin(); 
		
		while(racc != o->resultacc.end()){
			rout << racc->second << "\t";
			racc++;
		}
		
	}
	
	
	rout << "\t";
	
	if(hdr){
		
		rout << "Test(Alpha)->\t";
		map<string, double>::const_iterator alphp = pt->alphalist.begin();
		if(pt->alphalist.size() == 1 && alphp->first == "ALL"){ //special case for ALL print same alpha under each test
			map<string,string>::const_iterator racc = o->resultacc.begin(); 
			while(racc != o->resultacc.end()){
				rout << racc->first << "(" << alphp->second << ")\t";
				racc++;
			}
		}
		
		else {
			if(pt->alphalist.size() != o->resultacc.size()){
				cerr << "Error in rejector. List of alphas and list of tests not the same size! Aborting!" << endl;
				abort();
			}
			
			
			map<string,string>::const_iterator racc = o->resultacc.begin();
			
			
			while(racc != o->resultacc.end()){
				rout << racc->first << "(" << alphp->second << ")\t";
				alphp++;
				racc++;
			}
		}
		
	}
	
	
	
	else{
		rout << "\t";
		map<string,int>::const_iterator nacc = o->numacc.begin(); 
		
		while(nacc != o->numacc.end()){
			map<string,int>::const_iterator ntot = o->numtot.find(nacc->first);
			if(ntot != o->numtot.end()) rout << (double)nacc->second/(double)ntot->second << "\t";
			else rout  << "ERR\t";
			//if(nacc->second == -999)  rout  << "/\t";
			//else rout << nacc->second << "\t";
			nacc++;
		}
		
	}
	
	
	//Stats
	if(keepstats){
	
	if(hdr){
		rout << "******\tStats->\t";
		printstatsinrejout(rout, true);
	}
	
	else{
		rout << "\t\t";
		printstatsinrejout(rout, false);
		//end else if hdr
	}
	
}
	
	
	rout << endl;
	

	
	
}


//---------------------------------------------
//Function Name: printstats
//Function Definition: prints the summary stat values unto the current .rej.txt file
void world::printstatsinrejout(ostream &ksout, bool hdr){
	
	
	
	
	
	
	statmap::const_iterator avt = stats.begin();
	
	
	
	if(hdr){

		
		
		while(avt != stats.end()){
			
			
			
			sresult::const_iterator tp = avt->second.result.begin();
			
			while(tp != avt->second.result.end()){
				
				map<string, map<string, vector<double> > >::const_iterator tq = tp->second.begin();
				int outside = 1;
				while(tq != tp->second.end()){
					
					map<string, vector<double> >::const_iterator tr = tq->second.begin();
					int inside = outside+1;

					while(tr != tq->second.end()){
						
						vector<double>::const_iterator ts = tr->second.begin();
											int locus = 0;
						
						while(ts != tr->second.end()){
							
							if (avt->second.type == "Twoway" || avt->second.type == "TwowayXloci" || avt->second.type == "ByTwoDeme"){
								
								ksout << avt->first << "-" << tp->first << "-" << tq->first << "\t";
								
							}
							else if (avt->second.type == "TwoLoci") {
								
								ksout << avt->first << "-" << tp->first << "-" << tr->first << "<->" << locus << "\t";
								
							}

							
							else if (avt->second.type == "PerLocusXPop") {
								
								ksout << avt->first << "-"  << locus << "\t";
								
							}
							
							else if (avt->second.type == "PerLocusXPopAv") {
								
								ksout << avt->first   <<  "\t";
								
							}
							
							
							else if (avt->second.type == "Oneway") {
								
								//ksout << avt->first << "-"<< outside  << "\t";
								ksout << avt->first << "-"<< tq->first << "-" << locus  << "\t";
								
							}
							
							else if (avt->second.type == "OnewayAv" || avt->second.type == "PerSystem") {
								
								//ksout << avt->first << "-"<< outside  << "\t";
								ksout << avt->first << "-"<< tq->first  << "\t";
								
							}
							

							
							else {
								
								//ksout << avt->first << "-"<< outside  << "\t";
								ksout << avt->first << "-"<< tq->first << "-" << locus  << "\t";
								
							}
							
							locus++;
							ts++;
						}
						inside++;
						
						tr++;
					}
					outside++;
					tq++;
				}
				
				tp++;
			}
			avt++;
		}
		
		ksout << endl;
		
		//end if firsttime
		
		avt = stats.begin(); //reset so first results actually print!
	}
	
	
	
	

	else{
	
	//lets do the numbers
	
	while(avt != stats.end()){
		
		sresult::const_iterator tp = avt->second.result.begin();
		
		//cout << "--------" << endl << avt->first << endl;
		
		while(tp != avt->second.result.end()){
			//cout << "tp: " << tp->first << endl;
			
			map<string, map<string, vector<double> > >::const_iterator tq = tp->second.begin();
			
			while(tq != tp->second.end()){
				//cout << "tq: " << tq->first << endl;
				
				map<string, vector<double> >::const_iterator tr = tq->second.begin();
				
				while(tr != tq->second.end()){
					//cout << "tr: " << tr->first << endl;
					
					vector<double>::const_iterator ts = tr->second.begin();
					
					while(ts != tr->second.end()){
						//cout << "ts: " << *ts << endl;
						
						ksout << *ts << "\t";
						
						ts++;
					}
					tr++;
				}
				tq++;
			}
			tp++;
		}
		avt++;
	}
		
	}
	

	
	
	//end fxn
	
	
}



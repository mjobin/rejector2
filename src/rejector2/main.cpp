//------------------------------------------------
//
// REJECTOR2
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropology
// Stanford University, 2009
//
//------------------------------------------------


#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <set>
#include <map>
#include <cmath>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include "rdeme.h"
#include "CmdLine.h"
#include "rejector.h"
#include "MersenneTwister.h"
#include "public.h"
#include <unistd.h>




#define MAX_OUTFILE_SIZE 1000000

pthread_mutex_t mutex1, mutex2 = PTHREAD_MUTEX_INITIALIZER;
pthread_attr_t attr;

map<int, gridthread> xthreads;


int runs = 0;
int runlimit = 10000;
int brokenruns = 0;
int successcount = 0;
int cpuinuse = 0;
int initoutfilenum = 0; //starting outfile number
int outfilenum = 0; //number of outfile ie infilke.rej0.txt or infile rej1.txt
bool nswitch = false;
bool keepsnpfile = false;
bool printhdr = true;
bool arpread = false;
int simtype = 0; //0= simcoal2, 1=msHOT, 2=MaCS
bool genotypic = false; //if true will try to run with genotypic, not haplotypic data
int comparetype = 0; //only 0 for now simple rejection
bool acdc = false;
bool allornothing = true;
bool alwaysprintoutputline = false;
bool checkinfile = true;
bool suppressoutput = false;
bool keepstats = false;
string copyoutfiles = ""; //if empty copying never done
double rejseed = 0.0;



using namespace std;




//int main (int argc, char * const argv[]) {
int main(int argc, char **argv) {
	bool averagingrun = false;
	int averagingruns = 1000;
	bool averagedoutfile = false;
	bool multiproc = false;
	bool statsonly = false;

	

	
	
	string baseFile;
	string buf, abuf;
	
	
	
	int numcpus = 1; // will be altered by input. test value only
	
	
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
	
	
	
	time_t tstart, tend; //start and stop time of program run
	time (&tstart);
	
    cout << "Rejector" << endl;
	
	
	string inFile;
	
	//Command line parsing
	CCmdLine cmdLine;
	
	//Parse the line, and go to default case if no args given
	if (cmdLine.SplitLine(argc, argv) < 1){
		cout << "What is the name of the data file for the real (experimental) data?" << endl;
		cin >> inFile;
		cout << "What is the limit of number of runs?" << endl;
		cin >> buf;
		runlimit = atoi(buf.c_str());
	}
	
	//look for a supplied filename
	if(cmdLine.HasSwitch("-f")) {
		try {
			inFile = cmdLine.GetArgument("-f", 0);
		}
		catch (...) {
			
			cerr << "-f switch used, but no filename supplied! Aborting.";
			abort();
		}
	}
    
    //look for a supplied filename
	if(cmdLine.HasSwitch("-ao")) {
		arpread = true;
        try {
        
			inFile = cmdLine.GetArgument("-ao", 0);
		}
		catch (...) {
			
			cerr << "-ao switch used, but no filename supplied! Aborting.";
			abort();
		}
	}
	
	//Check if a runlimit is given, otherwise default to 10,000
	runlimit = atoi( cmdLine.GetSafeArgument( "-r", 0, "10000" ).c_str());
	
	//check if tests are to be accepted in any are in tolerance, or if all are
	if(cmdLine.HasSwitch("-s")) {
		allornothing = false;
	}
	
	//If -snp switch then don't kill sc2's modded -snp files
	if(cmdLine.HasSwitch("-snp")) {
		keepsnpfile = true;
	}
	
	//check if there will be an averaged outfile
	if(cmdLine.HasSwitch("-o")) {
		averagedoutfile = true;
	}
	
	//check if this will just be an averaging run, then get number iof runs, defaulting to 1000
	if(cmdLine.HasSwitch("-a")) {
		averagingrun = true;
		averagingruns = atoi( cmdLine.GetSafeArgument( "-a", 0, "1000" ).c_str());
		runlimit = averagingruns;
	}
	
	//This switch will print out one line per iteration of the simulation whether or not the iteration has been accepted
	if(cmdLine.HasSwitch("-p")) {
		alwaysprintoutputline = true;
	}
	
	if(cmdLine.HasSwitch("-rs")) {
		statsonly = true;
	}
	
	if(cmdLine.HasSwitch("-ks")) {
		keepstats = true;
	}
	
	//
	if(cmdLine.HasSwitch("-ms")) {
		simtype = 1;
	}
	
	//check if tests are to be accepted in any are in tolerance, or if all are
	if(cmdLine.HasSwitch("-macs")) {
		simtype = 2;
	}
	
	//check if tests are to be accepted in any are in tolerance, or if all are
	if(cmdLine.HasSwitch("-so")) {
		suppressoutput = true;
		cout.rdbuf(0);
	}
	
	
	//check if the user wants to strat numbering the outfiles with a different number than0, and get that number
	if(cmdLine.HasSwitch("-n")) {
		try {
			initoutfilenum = atoi(cmdLine.GetArgument( "-n", 0).c_str());
			outfilenum = initoutfilenum;
			nswitch = true;
		}
		catch (...) {
			
			cerr << "-x switch used, but could not find cpu number or hostname! Aborting.";
			abort();
		}
	}
	
	//check if this is an multiprocessor run
	if(cmdLine.HasSwitch("-m")) {
		multiproc = true;
		
		try {
			numcpus = atoi(cmdLine.GetArgument( "-m", 0).c_str());
		}
		catch (...) {
			
			cerr << "-m switch used, but could not find cpu number! Aborting.";
			abort();
		}
	}
	
	
	
	if(cmdLine.HasSwitch("-se")) {
		
		try {
			rejseed = atof(cmdLine.GetArgument( "-se", 0).c_str());
		}
		catch (...) {
			
			cerr << "-s switch used, but could not find seed! Aborting.";
			abort();
		}
	}
	
	
	//Abort if both a multiprocessor and an xgrid run specified. That makes no sense.
	if(multiproc && averagingrun){
		cerr << "Alert! Do not use both the -a and -m switch. They are mutually exclusive! Aborting.";
		abort();
	}
	
	
	//check if outfiles are to be copied to a separate location upon completion, and if so where
	if(cmdLine.HasSwitch("-c")) {
		try{
			copyoutfiles  = cmdLine.GetArgument( "-c", 0);
		}
		catch (...) {
			cerr << "-c switch used, but no destination directory given as argument. Aborting.";
			abort();
		}
	}
	
	
	
	//Check if genotypic data is requested in the sim
	if(cmdLine.HasSwitch("-g")) {
		genotypic = true;
		//But for now, this only works if using sc2
		/*
		if(simtype != 0){
			cerr << "Cannot yet use any sim other than simcoal2 with genotypic data." << endl;
			abort();
		}
		*/
	}
	
	//Print number of runs
	if(averagingrun){
		cout << "Running and averaging run for " << averagingruns << " runs." << endl << endl; 
	}
	else if(statsonly){
		cout << "Running to generate statistics only." << endl << endl; 
	}
    else if(arpread){
		cout << "Running to convert .arp file." << endl << endl; 
	}
	else {
		cout << "Running with a limit of " << runlimit << " runs." << endl << endl; 
	}
	
	
	//cout << "Type anything and hit return to GO!" << endl;
	//cin >> buf;
	
	
	
	
	//READ IN INPUT FILE HERE
	
	//load real data file
	
	ifstream in(inFile.c_str());
	if(!in) {
		cout << "Cannot find input file.\n";
		return -1;
	}
	
	//CREATE WORLD AND READ IN DATA
	world r;
    priordist pt;
    
    if(arpread){
        int arp = r.readarpfile(in);
        if(arp){
            cerr << "ERROR: Failed to read .arp file" << endl;
            abort();
        }
        pt = r.blankpriordist();
        
    }
    else pt = r.readinfile(in);
    
	in.close();
	
	/*
	ofstream lout("locusmemtest.txt");
	//TEST SIZE
	//int memtot = r.locinfosizechk(lout);
	int memtot = r.worldsizechk(lout);
	cout << "MEMORY TOTAL: " << memtot << " BYTES" << endl;
	 */
	
	//GET base file name to use for output
	baseFile = inFile;
	if(baseFile.find_last_of(".")){//if there is a period in the string
		baseFile.erase(baseFile.find_last_of("."));
	}
	
	
	
	
	
	
	
	//time setup
	timeval tvsc;
	//Stamp time
	gettimeofday(&tvsc, NULL);
	double bfstat = tvsc.tv_sec + (tvsc.tv_usec/1000000.0);
	
	
	//SEED RNG
	long seed;
	if(rejseed != 0.0)seed = rejseed;
	else seed = tvsc.tv_usec;
	cout << "Seeding with: " << -seed << endl;
	long lidum=-seed;
	ran3(&lidum);
	lidum=1;
	
	
	//cout << "Type anything and hit return to GO!" << endl;
	//cin >> buf;
	
	cout << "****************************TIMING LOOP***************************" << endl;
	
	//RUN STATS ON REAL DATA
	if(pt.testlistall()){
		r.doalltests();
	}
	else{
		r.dotests(&pt, true);
	}
	

	
	gettimeofday(&tvsc, NULL);
	double afstat = tvsc.tv_sec + (tvsc.tv_usec/1000000.0);
	double stattimetaken = afstat-bfstat;
	
	
	if(!arpread && !statsonly){
        cout << "*******************************TIMING SIM*******************************" << endl;
        

        
        
        //RUN PRIOR AND COAL SIM FOR TIME TRIAL
        
        //Stamp time
        gettimeofday(&tvsc, NULL);
        double bfsim = tvsc.tv_sec + (tvsc.tv_usec/1000000.0);
        
        int runchk = timingrun(&r, &pt, baseFile);
        if(runchk == -1){
            cerr << "Error in rejector: Timing sim run failed." << endl;
            abort();
        }
        

        
        gettimeofday(&tvsc, NULL);
        double afsim = tvsc.tv_sec + (tvsc.tv_usec/1000000.0);
        double simtimetaken = afsim-bfsim;
        
        cout << endl << "Estimated time for each run of simulation: " << simtimetaken << endl << endl;
        
        
        //Approximate one iteration
        double eachtime = stattimetaken + simtimetaken; //rough estimate ... time taken for sim and stats for this machine
        cout << "Estimated time for each iteration: " << eachtime << endl << endl;
        
        //Destroy TT. Don't need it anymore
	
    }
	
	cout << "**************************************************************" << endl;
	
	//SWAP IN AVERAGED RESULTS IF REQUESTED
	if(averagedoutfile){
		string rsoutFile = baseFile;
		rsoutFile += "-av.out"; //THIS ALLOWS THE AVERAGED OUTFILE TO BE SUBSTITUED AT THIS STAGE FOR TESTING, BUT REJECTOR STILL NEEDS A SUMMARY FILE

		ifstream oin(rsoutFile.c_str());
		if(!oin) {
			cerr << "Rejector: Cannot find " << rsoutFile << " generated by rejstats!";
			return -1;
		}
		
		r.readoutfile(oin);
		oin.close();
		
	}
	
    
	
	
	//****************************************************************	
	//---- AVERAGING RUN
	//****************************************************************	
	if(averagingrun) {
		
		int runchk = avrun(&r, &pt, baseFile, averagingruns, tstart);
		if(runchk == -1){
			cerr << "Error in rejector: Averaging run failed." << endl;
		}
		
		
	}
	
	else if (multiproc){
		
		cerr << "-m Not supported yet. Aborting." << endl;
		abort();
		
		/*
		 int runchk = multirun(&r, &pt, baseFile, tstart, eachtime, numcpus);
		 if(runchk == -1){
		 cerr << "Error in rejector: Multi run failed." << endl;
		 }
		 */
		
		
	}
	
	else if (statsonly){ //Just print here and end
		
		string outFile = inFile;
		if(outFile.find_last_of(".")){//if there is a period in the string
			outFile.erase(outFile.find_last_of("."));
		}
		
		outFile += ".out";
		ofstream out(outFile.c_str());
		
		r.printalltests(out, inFile);
		
	}
    
    else if (arpread){ //Just print here and end
        string outFile = inFile;
		if(outFile.find_last_of(".")){//if there is a period in the string
			outFile.erase(outFile.find_last_of("."));
		}
		
		outFile += ".txt";
		ofstream out(outFile.c_str());
		
		r.outputinfile(out, pt);
    }
	
	
	else { //the base assumption is a one processor run
		
		
		int runchk = onerun(&r, &pt, baseFile, tstart);
		if(runchk == -1){
			cerr << "Error in rejector: One CPU run failed." << endl;
		}
		
		
	}
	
	
	
	
	//END up to runlimit

	
	
	
	time (&tend);
	double timedif = difftime(tend,tstart);
	cout << endl << "The whole run took " << timedif << " seconds to complete." << endl;
	
	if(!statsonly && !arpread){
		
		
		if(brokenruns > 0) {
			cout << brokenruns << " did not complete due to errors." << endl;
			cout << runlimit - brokenruns << " did complete successfully." << endl;
		}
		else {
			cout << runlimit << " total iterations." << endl;
		}
	}
	
	
	
	
	return 0;
	//****************************************************************
	//END MAIN
	//****************************************************************	
}


int timingrun(world *r, priordist *pt, string baseFile){
//RUN CHOSEN SIM HERE
	
	//** For adding the -n number to the parfile in the sim
	stringstream pf;
	pf << baseFile;
	if(nswitch) pf << initoutfilenum;
	string parFile = pf.str();
	//**
	
world tt = runsim(r, pt, false, parFile);
	
	if(tt.broken){
		cerr << "Error! Simulaiton run failed during time trial." << endl;
		return -1;
		
	}
	return 0;
}

//****************************************************************
//AVERAGING RUN FXN
//****************************************************************	
int avrun(world *r, priordist *pt, string baseFile, int averagingruns, time_t tstart){
	cout << "AVERAGING MODE." << endl;
	
	//Start and stop times for touching files
	time_t touchstart, touchend; 
	time (&touchstart);
	
	string buf; 
	bool failedloop = false;
	
	list<statmap> avlist;
	runlimit = averagingruns;
	
	
	stringstream pf;
	pf << baseFile;

	if(nswitch) pf << initoutfilenum;
	
	string parFile = pf.str();
	
	for(int i=0; i<runlimit; i++){
		
		failedloop = false;
		

		
		
		world s = runsim(r, pt, false, parFile);
		if(s.broken) failedloop = true;
		
		else{
			if(pt->testlistall()){
				s.doalltests();
			}
			else{
				s.dotests(pt, false);
			}
			
		}
		
		avlist.push_back(s.stats);
		
		/*
		ofstream lout("slocusmemtest.txt");
		//TEST SIZE
		//int memtot = r.locinfosizechk(lout);
		int memtot = s.worldsizechk(lout);
		cout << "MEMORY TOTAL: " << memtot << " BYTES" << endl;
		 */
		
		if (failedloop) brokenruns++;
		runs++;
		
		 marktime(runs, runlimit, tstart);
		
		//Touch files so nasty awful sysadmin software doesn't delete it
		time(&touchend);
		double touchdif = difftime(touchend, touchstart);
		if(touchdif > 6000) {
			
			touchstuff(initoutfilenum, outfilenum, parFile);
			time (&touchstart);
			
			
		}
		
		
		//END up to runlimit

	}
	
	
	//PRINT INDIVIDUAL STATS IN LIST
	string avdoutfile = baseFile;
	avdoutfile += "-avdist.out";
	
	//ofstream avdout(avdoutfile.c_str());
	ofstream avdout(avdoutfile.c_str(), ios_base::app);
	if(!avdout) {
		cerr << "Cannot open " << avdoutfile << "in rejector!";
		abort();
	}
	avdistout(avlist, avdout);
	
	/*
	//AVERAGE
	statmap avout =  averagestatresults(avlist);
	
	//PRINT THE AVERAGE
	
	//REPLcing with average, but the quitting!
	r->stats = avout;
	
	string rsoutfile = baseFile;
	rsoutfile += "-av.out";
	
	ofstream rsout(rsoutfile.c_str());
	if(!rsout) {
		cerr << "Cannot open " << rsoutfile << "in rejector!";
		abort();
	}
	
	string inFile = baseFile + ".txt";
	
	r->printalltests(rsout, inFile);
	 */
	
	
	//END up to runlimit

	
	//END AVRUN
	return 0;
}


//****************************************************************
//MULTI RUN FXN
//****************************************************************	
int multirun(world *r, priordist *pt, string baseFile, time_t tstart, double eachtime, int numcpus){
	cout << "MULTIPROCESSOR MODE. " << numcpus << " CPUs." << endl;
	
	//Start and stop times for touching files
	time_t touchstart, touchend; 
	time (&touchstart);
	
	
	
	//figure out if outfgiles large enough
	int runspercpu = runlimit/numcpus;
	cout << "Each output file will be " << runspercpu << " maximum size (esp. if -p invoked.)." << endl;
	if(runspercpu > MAX_OUTFILE_SIZE){
		cerr << "Error: For multiprocessor run, each thread could have an output size larger than the defined MAX_OUTFILE_SIZE." << endl;
		cerr << "Either drop the total runs until each thread gets less than " << MAX_OUTFILE_SIZE << " runs or increase MAX_OUTFILE_SIZE." << endl;
		abort();
	}
	
	//get remander and if uneven add 1 to the first tghread
	int firstrunpercpu = runspercpu;
	if (runlimit%numcpus != 0) firstrunpercpu++;
	   
	
	
	//Now send each thread it's total chunk and let each run alone except for marktime!!!

	
	
	
	
	//Initialize counting variables
	int threadnum = 0;
	map<int, gridthread>::iterator xp = xthreads.begin();
	

	

	
	//Create data to be carried into each thread
	thread_data intd;
	intd.basefile = baseFile;
	intd.tstart = tstart;
	intd.r = r;
	intd.pt = pt;

	
	
	//Guesstimate time to sleep
	double sleeptime = eachtime;
	
	bool firstrun = true;
	
	//Run until all threads dead
	while(xthreads.size() > 0 || firstrun){
		
		firstrun = false;
		
		//If we haven't maxed out the number of threads, make a new one
		if(xthreads.size() < numcpus){
			intd.xid = threadnum;
			if(threadnum == 0) intd.threadrunlimit = firstrunpercpu;
			else intd.threadrunlimit = runspercpu;
	

			
			//SPAWN THREAD
			gridthread newthread;
			xthreads[threadnum] = newthread;
			if(!xthreads[threadnum].multithreadinit(intd)){
				cerr << "Could not create thread number " << threadnum << endl;
			}
			threadnum++;
			
			//End make new thread
		}
		
		
		//SLEEP THIS THREAD SO ITS NOT ALWAYS POLLING
		if(sleeptime > 1) {
			unsigned int sleepint = sleeptime;
			sleep(sleepint);
		}
		else {
			double sleepu = sleeptime*1000000;
			usleep(sleepu);
		}
		
		
		
		xp = xthreads.begin();
		while(xp != xthreads.end()){
			if(xp->second.isdead()){
				xthreads.erase(xp);
				//cout << "ERASED! number of active threads is now " << xthreads.size() << endl;
				break;
			}
			xp++;
		}
		
		
		
		//Touch files so nasty awful sysadmin software doesn't delete it
		time(&touchend);
		double touchdif = difftime(touchend, touchstart);
		if(touchdif > 6000) {
			
			touchstuff(initoutfilenum, outfilenum, baseFile);
			time (&touchstart);
			
			
		}
		
		//END whil
	}
	
	//END MULTIRUN
	return 0;
}


//****************************************************************
//ONE CPU RUN FXN
//****************************************************************	
int onerun(world *r, priordist *pt, string baseFile, time_t tstart){
	bool failedloop = false;
	
	cout << "ONE CPU MODE. " << endl;
	
	//Start and stop times for touching files
	time_t touchstart, touchend; 
	time (&touchstart);
	

	
	while(runs < runlimit){
		
		//Open rejector outfile
		
		stringstream ro;
		ro << baseFile << ".rej" << outfilenum << ".txt";
		string rejoutFile = ro.str();
		ofstream rout(rejoutFile.c_str());
		

			

		//
		
		//outfiletop(rout, inFile, runlimit);  ADD THIS IN LATER
		
		successcount = 0;
		
		printhdr = true;
		
		if(printhdr){
			//Print top
			rout << "Rejector2 Output" << endl;
			rout << "Infile: " << baseFile << ".txt" << endl;
			rout << "Run Limit: " << runlimit << endl;
			time_t fstart;
			time (&fstart);
			struct tm * timeinfo = localtime ( &fstart );
			rout << "Printing started " << asctime (timeinfo) << endl << endl;
			
			
		}
		
		
		//TEST
		//r->printalltests(cout, "R");
		
		while(successcount < MAX_OUTFILE_SIZE && runs < runlimit){ //nice excel-sized chunks
			failedloop = false;
			
			
			//cout << endl << endl << "****************************  " << runs << "   ***************************" << endl;
			
			//** For adding the -n number to the parfile in the sim
			stringstream pf;
			pf << baseFile;
			if(nswitch) pf << initoutfilenum;
			string parFile = pf.str();
			//**

			
			world s = runsim(r, pt, false, parFile);
			if(s.broken) failedloop = true;
			
			else{

				if(pt->testlistall()){
					s.doalltests();
				}
				else{
					s.dotests(pt, false);
				}
				
				//TEST
				//s.printalltests(cout, "S");
				
				//rout << "Compare results." << endl;
				
				if((compareworlds(r, &s, pt, rout, printhdr, keepstats)) == -1) failedloop = true;
				printhdr = false; //After first time through, turn header pritning off
				
			}
			

			

			
			
			runs++;
			if (failedloop) brokenruns++;
			
			//Inform user of how many runs are complete and estimate time to compleetion
			marktime(runs, runlimit, tstart);
			
			
			
			//Touch files so nasty awful sysadmin software doesn't delete it
			time(&touchend);
			double touchdif = difftime(touchend, touchstart);
			if(touchdif > 6000) {
				
				touchstuff(initoutfilenum, outfilenum, baseFile);
				time (&touchstart);
				
				
			}
			
			
			//End while successcout < max_outfile_size
		}
		rout.close();
		outfilenum++;
		//END while runs < runlimit
	}
	
	
	//END ONERUN
	return 0;
}


//MULTIPROCESSOR THREAD
void *multi_thread(void *threadarg){
	struct thread_data *targ;
	targ = (struct thread_data *) threadarg;
	bool failedloop = false;
	
	//Make infile name for printing on new headers
	stringstream rooi;
	rooi << targ->basefile << ".txt";
	string rejinFile = rooi.str();
	
	stringstream tfm;
	tfm << targ->basefile << targ->xid;
	string threadfilename = tfm.str();
	

	
	//Initialize control variables
	successcount = 0;
	int threadruns = 0;
	bool printhdr = true;
	
	

	//Start the outfile
	stringstream ro;
	ro << targ->basefile << ".rej" << targ->xid << ".txt";
	string rejoutFile = ro.str();
	
	
	
	if(printhdr){
		//Print top
		ro << "Rejector2 Output" << endl;
		ro << "Infile: " << targ->basefile << ".txt" << endl;
		ro << "Run Limit: " << runlimit << endl;
		time_t fstart;
		time (&fstart);
		struct tm * timeinfo = localtime ( &fstart );
		ro << "Printing started " << asctime (timeinfo) << endl << endl;
		
		
	}
	
	
	
	
	
	//EACH THREAD runs on its own number of runs
	
	while(threadruns < targ->threadrunlimit){
		failedloop = false;
		
		cout << "Iteration: " << runs << "  on thread " << targ->xid << " and my threadfilename is : " << threadfilename << endl;
		
		world s = runsim(targ->r, targ->pt, false, threadfilename);
		if(s.broken) failedloop = true;
		
		
		if(targ->pt->testlistall()){
			s.doalltests();
		}
		else{
			s.dotests(targ->pt, false);
		}
		
		
		if((compareworlds(targ->r, &s, targ->pt, ro, printhdr, keepstats)) == -1) failedloop = true;
		printhdr = false; //After first time through, turn header pritning off
		
		
		pthread_mutex_lock(&mutex1);
		if(runs < runlimit){ 
			

			
			cout << "I am thread " << targ->xid  << " and i completed run: " << runs << endl;
			runs++;
			if(failedloop) brokenruns++;
			marktime(runs, runlimit, targ->tstart);
			
		}
		pthread_mutex_unlock(&mutex1);
		
		threadruns++;
		
		/*
		world s = runsim(targ->r, targ->pt, false, threadfilename);
		if(s.broken) failedloop = true;
		
		//pthread_mutex_lock(&mutex2);
		if(targ->pt->testlistall()){
			s.doalltests();
		}
		else{
			s.dotests(targ->pt, false);
		}
		//pthread_mutex_unlock(&mutex2);
		
		
		
		//pthread_mutex_lock(&mutex1);
		if(runs < runlimit){ //Check to make sure other threads have not exceeded runlimit before 
			
			
			
			//DUMMY LINE FOR NOW 2/9/09
			if(true) {
				
				stringstream ro;
				ro << targ->basefile << ".rej" << outfilenum << ".txt";
				string rejoutFile = ro.str();
				ofstream rout(rejoutFile.c_str(), ios_base::app);
				
				if(!rout) {
					cerr << "Cannot open " << rejoutFile << "! Leaving thread early!";
					return 0;
				}
				
				cout << "I am thread " << targ->xid  << " and i completed run: " << runs << endl;
				successcount++;
				
				rout.close();
				
				if(successcount >= MAX_OUTFILE_SIZE){  // the outfile is full, advance to the next
					outfilenum++;
					successcount = 0;
					printhdr = true;
					
					stringstream roo;
					roo << targ->basefile << ".rej" << outfilenum << ".txt";
					string rejooutFile = roo.str();
					ofstream roout(rejooutFile.c_str());
					
					
					outfiletop(roout, rejinFile, runlimit);
					
				}
				
				//END Dummy if true
			}
			
			runs++;
			if(failedloop) brokenruns++;
			//marktime(runs, runlimit, targ->tstart);
			//END IF runs< runlimit
		}
		//pthread_mutex_unlock(&mutex1);	
		
		*/
		
		//END while runs<runlimit
	}
	
	
	/*
	 
	 struct thread_data *targ;
	 int syschk;
	 
	 
	 
	 
	 
	 
	 
	 
	 //Make sure everything named by thread number!
	 stringstream dn;
	 dn << targ->xid;
	 string tnum = dn.str();
	 
	 //Get and prepare the basic filename
	 string basename = targ->basefile;
	 if(basename.find_last_of("/") != string::npos){
	 basename.erase(0, (basename.find_last_of("/")+1));
	 
	 }
	 //cout << "basename " << basename << endl;
	 //cout << "I am thread " << targ->xid  << endl;
	 
	 
	 //Get rid of the  files in case some were left over from a previous run
	 string basefileforkill = targ->basefile;
	 basefileforkill += tnum;
	 killfiles(basefileforkill);
	 
	 //EACH THREAD is marking the global number of runs
	 while(runs < runlimit){
	 failedloop = false;
	 //PRIOR
	 string priorinFile = targ->basefile;
	 priorinFile += ".sum";
	 string priorcmd = "./prior -f ";
	 priorcmd += priorinFile;
	 priorcmd += " -o ";
	 priorcmd += tnum;
	 
	 
	 //cout << "Prior command for thread " << targ->xid << " " << priorcmd << endl;
	 
	 
	 syschk = system (priorcmd.c_str());
	 if (syschk==-1){
	 cerr << "Error executing prior on thread " << pthread_self() << "." << endl;
	 failedloop = true;
	 }
	 
	 
	 
	 //SIMCOAL
	 string simcoalinFile = targ->basefile;
	 simcoalinFile += tnum;
	 simcoalinFile += ".par";
	 
	 string simcoalcmd = "./simcoalrej2_1_2 ";
	 simcoalcmd += simcoalinFile;
	 simcoalcmd += " 1 0"; //one run, genotypic data
	 
	 //cout << "Simcoal command for thread " << targ->xid << " " << simcoalcmd << endl;
	 
	 
	 syschk = system (simcoalcmd.c_str());
	 if (syschk==-1){
	 cerr << "Error executing simcoal on thread " << pthread_self() << "." << endl;
	 failedloop = true;
	 }
	 
	 
	 //REJSTATS
	 string scmdpr = targ->basefile;
	 string scmd =  scmdpr;
	 scmd += tnum;
	 scmd +="_0.arp";
	 
	 string srejstatscmd = "./rejstats -sl ";
	 srejstatscmd += scmd;
	 srejstatscmd += " ";
	 srejstatscmd += priorinFile;
	 
	 //cout << "Rejstats command for thread " << targ->xid << " " << srejstatscmd << endl;
	 
	 
	 syschk = system (srejstatscmd.c_str());
	 if (syschk==-1){
	 cerr << "Error executing rejstats on thread " << pthread_self() << "." << endl;
	 failedloop = true;
	 }
	 
	 
	 
	 string rsoutfile = basename;
	 rsoutfile += tnum;
	 rsoutfile += "_0.out";
	 //cout << "Rejstats sim outfile " << rsoutfile << endl;	 
	 
	 
	 string rsparfile = basename;
	 rsparfile += tnum;
	 rsparfile += ".par";
	 //cout << "Parfile for thread " << tnum << " " << rsparfile << endl;	 
	 
	 
	 
	 if(runs < runlimit){ //Check to make sure other threads have not exceeded runlimit before 
	 pthread_mutex_lock(&mutex1);
	 //LOAD SIMULATED DATA FILE and PARAMETER FILE
	 rtestr s;
	 ifstream sin(rsoutfile.c_str());
	 if(!sin) {
	 cerr << "Cannot find " << rsoutfile << "." << endl;
	 failedloop = true;
	 }
	 else {
	 
	 s.readoutfile(sin);
	 
	 
	 //ENTER LINE
	 
	 
	 //cout << "runs: " << runs << endl;	
	 if(comparertestrs(targ->rptr, &s, alphalist) || alwaysprintoutputline) {  // Print an output line if the test accepts OR if instructed to print a line for each iteration
	 
	 ifstream sparin(rsparfile.c_str());
	 if(!sparin) {
	 cerr << "Cannot find " << rsparfile << "! Leaving thread early!";
	 return 0;
	 }
	 
	 s.readparams(sparin);
	 sparin.close();
	 
	 stringstream ro;
	 ro << targ->basefile << ".rej" << outfilenum << ".txt";
	 string rejoutFile = ro.str();
	 ofstream rout(rejoutFile.c_str(), ios_base::app);
	 
	 if(!rout) {
	 cerr << "Cannot open " << rejoutFile << "! Leaving thread early!";
	 return 0;
	 }
	 
	 s.outfileprint(rout, false, alphalist);
	 
	 successcount++;
	 
	 rout.close();
	 
	 if(successcount >= MAX_OUTFILE_SIZE){  // the outfile is full, advance to the next
	 outfilenum++;
	 successcount = 0;
	 printhdr = true;
	 }
	 
	 //END IF ACCEPTS
	 }
	 //END check runlimit before grabbing lock
	 //cout << targ->xid  << " ";
	 
	 //END else for reading in the outfile
	 }
	 
	 
	 runs++;
	 //cout << "I am thread " << targ->xid  << " and i completed run: " << runs << endl;
	 pthread_mutex_unlock(&mutex1);		
	 
	 if(failedloop) brokenruns++;
	 
	marktime(runs, runlimit, targ->tstart);
	 }
	 
	 
	 
	 
	 
	 
	 //Get rid of the  files in case of a problem with simcoal
	 
	 
	 killfiles(basefileforkill);
	 //END RUNS LOOP
	 }
	 
	 
	 //Clean up
	 finalkillfiles(basefileforkill);
	 
	 
	 */
	
	//HIT THE BUZZER ON THE WAY OUT THE DOOR
	xthreads[targ->xid].markdead();
	return 0;
	//END MULTIPROCESSOR THREAD
}




int gridthread::multithreadinit(thread_data t){
	threadalive = true;
	tdata.basefile = t.basefile;
	tdata.xid= t.xid;
	tdata.tstart = t.tstart;
	tdata.r = t.r;
	tdata.pt = t.pt;
	tdata.threadrunlimit = t.threadrunlimit;
	
	time (&threadstart);
	
	if(pthread_create( &xthread, &attr, &multi_thread, &tdata)) return 0;
	else return -1;
}


bool gridthread::isdead(){
	if(!threadalive)return true;
	return false;
}



void gridthread::markdead(){
	threadalive = false;
}

time_t gridthread::retstarttime(){
	return threadstart;
}

void gridthread::cancelthread(){
	pthread_cancel(xthread);	
}

int gridthread::whoami(){
	return tdata.xid;
}



//TIMING FUNCTION
void marktime(int runs, int runlimit, time_t tstart){
	
	time_t tnow;
	time (&tnow);
	double timenowdiff = difftime(tnow,tstart);
	double timeperrunsofar = timenowdiff/runs;
	double secstogo = (runlimit-runs) * timeperrunsofar;
	string timeword = "seconds";
	double timetogo = 0.0;
	
	if(secstogo > 86400) {
		timetogo = secstogo / 86400;
		timeword = "days";
	}
	else if(secstogo > 3600) {
		timetogo = secstogo / 3600;
		timeword = "hours";
	}
	else if(secstogo > 60) {
		timetogo = secstogo / 60;
		timeword = "minutes";
	}
	else {
		timetogo = secstogo;
	}
	
	
	cout.precision(4);
	cout << "Completed " << runs << "/" << runlimit << " iterations so far. " <<  "Estimate " << timetogo << " " << timeword << " to go. \r" << std::flush;
	//cout << "Completed " << runs << "/" << runlimit << " iterations so far. " <<  "Estimate " << timetogo << " " << timeword << " to go. \r" << endl;
	cout.precision(6);
	
	
}



//Touching function because dammit
void touchstuff(int initoutfilenum, int outfilenum, string baseFile){
	
	int syschk=0;
	
	syschk = system ("touch ./prior");
	if (syschk==-1) cerr << "Error touching prior" << endl;
	
	syschk = system ("touch ./rejstats");
	if (syschk==-1) cerr << "Error touching rejstats" << endl;
	
	syschk = system ("touch ./rejector");
	if (syschk==-1) cerr << "Error touching rejector" << endl;
	
	syschk = system ("touch ./simcoalrej2_1_2");
	
	
	syschk = system ("touch ./msHOT");
	
	
	for(int touchi=initoutfilenum; touchi<=outfilenum; touchi++){
		stringstream toucho;
		toucho << baseFile << ".rej" << touchi << ".txt";
		string rejtouch = toucho.str();
		string touchcmd = "touch ";
		touchcmd += rejtouch;
		
		syschk = system (touchcmd.c_str());
		if (syschk==-1) cerr << "Error touching " << rejtouch << endl;
		
	}
	
}

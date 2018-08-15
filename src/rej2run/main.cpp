#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <set>
#include <map>
#include <cmath>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>
#include "CmdLine.h"
#include "rej2run.h"
#include <unistd.h>

pthread_attr_t attr;
map<int, gridthread> xthreads;

using namespace std;


int main(int argc, char **argv) {
	stringstream os;
	string inthing;
	
	string runlimit;
	int numcpus;
	int initoutfilenum = 0;
	
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
	
	os << "./rejector2 ";
	
	//Command line parsing
	CCmdLine cmdLine;
	
	//Parse the line, and go to default case if no args given
	if (cmdLine.SplitLine(argc, argv) < 1){
		cerr << "Error. You must at least use the -m and -f switches." << endl;
		abort();
	}
	
	//look for a supplied filename
	if(cmdLine.HasSwitch("-f")) {
		try {
			inthing = cmdLine.GetArgument("-f", 0);
			os << "-f " << inthing << " ";
		}
		catch (...) {
			
			 		
			cerr << "-f switch used, but no filename supplied! Aborting.";
			abort();
		}
	}
	else {
		cerr << "Error. You must at least use the -m and -f switches." << endl;
		abort();
		
	}
	
	//Check if a runlimit is given, otherwise default to 10,000
	runlimit = cmdLine.GetSafeArgument( "-r", 0, "10000");
	os << "-r " << runlimit << " ";
	
	//check if this is an multiprocessor run
	if(cmdLine.HasSwitch("-m")) {
		
		try {
			numcpus = atoi(cmdLine.GetArgument( "-m", 0).c_str());
		}
		catch (...) {
			
			cerr << "-m switch used, but could not find cpu number! Aborting.";
			abort();
		}
	}
	
	else {
		cerr << "Error. You must at least use the -m and -f switches." << endl;
		abort();
		
	}

	
	
	
	
	
	
	//check if tests are to be accepted in any are in tolerance, or if all are
	if(cmdLine.HasSwitch("-s")) {
		os << "-s ";
	}
	
	//check if there will be an averaged outfile
	if(cmdLine.HasSwitch("-o")) {
		os << "-o ";
	}
	

	
	//This switch will print out one line per iteration of the simulation whether or not the iteration has been accepted
	if(cmdLine.HasSwitch("-p")) {
		os << "-p ";
	}
	
	if(cmdLine.HasSwitch("-rs")) {
		os << "-rs ";
	}
	
	if(cmdLine.HasSwitch("-ks")) {
		os << "-ks ";
	}
	
	if(cmdLine.HasSwitch("-ks")) {
		os << "-g ";
	}
	
	//
	if(cmdLine.HasSwitch("-ms")) {
		os << "-ms ";
	}
	
	//check if tests are to be accepted in any are in tolerance, or if all are
	if(cmdLine.HasSwitch("-macs")) {
		os << "-macs ";
	}
	
	//check if tests are to be accepted in any are in tolerance, or if all are
	if(cmdLine.HasSwitch("-so")) {
		//nothing. Always going to use -so
	}
	
	
	//check if the user wants to strat numbering the outfiles with a different number than0, and get that number
	if(cmdLine.HasSwitch("-n")) {
		try {
			initoutfilenum = atoi(cmdLine.GetSafeArgument( "-n", 0, "0").c_str());

		}
		catch (...) {
			
			cerr << "-n switch used, but could not find number! Aborting.";
			abort();
		}
	}
	

	
	
	
	if(cmdLine.HasSwitch("-se")) {
		
		try {
			inthing = cmdLine.GetArgument( "-se", 0);
			os << "-se " << inthing;
		}
		catch (...) {
			
			cerr << "-s switch used, but could not find seed! Aborting.";
			abort();
		}
	}
	
	

	

	
	
	os << "-so -n ";
	string baseline = os.str();
	
	
	
	
	
	//----RUUUn
	
	//Initialize counting variables
	int threadnum = 0;
	map<int, gridthread>::iterator xp = xthreads.begin();
	
	
	
	
	
	//Create data to be carried into each thread
	thread_data intd;
	

	
	
	
	while(xthreads.size() < numcpus){
		intd.xid = threadnum;
		
		
		
		stringstream to;
		to << baseline << initoutfilenum;
		string outline = to.str();
		
		cout << outline << endl;
		

		intd.outline = outline;
		

		
		cout << "Spawning " << threadnum << endl;
		
		
		//SPAWN THREAD
		gridthread newthread;
		xthreads[threadnum] = newthread;
		if(!xthreads[threadnum].multithreadinit(intd)){
			cerr << "Could not create thread number " << threadnum << endl;
		}
		threadnum++;
		initoutfilenum++;
			
		
		//End make new thread
	}
	

	
	bool firstrun = true;
	
	//Run until all threads dead
	while(xthreads.size() > 0 || firstrun){
		
		firstrun = false;
		

		
		//cout << "sleep" << endl;
		//SLEEP THIS THREAD SO ITS NOT ALWAYS POLLING
		sleep(10);
		
		
		
		xp = xthreads.begin();
		while(xp != xthreads.end()){
			if(xp->second.isdead()){
				xthreads.erase(xp);
				cout << "Number of active threads is now " << xthreads.size() << endl;
				break;
			}
			xp++;
		}
		
		
		

		
		//END whil
	}
	
	
    return 0;
}


int gridthread::multithreadinit(thread_data t){
	threadalive = true;
	tdata.outline = t.outline;
	tdata.xid= t.xid;
	
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


void gridthread::cancelthread(){
	pthread_cancel(xthread);	
}

int gridthread::whoami(){
	return tdata.xid;
}


//MULTIPROCESSOR THREAD
void *multi_thread(void *threadarg){
	struct thread_data *targ;
	
	targ = (struct thread_data *) threadarg;
	
	/*
	for (int i =0; i<1000; i++) {
		cout << targ->xid << " ";
		//	cout << i << " : I am thread " << targ->xid << " and my line is " << targ->outline  << std::flush;
	}
	 */
	
	cout << "I am thread " << targ->xid << " and my line is " << targ->outline  << endl;
	
	
	 //cout << "running on thread " << pthread_self() << endl;
	int syschk = system (targ->outline.c_str());
	if (syschk==-1){
		cerr << "Error executing rejector2 on thread " << pthread_self() << " leaving thread early!" << endl;
		xthreads[targ->xid].markdead();
		return 0;
		
	}
	 
	
	
	
	xthreads[targ->xid].markdead();
	
	return 0;
	//END MULTIPROCESSOR THREAD

}

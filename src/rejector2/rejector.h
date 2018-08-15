//------------------------------------------------
//
// REJECTOR
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropological Sciences
// Stanford University, 2006
//
//------------------------------------------------

#ifndef rejector_h
#define rejector_h

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <list>
#include "rdeme.h"



using namespace std;

extern bool allornothing;

typedef map<string, map<string, map<string, vector<double> > > > tresult;

typedef map<string, map<string, map<string, vector<vector<double> > > > > avtresult;


				   




//---------------------------------------------
//Class Name:
//Class Definition: For rejector: object to contain results
class rtestr {
public:
	map<string, tresult> results;
	map<string, string> testtype;
	double timetaken;
	map<string, string> resultacc; //each rtestr has its own results block now 3/16/06 MJJ
	map<string, int> numacc;
	map<string,int> numtot;
	rtestparams params;
	
	//for real deme only.. stores summary info
	vector<string> demelist;
	string vnaught;
	vector<string> syslist, typelist, mutlist, anclist, locinumlist;
	int numsys;
	
	
	//member fxns
	bool stralphaTest(const char* Str) const;
	bool anydigits(const char* Str) const;
	void readoutfile(ifstream &in);
	void readsumfile(ifstream &in);
	void readparams(istream &in);
	//void outfileprint(ostream &rout, bool hdr, priordist* pt);
	void rsoutfile(ofstream &out, string inFile);
	void onewaytestprint(ostream &out, tresult results, string name);
	void twowaytestprint(ostream &out, tresult results, string name);
	void persystestprint(ostream &out, tresult results, string name);
	void twolocitestprint(ostream &out, tresult results, string name);
	void onewayavtestprint(ostream &out, tresult results, string name);
	void twowayxlocitestprint(ostream &out, tresult results, string name);
	void perlocusxpoptestprint(ostream &out, tresult results, string name);
	void perlocusxpopavtestprint(ostream &out, tresult results, string name);

	//tresult outputresult(string name) { return results->second[name];}
};


struct thread_data{
	string basefile;
	int xid;
	int threadrunlimit;
	time_t tstart;
	world* r; //pointer to original world
	priordist* pt; // priordistribution 
	
	
	
};


class gridthread{
	bool threadalive;
	time_t threadstart;
	thread_data tdata;
	pthread_t xthread;

	
public:
		
	
	//member fxns
	int gridthreadinit(thread_data t);
	int multithreadinit(thread_data t);
	bool isdead();
	void markdead();
	time_t retstarttime();
	int retunfinishedruns();
	void cancelthread();
	int whoami();
};

//-------------------
//FUNCTIONS NOT IN CLASSES
//NEW
void outfiletop(ostream &rout, string inFile, int runlimit);

//RUN FXNS
int timingrun(world *r, priordist *pt, string baseFile);
int avrun(world *r, priordist *pt, string baseFile, int averagingruns, time_t tstart);
int multirun(world *r, priordist *pt, string baseFile, time_t tstart, double eachtime, int numcpus);
int onerun(world *r, priordist *pt, string baseFile, time_t tstart);
statmap averagestatresults(list<statmap> avlist);
void avdistout(list<statmap> avlist, ostream &dout);

//SIM FXNS
double dv8(double v1, double v2, string dtype, bool posit, bool zeroup);
//rtestparams makescinfile(world *r, priordist *pt, ofstream &pout, bool poplock);
//string makemsline(world *r, priordist *pt, bool poplock);
world runsim(world *r, priordist *pt, bool poplock, string basefile);
int compareworlds(world* r, world* s, priordist* pt, ostream &rout, bool hdr, bool keepstats);
rejresult rejectcompare(world *r, world *s, priordist* pt);
bool rejectcomparestats(world* r, world* s, double alpha, string testname, rejresult* o);



//OLD
//map<string, tresult> averagetresults(list<rtestr> avlist);
//void avdistout(list<rtestr> avlist, ostream &dout);
bool comparertestrs(rtestr *r, rtestr *s, map<string, double> alpha);
//bool compareresults(const tresult *r, const tresult *s, double alpha);
bool compareresults(rtestr *r, rtestr *s, double alpha, string testname);
//int compareresults(const tresult *r, const tresult *s, double alpha);
bool istresultempty(const tresult *r);
string xsimcoal (string baseFile, string xgridcmd, string xgridprefix);
bool xsimcoalchk(string xgridcmd, string scid);
string xrejstats (string baseFile, string xgridcmd, string priorinFile, string xgridprefix, string scid);
void xdeljob (string xgridcmd, string delid);
void xresults (string xgridcmd, string xid);
void *multi_thread(void *threadarg, world *r, priordist *pt, string baseFile);
//void *multi_thread(void *threadarg);
void *xgrid_thread(void *threadarg);
void marktime(int runs, int runlimit, time_t tstart);
void killscfiles(string basefile);
void finalkillfiles(string basefile);
void touchstuff(int initoutfilenum, int outfilenum, string baseFile);


#endif



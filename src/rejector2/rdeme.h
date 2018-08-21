//------------------------------------------------
//
// REJECTOR2
// 
// Part of the Rejector rejection analysis software
// by Matthew Jobin, Department of Anthropology
// Stanford University, 2009
//
//------------------------------------------------

#ifndef deme_h
#define deme_h

#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>
#include "MersenneTwister.h"


#include <deque>

#define LDP_ITERATIONS 1000


//some tests will give errors if this is set lower than 2!
#define DEME_CUTOFF 2 

//For runs of homozygosity
#define MIN_HOMOZYGOSITY_RUN 25
#define MIN_HOMOZYGOSITY_LENGTH 1000000
#define MIN_HOMOZYGOSITY_SAMPLE_PROP 0.95
#define MIN_SLIDING_WINDOW_LENGTH 5000000
#define MIN_SNPS_IN_WINDOW 50
#define MAX_HETZ_IN_WINDOW 1
#define MAX_MISS_IN_WINDOW 5
#define MIN_HOM_WIND_THRESH 0.05
#define LONGBLOCK_THRESHOLD 3000000 //and convert based on 1cM/Mb


#define MAX_LOCI_JUMP 500000
#define ADJ_RECOMB_PROB 0.00000001 //chance of recombination between adjacent base pairs over generation




//GLOBALS
using namespace std;
extern bool alwaysprintoutputline;
extern int successcount;
extern int simtype;
extern bool genotypic;
extern int comparetype;
extern bool acdc;


typedef map<string, map<string, map<string, vector<double> > > > sresult; 
typedef vector<double> histevent; 
typedef vector<vector<double> > migmatrix; 




struct statresult {
	string type;
	sresult result; 
};

typedef map<string, statresult> statmap;


//---------------------------------------------
//Class Name:
//Class Definition: sortying fucntion used to keep segmentpointersin maps iorganized
class lociCmp : public binary_function<string, string, bool> {
public:
    bool operator()(const string& x, const string& y) const {
        if(strIntTest(x.c_str()) && strIntTest(y.c_str())) { //they are both integers
            return atoi(x.c_str()) < atoi(y.c_str());
        }
        else {
            return x<y;
        }
    }
    bool strIntTest(const char* Str) const;
};


//Class Name:
//Class Definition: one window
struct slidingwindow {
	int start;
	int stop;
	bool hom;
};

//---------------------------------------------
//Class Name:
//Class Definition: locus data for once sc2 block
class sc2block {
public:
	string type;
	double mutrate;
	double recrate;
	int length;
	int numinblock;
};

//Class Name:
//Class Definition: locus data for once sc2 block
class sc2indchrom {
public:
	vector<sc2block> blocks;
};

//---------------------------------------------
//Class Name:
//Class Definition: a pair of values used for defining a prior distribution
class distn {
public:
	string name; //kind of distn
	double v1;
	double v2;
};

//---------------------------------------------
//Class Name:
//Class Definition: a pair of values used for defining a prior distribution - for migration matricews
class distnmm {
public:
	string name; //kind of distn
	vector<vector<double> > v1;
	vector<vector<double> > v2;
};

//---------------------------------------------
//Class Name:
//Class Definition: a pair of values used for defining a prior distribution - for historical events
class distnhe {
public:
	string name; //kind of distn
	vector<double> v1;
	vector<double> v2;
};



//---------------------------------------------
//Structure Name: Prior
//Class Definition: object containing all the prior distribution parameters
struct priordist {
	
	int numpops;
	map<string, int> demesampsizes;
	vector<distn> popsizes; //type of curce used for this prior
	vector<distn> sampsizes; //type of curce used for this prior
	vector<distn> growthrates;
	vector <distnmm> migrmat;
	vector <distnhe> histev;
	distn msrc;
	map<string, double> alphalist; 
	
	
	//Member functions
	bool testlistall();

};


struct rtestparams{
	map<string, double> demesize;
	map<string, int> samplesize;
	map<string, double> growthrate;
	list<vector<vector <double> > > migmat;
	list<vector<double> > histev;
	
};

//---------------------------------------------
//Class Name:
//Class Definition: Structure to hold results of rejection comparison
class rejresult {
public:
	bool anyaccepted;
	bool anyrejected;
	map<string, string> resultacc; //each rtestr has its own results block now 3/16/06 MJJ
	map<string, int> numacc;
	map<string,int> numtot;
	
	bool rejectordecision();
	//void outfileprint(ostream &rout, bool hdr, priordist* pt);
};

//---------------------------------------------
//Class Name:
//Class Definition: a single stretch of DNA, named and with a list of loci
class locinfo {
public:
    string type;
    double mutrate;
	vector<string> ancstate;
	int length;
	int offset;
	double recrate;
	vector<int> lmap;
	double rsquared;
	int phdist; //position on chromosome

	
	

};
	



struct loopav{
	double tot;
	double count;
};



//---------------------------------------------
//Class Name:
//Class Definition: A collection of demes and the attentdant tests results
class world {
public:
	//ATTRIBUTES
	bool summaryprint;
	bool broken;
	
	//DATA
    double Vnaught;
	vector<locinfo> loci;
	char ** data;
	int datalength, datawidth;
	vector<string> nmap; //Referenced by the vec #
	map<string, vector<int> > dmap;
	map<string, map<int, int> > dblocks;
	map<int, int> pblocks;
	map<string, vector<bool> > droh;
	map<string, map <int, vector<bool> > > dslideroh;
	map<string, statresult> stats;
	rtestparams tparams;

	


	

	


	
	//UTILITY AND CALC FXNS
	int locinfosizechk(ostream &out);
	int worldsizechk(ostream &out);
	int hcounter(vector<int>* y, int j);
	double loctot(vector<int>* y, int j);
	double locmean(vector<int>* y, int j) {return loctot(y,j)/y->size();}
	double locvar(vector<int>* y, int j);
	double difnuc(string w, string x);
	double averageheterozygosity(int j);
	double asd(vector<int>* a, vector<int>* b, int j);
	double vx(vector<int>* y,int j);
	bool missingchk(int match, vector <int> locs, int j);
	vector<int> getnomissing(vector<int>* y, vector <int> locs);
	map<string, map <string, double> > twowaycount(vector<string> a, vector<string> b);
	double chisqrld (map<string, map<string, double> > a);
	double capd(map<string, double>* a, map<string, double>* b);
	double dw(map<string, double>* a, map<string, double>* b);
	map<int, int> rsquaredblocks(vector<int> invec);
	void rsquaredblocksbydeme();
	vector<bool> runsofhomozyg(vector<int>* y);
	void rohbydeme();
	map <int, vector<bool> > rohsliding(vector<int>* y);
	void rohslidebydeme();
	set<int> segsites(vector<int>* y);
	int countsegsites(vector<int>* y);
	bool datamatch(int x, int y, int i, bool missingchk);
	double countpairwise(int x, int y);
	double meanpairwise(vector<int>* y);
	string recrateandhotspots(double n0);
    priordist blankpriordist();


	//STATISTICS
	double tc_haplo(vector<int>* y);
	double tc_heterozygosity(vector<int>*, int);
	double tc_allhet(vector<int>*);
	double tc_allhom(vector<int>* y);
	double tc_ndalleles(vector<int>* y, int j);
	double tc_sampsz(vector <int>* y) {return (double) (y->size());}
	double tc_meanl(vector<int>* y, int j);
	double tc_minl(vector<int>* y, int j);
	double tc_maxl(vector<int>* y, int j);
	double tc_range(vector<int>* y, int j) {return ((tc_maxl(y,j)-tc_minl(y,j))+1.0);}
	double tc_nucdiv(vector<int>* y, int j);
	double tc_fst(int j);
	double tc_beta(vector<int>* y, int j);
	double tc_micsatvar(vector<int>* y, int j);
	double tc_derfrac(vector<int>* y);
	double tc_jstat(vector<int>* y, int j);
	double tc_istat(vector<int>* y, int j);
	double tc_delmusqrd(vector<int>* a, vector<int>* b, int j);
	double tc_ubdelmusqrd(vector<int>* a, vector<int>* b, int j);
	double tc_td(vector<int>* a, vector<int>* b, int j);
	double tc_nei(vector<int>* a, vector<int>* b,  string locustype);
	double tc_dsw(vector<int>* a, vector<int>* b,  string locustype);
	double tc_multldp(vector<int>* y, int j, int k);
	double tc_multld(vector<int>* y, int j, int k);
	double tc_ld(vector<int>* y, int j, int k);
	double tc_hamming(vector<int>* y);
	double tc_numhaploblocks(string dname);
	double tc_blocklength(string dname);
	double tc_blockboundary(string dnamex, string dnamey);
	double tc_numroh(string dname);
	double tc_rohcorr(string dx, string dy);
	double tc_froh(string dname);
	double tc_frohvar(string dname);
	double tc_hammingbtwn(string dnamex, string dnamey);
	double tc_sharedhaplo(string dx, string dy);
	double tc_maf(vector<int>* y);
	double tc_numsegsites(string dname);
	double tc_tajimasd(string dname);
	double tc_longsharedbetween(string dx, string dy);
	double tc_longsharedwithin(string d);
    double tc_varrohlength(string dname);


	
	//STAT LOOPS	
	statresult persystemloop (string locustype, string ancdec, double (world::*tc) (vector<int>*));
	statresult perlocusxpoploop (string locustype, string ancdec, double (world::*tc) (int));
	statresult perlocusxpopavloop (string locustype, string ancdec, double (world::*tc) (int));
	statresult onewayloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, int));
	statresult onewayavloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, int));
	statresult modonewayavloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, int));
	statresult twowayloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, vector<int>*, int));
	statresult twowayxlociloop (string locustype, string ancdec, double (world::*tc) (vector<int>*, vector<int>*, string));
	statresult twolociloop (string locustype, string ancdec, double (world::*tc) (vector<int>* , int , int ));
	statresult twolociadjloop (string locustype, string ancdec, double (world::*tc) (vector<int>* , int , int ));
	statresult bydemeloop (string locustype, string ancdec, double (world::*tc) (string));
	statresult bytwodemeloop (string locustype, string ancdec, double (world::*tc) (string, string));
	
	//PRINT LOOPS
	void persystemprint(statresult* sr, ostream &out);
	void perlocusxpopprint(statresult* sr, ostream &out);
	void perlocusxpopavprint(statresult* sr, ostream &out);
	void onewayprint(statresult* sr, ostream &out);
	void onewayavprint(statresult* sr, ostream &out);
	void modonewayavprint(statresult* sr, ostream &out);
	void twowayprint(statresult* sr, ostream &out);
	void twowayxlociprint(statresult* sr, ostream &out);
	void twolociprint(statresult* sr, ostream &out);
	void twolociadjprint(statresult* sr, ostream &out);
	void bydemeprint(statresult* sr, ostream &out);
	void bytwodemeprint(statresult* sr, ostream &out);
	
	
	int ancdecchk(int i);
	vector<int> getancdecvec(vector<int>* orig, bool anc);
	
	string getdata(int i, int lp);
	double getdataasdbl(int i, int lp);
	map<string, int> getdatacounts(vector<int>* y, int j);
	map<string, double> getdatafreqs(vector<int>* y, int j);
	set<string> getalleles(vector<int>* y, int j);
	vector<string> getallelesasvec(vector<int>* y, int j);
	vector<int> getsubset(vector<int>* y, int j, string match);
	vector<int> getintersection(vector<int>* y, vector<int>* z);
	string makehaplo(int indiv, int low, int high);
	
	
	void printstatresult(statresult* sr, ostream &out);
	void outfileprint(ostream &rout, bool hdr, priordist* pt, rejresult* o, bool keepstats);
	
	
	int makescinfile(world *r, priordist *pt, ofstream &pout, bool poplock);
	string makemsline(world *r, priordist *pt, bool poplock, bool macs);
	
		

    //Member functions
	world();
	~world();
	priordist readinfile(ifstream &in);
	int readarpfile(ifstream &in, world* r);
    int readarpfile(ifstream &in);
	//void readinfile(ifstream &in);
	void readinfile(ifstream &in, ifstream &pin);
	int readmsinfile(world* r, string msline);
	int readmacsinfile(world* r, string msline);
	int readmacsunformatted(world* r, string msline);
	int readoutfile(istream &in);
	bool isempty();
	//world(ifstream &in);
    void heterozygosity();
	void heterozygosityav();
	void heterozygosityanc();
	void heterozygositydec();
	void heterozygosityavanc();
	void heterozygosityavdec();
	void delmusquared();
	void delmusquaredanc();
	void delmusquareddec();
    void t_d();
	void t_danc();
	void t_ddec();
    void DeeSW();
	void NeiDist();
	void betaimbalance();
	void betaimbalanceav();
	void betaimbalanceanc();
	void betaimbalancedec();
	void betaimbalanceavanc();
	void betaimbalanceavdec();
	void Fstdist();
	void Fstdistav();
	void multalleleLD();
	void LD();
	void LDadj();
	void meanlength();
	void maxlength();
	void minlength();
	void samplesize();
	void samplesizeanc();
	void samplesizedec();
	void numdiffalleles();
	void numdiffallelesav();
	void multalleleLDp();
	void rangeoflocus();
	void rangeoflocusav();
	void allheterozygotes();
	void allhomozygotes();
	void haplotypes();
	void jstatistic();
	void jstatisticav();
	void modjstatisticav();
	void istatistic();
	void istatisticav();
	void modistatisticav();
	void microsatvariance();
	void microsatvarianceav();
	void derivedfraction();
	void nucleotidediversity();
	void nucleotidediversityav();
	void ubdelmusquared();
	void ubdelmusquaredanc();
	void ubdelmusquareddec();
	void hamming();
	void numhaploblocks();
	void blocklength();
	void blockboundary();
	void numroh();
	//void rohcorr();
	void froh();
	void frohvar();
	void sharedhaploblocks();
	void hammingbtwn();
	void minorallelefreq();
	void tajimasd();
	void numsegsites();
	void longsharedbetween();
	void longsharedwithin();
    void varrohlength();




    void counts(ostream &out);
	void outputinfile(ostream &out, priordist pt);
	void outputsummary(ostream &sout);
	void dotests(priordist* pt, bool timetests);
	void doalltests();
	bool allelenumcheck();
	bool allelematch(set<string> allele, string match);
	bool allelethreshold(int loc, int thresh, bool over);
	bool checksnpdiversity(int thresh, bool over);
	void printalltests(ostream &out, string inFile);
	bool anydigits(const char* Str) const;
	void printstatsinrejout(ostream &ksout, bool hdr);
		
	
	


};

#endif


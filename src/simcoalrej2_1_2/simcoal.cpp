/*//////////////////////////////////////////////////////////////////////////////

Copyright Laurent Excoffier, Genetics and Biometry Lab, University of Geneva
E-mail: laurent.excoffier@anthro.unige.ch
URL: http://anthropologie.unige.ch/~laurent

The use of the following source code in only allowed if it is kept intact
with its original banner

//////////////////////////////////////////////////////////////////////////////*/

#pragma hdrstop
//#include <condefs.h>
                          
#include "cstring.h"
#include "cond_var.h"
#include "genealgy.h"
#include "public.h"
#include "time.h"
#include "migrmat.h"
#include "deme.h"


#include "simcoal.h"

#ifdef  _LINUX_
   #include <sys/stat.h>
   #include <sys/time.h>
   #include <sys/types.h>
   #include <unistd.h>
#else
   #include <Filectrl.hpp>
   #include "dir.h"
#endif


#ifdef _PRINT_ARLEQUIN_OUTPUT_
   ofstream ArlSimulConditions;
#endif

//---------------------------------------------------------------------------
int
TMeanMat::compute_mean(const int& n) {
	_num_updates=n   ;
	if (!_num_updates) return 0;
	for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]/=_num_updates;
      }
   }
   return 1;
}

//Compute the s.d. from the values present in the matrix, assuming that
//those are the sum of square values
int
TMeanMat::compute_sd(const int& n, const TMeanMat& mean) {
   if (n<2) return 0;
	_num_updates=n   ;
	if (!_num_updates) return 0;
	for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	double s2=_data[i][j],
         		 m=mean._data[i][j]; //assumes that mean contains mean values
      	_data[i][j]=sqrt( (s2-n*m*m)/(n-1) );
      }
   }
   return 1;
}

int
TMeanMat::update_with(const TMigrationMatrix& MM) {
	if (!MM.size() || _size!=MM.size()) return 0;
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]+=MM(i,j);
      }
   }return 1;
}

int
TMeanMat::update_with_square(const TMigrationMatrix& MM) {
	if (!MM.size() || _size!=MM.size()) return 0;
   for (int i=0; i<_size; ++i) {
   	for (int j=0; j<_size; ++j) {
      	_data[i][j]+=MM(i,j)*MM(i,j);
      }
   }return 1;
}

//------------------------------------------------------------------------------
void my_strrev(char * str) {
   int size=strlen(str);
   char *buf = new char[size+1];
   for (int i=0; i<size; ++i) {
      buf[i]=str[size-1-i];
   }
   for (int i=0; i<size; ++i) {
      str[i]=buf[i];
   }
   delete[] buf;
}
//------------------------------------------------------------------------------
my_string remove_extension(const my_string& s) {
  	my_string sout;
   char sin1[400], sin2[400];
   strcpy(sin1,s.c_str());
   //strrev(sin1);
   int i, len = strlen(sin1);
   //First: invert my_string;
   my_strrev(sin1);
   for (i=0; i<len ; ) {
      sin2[i]=sin1[i];
   	if (sin1[i]=='.') break;
      ++i;
   }
   sin2[i+1]='\0';
   my_strrev(sin2);
   my_strrev(sin1);
   char* c=strstr(sin1,sin2);
   if (!c) return "";
   c[0]='\0';
   sout=sin1;
   return sout;
}
//------------------------------------------------------------------------------
my_string extract_path(const my_string& s) {
	my_string sout;
   char sin1[400], sin2[400];
   strcpy(sin1,s.c_str());
   int i, len=strlen(sin1);
   //First: invert my_string;
   my_strrev(sin1);
   for (i=0; i<len ; ) {
      sin2[i]=sin1[i];
   	if (sin1[i]=='\\') break;
      ++i;
   }
   sin2[i+1]='\0';
   my_strrev(sin2);
   my_strrev(sin1);
   if (strlen(sin1)==strlen(sin2)) {   //it means that no backslash was detected
   	sout="";
      return sout;
   }
   char* c=strstr(sin1,sin2);
   if (!c) return "";
   c[1]='\0';

   sout=sin1;
   return sout;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

typedef MY_TArrayAsVector<my_string> NameArray;

float UNIT_TIME;



int main(int argc, char *argv[]) {



	my_string usage;

   usage ="\nsim_coale.exe usage:\n";
   usage+="\n [1] : sim_coal <no arguments>";
   usage+="\n        The program will prompt you for an input file name to process and the";
   usage+="\n        number of simulations to perform";
   usage+="\n [2] : sim_coal <batch file name>";
   usage+="\n        The program will open the batch file that should contain, on the first";
   usage+="\n        line, the number of simulations to perform, and a flag to  produce ";
   usage+="\n        either haplotypic (0) or genotypic (1) data";
   usage+="\n        On the following lines,the names of the input files to process";
   usage+="\n        should be listed";
   usage+="\n [3] : sim_coal <file name> <number of simulations to perform>";
   usage+="\n [3]    <haplotypic[0]/genotypic[1] data flag>\n";


   NameArray infile_name(1,5);

   my_string batchFile;
   ifstream ifs_batch;

   my_string Path="";


   int num_rand_samples=1, tot_num_rand_samples;
   int num_files=0;
   int genotypic_data=0;

   if (argc==2) {
   	//Then one supposes that one receives a file containing the name of
      //the files to be processed
      batchFile=argv[1];
      Path=extract_path(batchFile);
      ifs_batch.open(batchFile.c_str());
      if (!ifs_batch) {
        cout << "Unable to open input file (" << batchFile << ")\n";
         cout << usage;
			return -1;
      }
      else {
         	my_string line;
      	//Number of random samples to draw must be on first line
      	ifs_batch >> num_rand_samples >> genotypic_data;
         line.read_to_delim(ifs_batch);
      	while (ifs_batch) {
         	//line.read_to_delim(ifs_batch);
            line.read_line(ifs_batch);
            if (!line.is_null()) {
      			++num_files;
               infile_name.Add(Path+remove_extension(line));
            }
      	}
         ifs_batch.close();
      }
   }
   else
   if (argc==4) {
   	//Then one supposes that one receives a file name to be processed
      //and the number of permutations to be performed
      num_files=1;
      batchFile=argv[1];
      Path=extract_path(batchFile);
      ifs_batch.open(batchFile.c_str());
      if (!ifs_batch) {
      	cout << "Unable to open input file (" << batchFile << ")" << endl;
         cout << usage;
			return -1;
      }
      else {
      	num_rand_samples=atoi(argv[2]);
         genotypic_data=atoi(argv[3]);
         ifs_batch.close();
         infile_name.Add(Path+remove_extension(batchFile));
         if (!num_rand_samples) cout << usage;
      }
   }
   else {
      num_files=1;
      my_string curName;
      cout << "\n\nGeneric input file name  (*.par)  : ";
      cin >> curName;
      infile_name.Add(curName);
      cout << "\n\nNo. of random samples to generate : ";
      cin >> num_rand_samples;
      cout << "\n\nGenotypic [1] or haplotypic [0] data output for Arlequin : ";
      cin >> genotypic_data;
      if (genotypic_data <0 || genotypic_data>1 ) {
         cout << usage;
         return -1;
      }
   }

   //Reads and generate genealogies for the required number of files

   for (int nf=0; nf<num_files; ++nf) {
   	//Creating input and output file names
   	my_string	outfile_name,
      				arl_batch_file_name,
               	arl_generic_name,
      				paup_batch_file_name,
               	paup_generic_name,
               	paup_name, log_file_name,
               	paup_true_trees_name, paup_mut_trees_name,
                  simul_params_file_name;;


		outfile_name			= infile_name[nf]+".gen";
      arl_batch_file_name 	= infile_name[nf]+".arb";
      arl_generic_name		= infile_name[nf];
      paup_generic_name		= infile_name[nf];
      paup_batch_file_name = infile_name[nf]+".bat";
      paup_name				= infile_name[nf]+".paup";
      paup_true_trees_name	= infile_name[nf]+"_true_trees.trees";
      paup_mut_trees_name	= infile_name[nf]+"_mut_trees.trees";
      log_file_name 			= infile_name[nf]+".log";

      simul_params_file_name = infile_name[nf]+".simparam";

   	//Initialization of a bunch of object and vectors

   	TDemeCollection 	myDemes;
      TEventArray 	   myEvents(10,10);
      TDemeSizeArray  	myDemeSizes(10,10);
      TGrowthArray 		myGrowthRates(10,10);
      TSampSizeArray 	mySampleSizes(10,10);
      TMigrMatArray     myMigrationMatrices(1,10);  //Loro_16_9_99
	  
	

   	//Initialize input file
   	ifstream ifs;

      my_string dir        = infile_name[nf];
      const char * pdir          = dir.c_str();
   	infile_name[nf]+=".par";       
      ifs.open(infile_name[nf].c_str());
      
   	if (!ifs) {
      	cout << "Unable to open input file (" << infile_name[nf] << ")" << endl;
			continue;
   	}

       /*
      #ifdef  _LINUX_
        mode_t mode=0755;

        //if (DirectoryExists(pdir)){
         mkdir(pdir,mode);
        //}
        //if (DirectoryExists(pdir)){
         chdir(pdir);
        //}
      #else
        //char old_dir[MAXDIR],new_dir[MAXDIR];

        if (!DirectoryExists(pdir)){
                if (!mkdir(pdir)){
                throw Exception("Cannot create c:\\temp Results.");
                }
        }
        if (DirectoryExists(pdir)){
                chdir(pdir);
        }
      #endif
       */

	// ***********************
	//edited out so a directory is not made mjj 3/6/06
/*
      #ifdef  _LINUX_
        mode_t mode=0755;

        //if (DirectoryExists(pdir)){
                mkdir(pdir,mode);
        //}
        //if (DirectoryExists(pdir)){
                chdir(pdir);
        //}
      #else
        //char old_dir[MAXDIR],new_dir[MAXDIR];mkdir(pdir);
        mkdir(pdir);
        chdir(pdir);
      #endif
*/

      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions.open(simul_params_file_name.c_str());
      #endif

   	int num_pop;
   	my_string line;

   	line.read_to_delim(ifs);   //Reads a blank line

   	ifs >> num_pop;        //Reads the number of populations to simulate
		line.read_to_delim(ifs);

   	//Read deme sizes
   	line.read_to_delim(ifs);   //Reads a blank line
   		long cur_size;
   	for (int i=0; i<num_pop; ++i) {
      	ifs >> cur_size;
      	line.read_to_delim(ifs);
   		myDemeSizes.Add(cur_size);
   	}

      #ifdef _VERBOSE_
   	cout << "\nDeme sizes\n";
   	for (int i=0; i<num_pop; ++i) {
   		cout << "Deme " << i << "\t" << myDemeSizes[i] << "\n";
   	}
   	cout << "\n";
      #endif


      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions
         << "Simulation parameters for current simulations\n"
         << "=============================================\n\n";
   	ArlSimulConditions << "\nDeme sizes\n";
   	for (int i=0; i<num_pop; ++i) {
   		ArlSimulConditions << "Deme " << i << "\t" << myDemeSizes[i] << "\n";
   	}
   	ArlSimulConditions << "\n";
      #endif


   	//Reads the sample sizes
   	line.read_line(ifs);  //Reads a blank line
      int sum_sample_size=0;
   	for (int i=0; i<num_pop; ++i) {
      	ifs >> cur_size;
         sum_sample_size+=cur_size;
      	line.read_to_delim(ifs);
   		mySampleSizes.Add(cur_size);
   	}

		int tot_num_nodes=0;
      #ifdef _VERBOSE_
   	cout << "\nSample sizes\n";
      #endif

      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions << "\nSample sizes\n";
      #endif

   	for (int i=0; i<num_pop; ++i) {
      	tot_num_nodes+=mySampleSizes[i];
         #ifdef _VERBOSE_
   		cout << "Deme " << i << "\t" << mySampleSizes[i] << "\n";
         #endif
         
         #ifdef _PRINT_ARLEQUIN_OUTPUT_
         ArlSimulConditions << "Deme " << i << "\t" << mySampleSizes[i] << "\n";
         #endif
   	}
      #ifdef _VERBOSE_
   	cout << endl;
      #endif
      
      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions << endl;
      #endif

   	line.read_line(ifs);  //Reads a blank line
   		float cur_growth;
   	for (int i=0; i<num_pop; ++i) {
   		ifs >> cur_growth;
      	line.read_to_delim(ifs);
   		myGrowthRates.Add(cur_growth);
   	}

      #ifdef _VERBOSE_
   	cout << "Growth rates\n";
   	for (int i=0; i<num_pop; ++i) {
   		cout << i << "\t" << myGrowthRates[i] << "\n";
   	}
   	cout << endl;
      #endif
      
      #ifdef _PRINT_ARLEQUIN_OUTPUT_
   	ArlSimulConditions << "Growth rates\n";
   	for (int i=0; i<num_pop; ++i) {
   		ArlSimulConditions << i << "\t" << myGrowthRates[i] << "\n";
   	}
   	ArlSimulConditions << endl;
      #endif

      //Loro_16_9_99
      //Reading how many migration matrices to read
      int numMigMat=0; //Default
   	line.read_line(ifs);  //Reads a blank line
      ifs >> numMigMat;
      line.read_line(ifs);

      //Loro_16_9_99
      if (numMigMat)
      for (int m=0; m<numMigMat; ++m) {
   		//Reading migration matrix
   		TMigrationMatrix myMigrationRates(num_pop);
   		line.read_line(ifs);  //Reads a blank line
   		for (int i=0; i<num_pop; ++i) {
   			for (int j=0; j<num_pop; ++j) {
      			ifs >> myMigrationRates(i,j);
            }
      	}
         myMigrationMatrices.Add(myMigrationRates);
         //	  TMigrationMatrix test_mat(2, STEPPING_STONE_3D, 0.01);

         #ifdef _VERBOSE_
   		cout << "\nMigration matrix\n";
   		cout << myMigrationMatrices[m];
         #endif

         #ifdef _PRINT_ARLEQUIN_OUTPUT_
   		ArlSimulConditions << "\nMigration matrix\n";
   		ArlSimulConditions << myMigrationMatrices[m];
         #endif

      	line.read_to_delim(ifs); //Reads until the end of line
   	}
      else { //Then a model without migration is assumed
      	TMigrationMatrix myMigrationRates(num_pop); //Enough to create a matrix filled with zeroes
         myMigrationMatrices.Add(myMigrationRates);
      }

   	int num_events;
   	//Read historical events
  		line.read_line(ifs); //Reads a blank line
   	ifs >> num_events;   //Reads the number of historical events to read
   	line.read_line(ifs);
   	for (int i=0; i<num_events; ++i) {
   		THistoricalEvent curevent;
   		ifs >> curevent;
   		myEvents.Add(curevent);
   	}

      #ifdef _VERBOSE_
   	cout << "\nHistorical events\n";
   	for (int i=0; i<num_events; ++i) {
   		cout << "Event " << i << "\n" << myEvents[i] << "\n";
   	}
      if (!num_events) cout << "No historical events defined in input file ...\n";
   	cout << endl;
      #endif
                  
      #ifdef _PRINT_ARLEQUIN_OUTPUT_
   	ArlSimulConditions << "\nHistorical events\n";
   	for (int i=0; i<num_events; ++i) {
   		ArlSimulConditions << "Event " << i << "\n" << myEvents[i] << "\n";
   	}
      if (!num_events) ArlSimulConditions << "No historical events defined in input file ...\n";
   	ArlSimulConditions << endl;
      #endif


      //Loro_26_2_99 to avoid a bug when no historical event is defined
      if (num_events) line.read_to_delim(ifs); //Reads until the end of line

      //Loro_19_06_01 : Modification to allow for the computation of independent loci
      //Read number of independent loci
      int num_indep_loci;
      line.read_line(ifs);


      ifs >> num_indep_loci;
      line.read_line(ifs);
      line.remove_blanks();
      //Loro_02_03_04 Add one parameter to allow the simulation of chromosomes
      //with different number and types of loci
      int diff_chrom_struct=0;  //A flag to say if we want to implement loci with
                              //different structure and marker types...


      if ( ! line.is_null()) {
        diff_chrom_struct=atoi(line.c_str());
      }
      //ifs >> diff_chrom_struct;




      #ifdef _VERBOSE_
      cout << "\nNumber of independent loci to simulate : " << num_indep_loci;
      if (diff_chrom_struct) cout << "\nEach with a different chromosomal structure\n";
      else cout << "\nwith the same chromosomal structure\n";
      #endif    

      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions << "\nNumber of independent loci to simulate : " << num_indep_loci;
      if (diff_chrom_struct) ArlSimulConditions << "\nEach with a different chromosomal structure\n";
      else ArlSimulConditions << "\nwith the same chromosomal structure\n";
      #endif

      //line.read_to_delim(ifs); //Reads until the next item

      int num_chrom_structs;
      if (diff_chrom_struct) num_chrom_structs=num_indep_loci;
      else num_chrom_structs=1;

      double gamma_a=0.0;
      Mut_Type data_type;

      double * mut_ratio            =  new double[num_chrom_structs];
      double * prop_sites           =  new double[num_chrom_structs];
      double * num_rate_categories  =  new double[num_chrom_structs];
      int    * num_linked_loci      =  new int   [num_chrom_structs];
      double * rec_rate             =  new double[num_chrom_structs];
  		double * mut_rate             =  new double[num_chrom_structs];
      double transition_rate;
      double rec_rate_perBlock;
  		double mut_rate_perBlock;
  		double geom_param_perBlock;
      int    rangeConstraint;

      TNumLociArray  * rangeConstraintPerLoc    = new TNumLociArray [num_chrom_structs]; //A value of zero assumes no range constraint
      TNumLociArray  * num_linked_loci_perBlock = new TNumLociArray [num_chrom_structs];
      TRecRateArray  * rec_rate_perLoc          = new TRecRateArray [num_chrom_structs];
      TMutRateArray  * mut_rate_perLoc          = new TMutRateArray [num_chrom_structs];
      TMutRateArray  * geom_param_perLoc        = new TMutRateArray [num_chrom_structs];
      TMutRateArray  * transRatePerLoc          = new TMutRateArray [num_chrom_structs];
      double*        * freq_SNP_min             = new double*       [num_chrom_structs];
      TDataTypeArray * data_type_perLoc         = new TDataTypeArray[num_chrom_structs];
      bool           * Data_mix                 = new bool          [num_chrom_structs];
      int            * num_linkage_blocks       = new int           [num_chrom_structs];

      rec_rate_perBlock       =  0;
      rangeConstraint         =  0;
      mut_rate_perBlock       =  0;
      geom_param_perBlock     =  0;
      transition_rate         =  0;

      for (int k=0; k<num_chrom_structs; ++k) {
         mut_ratio                [k]   =  1.0;
         prop_sites               [k]   =  1.0;
         num_rate_categories      [k]   =  0;
         num_linked_loci          [k]   =  0;
         rec_rate                 [k]   =  0;
  			mut_rate                 [k]   =  0;
         Data_mix                 [k]   =  false;
      }

      int num_tot_loci=0;
      int num_tot_interval=0;
      for (int num_structs=0; num_structs<num_chrom_structs; ++num_structs) {
         //Read number of contiguous blocks of markers (linkage block)
         //int num_linkage_blocks;
         line.read_line(ifs);
         ifs >> num_linkage_blocks[num_structs];

         #ifdef _VERBOSE_
         cout << "\nNumber of linkage blocks to simulate in structure " << (num_structs+1) << ": " << num_linkage_blocks[num_structs] << "\n";
         #endif
       
         #ifdef _PRINT_ARLEQUIN_OUTPUT_
         ArlSimulConditions << "\nNumber of linkage blocks to simulate in structure " << (num_structs+1) << ": " << num_linkage_blocks[num_structs] << "\n";
         #endif

         line.read_to_delim(ifs); //Reads until the end of line

         /*
         double gamma_a=0.0, mut_ratio=1.0, prop_sites=1.0, num_rate_categories=0;

         int num_linked_loci=0;
         double rec_rate=0;
         double rec_rate_perBlock;
  			float mut_rate=0;
  			float mut_rate_perBlock;
         float transition_rate=0.0;
         int range_constraint=0; //A value of zero assumes no range constraint
         TNumLociArray  num_linked_loci_perBlock;
         TRecRateArray rec_rate_perLoc;
         TMutRateArray mut_rate_perLoc;
         double* freq_SNP_min= new double[num_linkage_blocks];
         Mut_Type data_type;
         TDataTypeArray data_type_perLoc;
         bool Data_mix=false;
         */

         freq_SNP_min[num_structs] =  new double [num_linkage_blocks[num_structs]];

         line.read_line(ifs);  //Read comments
         for (int i=0, k=0; i<num_linkage_blocks[num_structs]; ++i) {
            //Read number of loci per Block, per generation recombination rate (between two adjacent loci),
            //mutation rate (for the whole block), type of sequence
            line.read_to_delim(ifs, ' ');
            line.to_upper();
            //Loro_01_03_04 Modif to avoid bug...
            int numLinkedLociInBlock;
            ifs >> numLinkedLociInBlock;
            num_linked_loci_perBlock[num_structs].Add(numLinkedLociInBlock);
            if (line.contains("MICROSAT")) {
         	   data_type=MICROSAT;
               ifs
                  >> rec_rate_perBlock
                  >> mut_rate_perBlock
                  >> geom_param_perBlock
                  >> rangeConstraint; //Extract range constraint for microsat data
         	   //line.read_to_delim(ifs);
            }
            else
            if (line.contains("DNA")) {
         	   data_type=DNA;
               ifs >> rec_rate_perBlock >> mut_rate_perBlock;
               ifs >> transition_rate; //Extract transition rate for DNA
         	   //line.read_to_delim(ifs);
         	   /*
         	   //Reads gamma parameter
               ifs >> gamma_a;
         	   //cout << "\nGamma parameter : " << gamma_a << endl;
         	   //If gamma parameter is negative, then it means that we have a
         	   //two mutation rate model. So read the mutation ratio and the proportion
         	   //of sites with the high mutation rate
         	   if (gamma_a<0.0) ifs >> mut_ratio >> prop_sites;
         	   else
         	   if (gamma_a>0.0) ifs >> num_rate_categories;
               */
            }
            else
            if (line.contains("SNP")) {
         	   data_type=SNP;
               ifs >> rec_rate_perBlock >> freq_SNP_min[num_structs][i]; //Extract transition rate for DNA
         	   mut_rate_perBlock=0;
               //line.read_to_delim(ifs);
            }
            else
            if (line.contains("RFLP")) {
               data_type=RFLP;
               ifs >> rec_rate_perBlock >> mut_rate_perBlock;
            }

            for(int j=0; j<num_linked_loci_perBlock[num_structs][i]; ++j) {
               rec_rate_perLoc[num_structs].Add(rec_rate_perBlock);
               mut_rate_perLoc[num_structs].Add(mut_rate_perBlock);
               geom_param_perLoc[num_structs].Add(geom_param_perBlock);
               rangeConstraintPerLoc[num_structs].Add(rangeConstraint);
               transRatePerLoc[num_structs].Add(transition_rate);
               data_type_perLoc[num_structs].Add(data_type);
               mut_rate[num_structs]+=mut_rate_perBlock;
               ++k;
            }
            num_linked_loci[num_structs]+=num_linked_loci_perBlock[num_structs][i];

            line.read_to_delim(ifs); //Reads until the end of line

            #ifdef _VERBOSE_
            cout << "\n   " << num_linked_loci_perBlock[num_structs][i] << " partially linked ";
            if (data_type==MICROSAT) {
               cout << "MICROSAT:\n      recombination, mutation rates, geometric parameter, and range constraint : ";
               cout << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               cout << fixed << setprecision(7) << mut_rate_perBlock << "   ";
               cout << fixed << setprecision(7) << geom_param_perBlock << "   ";
               cout << fixed << setprecision(7) << rangeConstraint;
            }
            else if (data_type==RFLP) {
               cout << "RFLP    :\n      recombination and mutation rates                  : ";
               cout << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               cout << fixed << setprecision(7) << mut_rate_perBlock;
            }
            else if (data_type==DNA) {
               cout << "DNA     :\n      recombination, mutation rates and transition rates :";
               cout << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               cout << fixed << setprecision(7) << mut_rate_perBlock << "   ";
               cout << fixed << setprecision(7) << transition_rate ;
            }
            else if (data_type==SNP) {
               cout << "SNP     :\n      recombination rate and minimum frequency      : ";
               cout << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               cout << fixed << setprecision(7) << freq_SNP_min[num_structs][i] ;
            }
            cout << "\n";

            #endif  //_VERBOSE_ 

            #ifdef _PRINT_ARLEQUIN_OUTPUT_
            ArlSimulConditions << "\n   " << num_linked_loci_perBlock[num_structs][i] << " partially linked ";
            if (data_type==MICROSAT) {
               ArlSimulConditions << "MICROSAT loci:\n      recombination, mutation rates, geometric parameter, and range constraint : ";
               ArlSimulConditions << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << mut_rate_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << geom_param_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << rangeConstraint;
            }
            else if (data_type==RFLP) {
               ArlSimulConditions << "RFLP loci:\n      recombination and mutation rates                  : ";
               ArlSimulConditions << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << mut_rate_perBlock;
            }
            else if (data_type==DNA) {
               ArlSimulConditions << "DNA loci:\n      recombination, mutation rates and transition rates :";
               ArlSimulConditions << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << mut_rate_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << transition_rate ;
            }
            else if (data_type==SNP) {
               ArlSimulConditions << "SNP loci:\n      recombination rate and minimum frequency      : ";
               ArlSimulConditions << fixed << setprecision(7) << rec_rate_perBlock << "   ";
               ArlSimulConditions << fixed << setprecision(7) << freq_SNP_min[num_structs][i] ;
            }
            ArlSimulConditions << "\n";
            #endif  //_PRINT_ARLEQUIN_OUTPUT_
            

            num_tot_loci+=num_linked_loci_perBlock[num_structs][i];
         }
          
        num_tot_interval+=(num_linked_loci[num_structs]-1);

         for(int i=0; i<(num_linked_loci[num_structs]-1); ++i) {
            rec_rate[num_structs]+=rec_rate_perLoc[num_structs][i];
            if(i!=0) {
               if(data_type_perLoc[num_structs][i] != data_type_perLoc[num_structs][i-1] ){
                  Data_mix[num_structs]=true;
               }
            }
         }



      }//End of reading chromosome structure

      if(num_chrom_structs==1){
       num_tot_loci=num_indep_loci*num_linked_loci[0];
       num_tot_interval=num_indep_loci*(num_linked_loci[0]-1);
      }


      #ifdef _PRINT_ARLEQUIN_OUTPUT_
         ArlSimulConditions << endl;
      #endif  //_PRINT_ARLEQUIN_OUTPUT_

      /*************************************************************************

                       Call a new function for recombination
       which set once for all the number of linked loci the recombination rate


      *************************************************************************/

      //loro_27_2_99
      ifs.close();





   	//Opening result files
   	//ofstream ofs(outfile_name.c_str());

      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ofstream ArlBatch(arl_batch_file_name.c_str());
      #endif

      #ifdef _PRINT_SEQ_PAUP_BATCH_
      ofstream PaupBatch(paup_batch_file_name.c_str());
      #endif

      // Loro_28_2_99 : Be careful: only one node allowed to be created before
      // building the tree
      TNode first_node;

      first_node.reset_node_count(); //Do not forget to reset the node count...

      first_node.initialise_first_Node(num_linked_loci[0]);

   	//Initialization of the demes  or of the TRecDemes
      myDemes.create_demes(&myMigrationMatrices);
      myDemes.initialize_deme_size(myDemeSizes);
      myDemes.initialize_growth_rates(myGrowthRates);
      myDemes.initialize_events(myEvents);
      myDemes.create_lineages(mySampleSizes);
      myDemes.set_num_indep_loci(num_indep_loci);

      myDemes.set_num_linked_loci(num_linked_loci[0]);
      myDemes.set_rec_rate(rec_rate[0]);
      myDemes.set_rec_rate_perLoc(&rec_rate_perLoc[0]);
      myDemes.set_mut_rate_perLoc(&mut_rate_perLoc[0]);
      myDemes.set_geom_param_perLoc(&geom_param_perLoc[0]);
      myDemes.set_rangeConstraintPerLoc(&rangeConstraintPerLoc[0]);
      myDemes.setTransitionRatePerLoc(&transRatePerLoc[0]);

      const TNumLociArray* pnum_linked_loci_perBlock=&num_linked_loci_perBlock[0];
      myDemes.set_num_linked_loci_perBlock(pnum_linked_loci_perBlock);
      myDemes.set_freq_SNP_min(freq_SNP_min[0]);
      myDemes.set_num_linkage_blocks(num_linkage_blocks[0]);
      myDemes.set_sum_sample_size(sum_sample_size);
      //myDemes.ind_loci=num_indep_loci;

      #ifdef _COMPUTE_MOMENTS_
         //myDemes.allocate_recombHits();
         myDemes.allocate_recombHits(num_tot_interval);
      #endif

    	//Output initial conditions
   	//ofs 	<< "\nDemes initial conditions"
      //   	<< "\n========================\n"
      //   	<< myDemes
      //   	<< endl;

   	//Initialize the generator of random numbers once for all
      ///////////////////////////////////////////////////////////////////////////
      //Guillaume 24 06 03: I have troubles with the function clock()
      //I prefere using the function time() to initialize the random generator


      //clock_t start=clock();
   	time_t start=clock();
   	//long lidum=-start;
      //ran3(&lidum);
      //lidum=1;
   	//cout << "Random generator initialized with : " << start << endl;
      //ofs  << "Random generator initialized with : " << start << "\n\n";
      //long seed = time(NULL) - ( time(NULL)/1000 ) * 1000;

      #ifdef  _LINUX_
      struct timeval tv ;
      struct timezone tz;
      gettimeofday(&tv, &tz);

      long seed = tv.tv_usec;
      #else
      //long seed = (long) time(NULL);
      long seed = time(NULL) - ( time(NULL)/1000 ) * 1000;
      
      //Debug
      //seed = 530;
      #endif

   	long lidum=-seed;
      ran3(&lidum);
      lidum=1;

      #ifdef _VERBOSE_
   	cout << "\nRandom generator initialized with : " << seed << endl;
      //ofs  << "Random generator initialized with : " << seed << "\n\n";
   	cout << "\nBuilding " << num_rand_samples << " genealogies ...\n";
      #endif

      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions << endl ; 
      ArlSimulConditions << "\nRandom generator initialized with : " << seed << endl;
      #endif


      #ifdef _COMPUTE_MOMENTS_
   	//A matrix to store mean coalescent times within and among populations
   	TMeanMat MeanCoalMat(myDemes.num_demes(), 0.0);
   	//A matrix to store s.d. coalescent times within and among populations
   	TMeanMat sdCoalMat(myDemes.num_demes(), 0.0);
      //A matrix to store mean number of mutations within and among populations
   	TMeanMat MeanPairDiff(myDemes.num_demes(), 0.0);
   	//A matrix to store s.d. number of mutations within and among populations
   	TMeanMat sdPairDiff(myDemes.num_demes(), 0.0);

      //We write results for every chromosome in different files
      ofstream* CoalTimesFile;
      my_string* CoalTimes_out;
      if (num_chrom_structs==1) { //Only compute moments if there is a single chromosome structure...
         CoalTimesFile= new ofstream[num_indep_loci];
         CoalTimes_out= new my_string[num_indep_loci];

         for(int i=0;i<num_indep_loci;++i) {
            int deca, signa, ndiga = 0;
            CoalTimes_out[i]=paup_generic_name + "_Chr" + fcvt(i+1, ndiga, &deca, &signa) + "_CoalTimes.txt";
            CoalTimesFile[i].open(CoalTimes_out[i].c_str());
         }
      }
      double* T_MRCA;
      double* mean_T_MRCA;
      double* Num_recomb_hits;
      double* Num_effective_recomb_hits;
      double* mean_Num_recomb_hits;
      double* mean_Num_effective_recomb_hits;

      int upBoundary;
      int** recomb_distrid;
      int** eff_recomb_distrid;


         int numMAX_stat_computed=10;
         char** list_name_stat= new char*[numMAX_stat_computed];
         double** list_stat_time=new double*[numMAX_stat_computed];
         for(int i=0;i<numMAX_stat_computed;++i) {
            list_name_stat[i]="empty";
            list_stat_time[i]=NULL;
         }

         cout << num_tot_loci << num_tot_interval;

         T_MRCA=new double[num_tot_loci];
         mean_T_MRCA=new double[num_tot_loci];
         Num_recomb_hits                =  new double[num_tot_interval];
         Num_effective_recomb_hits      =  new double[num_tot_interval];
         mean_Num_recomb_hits           =  new double[num_tot_interval];
         mean_Num_effective_recomb_hits =  new double[num_tot_interval];

         upBoundary=1000;
         recomb_distrid=new int*[upBoundary];
         eff_recomb_distrid=new int*[upBoundary];

         for(int i=0; i<upBoundary; ++i) {
            recomb_distrid[i]=new int[num_tot_interval];
            eff_recomb_distrid[i]=new int[num_tot_interval];
            for(int j=0; j<num_tot_interval; ++j) {
               recomb_distrid[i][j]=0;
               eff_recomb_distrid[i][j]=0;
            }
         }
         for(int curLocus=0;curLocus<num_tot_loci;++curLocus) {
            T_MRCA[curLocus]=0;
            //within_loci_mean_time[curLocus]=0;
            //within_loci_sd_time[curLocus]=0;
            //sumTime[curLocus]=0;
            //ProdSumTime[curLocus]=0;

            mean_T_MRCA[curLocus]=0;
            //mean_within_loci_mean_time[curLocus]=0;
            //mean_within_loci_sd_time[curLocus]=0;
            //mean_within_loci_sum_time[curLocus]=0;
            //covTime_first_others[curLocus]=0;
            //Between_loci_mean_T[curLocus]=0;
            //Between_loci_sd_T[curLocus]=0;
            //THudson[curLocus]=0;
         }
         for(int curLocus=0;curLocus<num_tot_interval;++curLocus) {
            Num_recomb_hits[curLocus]=0;
            Num_effective_recomb_hits[curLocus]=0;
            mean_Num_recomb_hits[curLocus]=0;
            mean_Num_effective_recomb_hits[curLocus]=0;
         }

         list_stat_time[0]=T_MRCA;
         //list_stat_time[1]=within_loci_mean_time;
         //list_stat_time[2]=within_loci_sd_time;
         list_stat_time[1]=Num_recomb_hits;
         list_stat_time[2]=Num_effective_recomb_hits;
         //list_stat_time[5]=ProdSumTime;
         //list_stat_time[6]=sumTime;
         //list_stat_time[2]=Between_loci_mean_T;
         //list_stat_time[3]=Between_loci_sd_T;
         //list_stat_time[6]=THudson;
     
      #endif //_COMPUTE_MOMENTS_

      #ifdef _PRINT_SEQ_PAUP_
      ofstream PaupFile;
      if (num_chrom_structs==1) PaupFile.open(paup_name.c_str());
      #endif

      #ifdef _PRINT_TREE_PAUP_
   	//ofstream
      //   PaupTrueTreesFile(paup_true_trees_name.c_str()),
      //   PaupMutTreesFile(paup_mut_trees_name.c_str());
   	ofstream  PaupTrueTreesFile;
      if (num_chrom_structs==1) PaupTrueTreesFile.open(paup_true_trees_name.c_str());
      if (data_type != SNP){
      }
      else{
         //paup_mut_trees_name	= infile_name[nf]+"_mut_SNPtrees.trees";
         paup_mut_trees_name	= " ";
      }
      ofstream  PaupMutTreesFile;
      if (num_chrom_structs==1) PaupMutTreesFile.open(paup_mut_trees_name.c_str());
      #endif //_PRINT_TREE_PAUP_

      #ifdef _PRINT_SEQ_PAUP_
      if (num_chrom_structs==1) {
	      PaupFile
            << "#NEXUS\n"
         	<< "\n[Simulated data generated by program sim_coal.exe (Laurent Excoffier)"
            << "\n\tSimulation conditions:"
            << "\n\tNo. of taxa   : " << tot_num_nodes
            << "\n\tNo. of loci   : " << num_linked_loci
            << "\n\tNo. of demes  : " << num_pop
            << "\n\tMutation rate : " << (mut_rate/num_linked_loci) << " per locus (site) per generation"
            << "\n\tAlpha         : " << gamma_a << " (0=homogeneity)"
            << "\n\t%Transitions  : " << (100*transition_rate) << " (ts/tv=";
         if (transition_rate!=1)
         	PaupFile << (transition_rate/(1-transition_rate));
         else
				PaupFile << "inf";
         PaupFile << ")\n]\n\n";
         write_PAUP_header(PaupFile,log_file_name.c_str());
      }
  	   #endif  //_PRINT_SEQ_PAUP_

      #ifdef _PRINT_TREE_PAUP_
      if (num_chrom_structs==1) {
         write_PAUP_trees_section_header(PaupTrueTreesFile);
         if(data_type != SNP){
            write_PAUP_trees_section_header(PaupMutTreesFile);
         }
      }
  	   #endif //_PRINT_TREE_PAUP_

      //float mut_rate_per_locus=0.0;
      //if (num_linked_loci) mut_rate_per_locus=mut_rate/num_linked_loci;

      Mut_Model mut_mod;

      if (num_chrom_structs==1) {
         if (fabs(gamma_a)>1e-7) mut_mod=K80_GAMMA;
         else mut_mod=K80_NOGAMMA;
      }
      //#endif

      tree_output_type tree_type;


      //Build unequal mutation rates if needed
      //Create the vector of gamma distributed mutation rates once for all

   	//Loro_26_2_99
      if (num_chrom_structs==1) {
   		if (!gamma_a==0.0 && (data_type==DNA || data_type==SNP)) {
            #ifdef _VERBOSE_
         	cout << "\nInitializing unequal mutation rates...";
            #endif
            
            double mjj_nrc = *num_rate_categories; //MJJ 10/1/12
            
			   const int inum_rate_categories= (int) mjj_nrc; //MJJ 10/1/12
   			first_node.mut_rates.get_new_rates( num_linked_loci[0],gamma_a, mut_ratio[0], prop_sites[0],
            											   inum_rate_categories);
            #ifdef _VERBOSE_
         	cout << "done\n\n";
            #endif
   		}
      }

      /* TODO -olaurent -cMutation rates : Compute new mutation rates and put it into a dedicated vector, also allowing for gamma distributed mutation rates */

     	int tot_cases=0;
      my_float tree_mut_length=0.0, tree_mut_length2=0.0,
      			tree_exp_length=0.0, tree_exp_length2=0.0;

      tot_num_rand_samples= num_rand_samples * num_indep_loci;

   	for ( ; tot_cases<num_rand_samples; ++tot_cases) {

         #ifdef _COMPUTE_MOMENTS_
         if (num_chrom_structs==1) {
            for(int curLocus=0;curLocus<num_tot_interval;++curLocus) {
               Num_recomb_hits[curLocus]=0;
               Num_effective_recomb_hits[curLocus]=0;
            }
            myDemes.initialize_recombHits();
         }
         #endif  //_COMPUTE_MOMENTS_
        // cout << "\t\tsample simulated : " << tot_cases << endl;


         #ifdef _PRINT_ARLEQUIN_OUTPUT_
         ofstream ArlInFile;
         #endif

         for (int num_ind_loci=0, curChromStruct=0; num_ind_loci<num_indep_loci; ++num_ind_loci ) {
               
			 //  cout << "num_ind_loci: " << num_ind_loci << " curChromStruct: " << curChromStruct << endl;
			   
			   

			   
            myDemes.ind_loci=num_ind_loci;
            if (diff_chrom_struct) curChromStruct=num_ind_loci;

            /*
            //Loro_02_03_04 potential bug because each microsat loci should be checked for potential range constraint...
            if (num_chrom_structs==1) {
                if (data_type==MICROSAT && range_constraint[curChromStruct]) {//Fix minimum and maximum size for microsats
							double y=2.0, x=range_constraint[curChromStruct];
                   if (fmod(x,y)==0.0) { //Even number
                   	first_node.min_mic=-(range_constraint[curChromStruct]/2)+1;
                      first_node.max_mic=range_constraint[curChromStruct]/2;
                   }
                   else {
                   	first_node.min_mic=(-range_constraint[curChromStruct]+1)/2;
                      first_node.max_mic=(range_constraint[curChromStruct]-1)/2;
                   }
                }
            }
            */

            if(tot_cases==1){
                 int stop=0;
            }


            //if (curChromStruct) {
               if (num_linked_loci[curChromStruct]!=num_linked_loci[curChromStruct-1]) {
                  myDemes.set_num_linked_loci(num_linked_loci[curChromStruct]);
               }
               myDemes.set_rec_rate(rec_rate[curChromStruct]);

               myDemes.set_rec_rate_perLoc(&rec_rate_perLoc[curChromStruct]);
               myDemes.set_mut_rate_perLoc(&mut_rate_perLoc[curChromStruct]);
               myDemes.set_geom_param_perLoc(&geom_param_perLoc[curChromStruct]);
               myDemes.set_rangeConstraintPerLoc(&rangeConstraintPerLoc[curChromStruct]);
               myDemes.setTransitionRatePerLoc(&transRatePerLoc[curChromStruct]);

               pnum_linked_loci_perBlock=&num_linked_loci_perBlock[curChromStruct];

               myDemes.set_num_linked_loci_perBlock(pnum_linked_loci_perBlock);

               myDemes.set_freq_SNP_min(freq_SNP_min[curChromStruct]);
               myDemes.set_num_linkage_blocks(num_linkage_blocks[curChromStruct]);
            //}

            if (myDemes.build_tree()) {

            	//Add mutations
      			//long mut=myDemes.sprinkle_mutations(mut_rate, num_linked_loci,
					//		data_type, gamma_a, mut_ratio, prop_sites, transition_rate,
               //      range_constraint);
               //Loci with different data types

               if(tot_cases==1){
                 int stop=0;
               }

				
				//MJJ 2/1/10
				ofstream snpfs;
				my_string snpout        = dir;
				snpout +="-snp.txt";
				
				/*
				my_string killfile = "rm -f ";
				killfile += snpout;
				int syschk = system (killfile.c_str());
				if (syschk==-1){
					cerr << "Error executing " << killfile << " !";
				}
				 */
				
				snpfs.open(snpout.c_str());
				
				if (!snpfs) {
					cout << "Unable to open snp out file (" << snpout << ")" << endl;
					continue;
				}
				
				
				
				
				//MJJ 2/1/10
				
				
               long mut=myDemes.sprinkle_mutations(     mut_rate[curChromStruct],
                                                        num_linked_loci[curChromStruct],
					 		&data_type_perLoc[curChromStruct],
                                                        gamma_a, mut_ratio[curChromStruct],
                                                        prop_sites[curChromStruct],
                                                        transition_rate,
                                                        0, snpfs);

               //cout << "mutation completed" << endl;
			   


               tree_mut_length+=mut;
               tree_mut_length2+=mut*mut;
               tree_exp_length+=first_node.tree_exp_length();
               tree_exp_length2+=first_node.tree_exp_length()*first_node.tree_exp_length();

               #ifdef _PRINT_SEQ_PAUP_
               if (num_chrom_structs==1) {
							//Writing data to PAUP file
                   write_PAUP_replicate_number((tot_cases+1), PaupFile);
                   write_PAUP_data_header(tot_num_nodes, num_linked_loci, data_type,PaupFile);
                   myDemes.write_samples_to_PAUP_file(PaupFile, data_type);
							write_PAUP_end_matrix(PaupFile);
							write_PAUP_tree_header(PaupFile, (tot_cases+1));
                   tree_type=MUT_RATE;
                   myDemes.print_gene_tree(PaupFile, tree_type, mut_rate_per_locus);
                   write_PAUP_end(PaupFile);
                   if (num_linked_loci) {
                   	PaupFile
                   		<< "\n\t[Tree length = " << (first_node.tree_exp_length()/num_linked_loci)
                         << "\n\tNumber of sites hit by mutations = "  << first_node.count_polym_sites()
                         << "]\n";
                   }
                   write_PAUP_block(PaupFile,mut_mod);
               }
               #endif  //_PRINT_SEQ_PAUP_

               #ifdef _PRINT_TREE_PAUP_
               if (num_chrom_structs==1) {
                   //No Recombination
                   if(rec_rate[curChromStruct]==0){
                      tree_type=GENERATIONS;
                      write_PAUP_tree_name(PaupTrueTreesFile, (tot_cases+1), tree_type);
                      myDemes.print_gene_tree(PaupTrueTreesFile, tree_type, mut_rate[curChromStruct]);
                      tree_type=NUM_MUT;
                      write_PAUP_tree_name(PaupMutTreesFile, (tot_cases+1), tree_type);
                      myDemes.print_gene_tree(PaupMutTreesFile, tree_type, mut_rate[curChromStruct]);
                   }
                   else { //Recombination case
                      for(int curLocus=0; curLocus<num_linked_loci[curChromStruct];++curLocus) {
                         tree_type=GENERATIONS;
                         write_PAUP_tree_name_perLoc(PaupTrueTreesFile, (tot_cases+1), curLocus, tree_type);
                         myDemes.print_gene_tree(PaupTrueTreesFile, tree_type, mut_rate_perLoc[curChromStruct][curLocus], curLocus);
                         if(data_type == SNP){
                         }
                         else{
                             //tree_type=NUM_MUT;   //Guillaume 15 11 2004
                            tree_type=MUT_RATE;     //TO AVOID the storing for every TNode
                                                    //of the number of mutations per every locus

                            write_PAUP_tree_name_perLoc(PaupMutTreesFile, (tot_cases+1), curLocus, tree_type);
                            myDemes.print_gene_tree(PaupMutTreesFile, tree_type, mut_rate_perLoc[curChromStruct][curLocus], curLocus);
                         }
                      }
                   }
               }
               #endif //_PRINT_TREE_PAUP_

               #ifdef _PRINT_SEQ_PAUP_BATCH_
               if (num_chrom_structs==1) {
							//Build paup file name
                   my_string PaupOut, PaupLog;
  							int dec, sign, ndig = 0;
                   PaupOut = paup_generic_name + "_" + fcvt(tot_cases,ndig, &dec, &sign) + ".pau";
                   PaupLog = paup_generic_name + "_" + fcvt(tot_cases,ndig, &dec, &sign) + ".log";
                   //Open paup output file
      					ofstream PaupOutFile(PaupOut.c_str());

                   //Write name in paup batch file
                   PaupBatch << "paup " << PaupOut << "\n";
                   PaupOutFile
            	 		<< "#NEXUS\n"
      						<< "\n[Simulated data generated by program sim_coal.exe (Laurent Excoffier)"
            	 		<< "\n\tSimulation conditions:"
            	 		<< "\n\tNo. of taxa   : " << tot_num_nodes
            	 		<< "\n\tNo. of loci   : " << num_loci
            	 		<< "\n\tNo. of demes  : " << num_pop
            	 		<< "\n\tMutation rate : " << (mut_rate/num_loci) << " per locus (site) per generation"
            	 		<< "\n\tAlpha         : " << gamma_a << " (0=homogeneity)"
            	 		<< "\n\t%Transitions  : " << (100*transition_rate) << " (ts/tv=";
                   if (transition_rate!=1) PaupOutFile << (transition_rate/(1-transition_rate));
      					else	PaupOutFile << "inf";
      					PaupOutFile << ")\n]\n\n";

                   write_PAUP_header(PaupOutFile,PaupLog.c_str());
                   write_PAUP_replicate_number((tot_cases+1), PaupOutFile);
                   write_PAUP_data_header(tot_num_nodes, num_loci, data_type,PaupOutFile);
                   myDemes.write_samples_to_PAUP_file(PaupOutFile, data_type);
							write_PAUP_end_matrix(PaupOutFile);
							write_PAUP_tree_header(PaupOutFile, (tot_cases+1));
                   tree_type=MUT_RATE;
                   myDemes.print_gene_tree(PaupOutFile, tree_type, mut_rate_per_locus);
                   write_PAUP_end(PaupOutFile);
                   	PaupOutFile
                   		<< "\n\t[Tree length = " << (first_node.tree_exp_length()/num_loci)
                         << "\n\tNumber of sites hit by mutations = "  << first_node.count_polym_sites()
                         << "]\n";
							write_PAUP_block(PaupOutFile,mut_mod);
							write_PAUP_footer(PaupOutFile);
               }
      		   #endif  //_PRINT_SEQ_PAUP_BATCH_

               #ifdef _PRINT_ARLEQUIN_OUTPUT_
            	if (!num_ind_loci) {
                  //Build Arlequin file name
                  my_string ArlOut;
						int deca, signa, ndiga = 0;
						char * tot_cases_str = (char *) malloc(10);
						sprintf(tot_cases_str,"%d", tot_cases);
						ArlOut = paup_generic_name + "_" + tot_cases_str + ".arp";
                //  ArlOut = paup_generic_name + "_" + fcvt(tot_cases, ndiga, &deca, &signa) + ".arp";

      				//Open Arlequin output file
                  ArlInFile.open(ArlOut.c_str());

                  //Write name in Arlequin batch file
                  ArlBatch << ArlOut << "\n";

                  ArlInFile
                  	<< "#Arlequin input file written by program the simulation program sim_coal.exe\n\n";

                  /*
                  ArlInFile
      					<< "#Simulation parameters:\n"
                     << "#======================\n"
                     << "\n#Deme sizes\n";

                  for (int i=0; i<num_pop; ++i) {
   						ArlInFile << "#Deme #" << i << "\t" << myDemeSizes[i] << "\n";
   					}
                  ArlInFile << "\n#Sample sizes\n";
   					for (int i=0; i<num_pop; ++i) {
   						ArlInFile << "#Deme #" << i << "\t" << mySampleSizes[i] << "\n";
   					}
                  ArlInFile << "\n#Growth rates\n";
   					for (int i=0; i<num_pop; ++i) {
   						ArlInFile << "#Deme #" << i << "\t" << myGrowthRates[i] << "\n";
   					}
                  for (int i=0; i<numMigMat; ++i) {
                  	ArlInFile
                  		<< "\n#Migration matrix " << i << "\n"
   							<< myMigrationMatrices[i];
                  }
                  ArlInFile << "\n\n#Historical events\n";
   					for (int i=0; i<num_events; ++i) {
   						ArlInFile << "#Event #" << i << "\n" << myEvents[i] << "\n";
   					}
      				ArlInFile
                  	<< "\n#Data type : " ;

                  if (Data_mix){
                     ArlInFile << "Mixture of different data types\n";
                  }
                  else{
                        switch (data_type) {
                  	   case DNA     	: ArlInFile << "DNA\n";
                     					   break;
                  	   case SNP     	: ArlInFile << "SNP (by default used as STANDARD data in arlequin)";
                     					   break;
                        case RFLP      : ArlInFile << "RFLP\n"; break;
                        case MICROSAT  : ArlInFile << "MICROSAT\n"; break;
                     }
                  }

                  if (!num_ind_loci) //only print the following if this is the first independent locus we simulate
                  ArlInFile
                  	//<< setprecision(10)
                  	<< "\n#Mutation rate per generation                     : " << mut_rate[0]
                     << "\n#Number of partially linked loci to simulate      : " << num_linked_loci[0]
                     << "\n#Recombination rate per generation                : " << rec_rate[0]
                     << "\n#Number of indep. set of loci to simulate         : " << num_indep_loci
                   	<< "\n#Gamma parameter                                  : " << gamma_a
                  	//<< "\n#Mutation ratio                                   : " << mut_ratio
                     << "\n\n";
                  
                  */

                  //Loro_20_9_99: Get the number of demes with non-empty samples
                  int num_real_samples=0;
                  for (int d=0; d<num_pop; ++d) {
                  	if (myDemes[d].sample_size()) ++num_real_samples;
                  }

                  if (Data_mix[curChromStruct]){
                     Mut_Type data_type_mix;
                     data_type_mix = STANDARD;
                     write_Arlequin_header(num_real_samples, data_type_mix, ArlInFile, genotypic_data);
                  }
                  else{
                     write_Arlequin_header(num_real_samples, data_type, ArlInFile, genotypic_data);
                  }
               }

               //myDemes.write_samples_to_Arlequin_file(ArlInFile, data_type);
               //myDemes.write_locus_data_to_array(data_type);
               //to output a mix of different data types
               myDemes.write_locus_data_to_array(&data_type_perLoc[curChromStruct]);

               #endif //_PRINT_ARLEQUIN_OUTPUT_

               #ifdef _COMPUTE_MOMENTS_
               if (num_chrom_structs==1) {
            	 //Compute the first two moments of coalescence times and mean number
						//of mutation differences over the whole tree
						//myDemes.compute_moments_of_coalescence_times();
                  //Recombination:
                  //myDemes.compute_moments_of_coalescence_times(list_stat_time);


                  for(int curLocus=0;curLocus<num_linked_loci[0];curLocus++) {
                       T_MRCA[num_ind_loci*num_linked_loci[0] + curLocus]= (double) myDemes.GeneTree.MRCA_list[curLocus]->time;
                  }



                  for(int i=0;i<num_pop;++i) {
                     for(int curLocus=0;curLocus<(num_linked_loci[0]-1);curLocus++) {
                        Num_recomb_hits[num_ind_loci*(num_linked_loci[0]-1) + curLocus] +=
                           myDemes.return_recombHits(myDemes.get_deme(i).get_id()*num_indep_loci*(num_linked_loci[0]-1) + (num_ind_loci*(num_linked_loci[0]-1) + curLocus));

                        Num_effective_recomb_hits[num_ind_loci*(num_linked_loci[0]-1) + curLocus] +=
                        (
                           myDemes.return_recombHits(myDemes.get_deme(i).get_id()*num_indep_loci*(num_linked_loci[0]-1) + (num_ind_loci*(num_linked_loci[0]-1) + curLocus))
                        -
                           myDemes.return_no_effective_recombHits(myDemes.get_deme(i).get_id()*num_indep_loci*(num_linked_loci[0]-1) + (num_ind_loci*(num_linked_loci[0]-1) + curLocus))
                        );

                     }
                  }

            	 //Compute mean coalescence times and mean number of mutation
                  //differences both within and among demes


            	 //myDemes.compute_moments_of_demes_coalescence_times();

             	  //MeanCoalMat.update_with(myDemes.coal_time_mat());
            	 //Adds up the square of the coalescent times
            	 //sdCoalMat.update_with_square(myDemes.coal_time_mat());
                  //MeanPairDiff.update_with(myDemes.mean_pair_diff_mat());
                  //sdPairDiff.update_with_square(myDemes.mean_pair_diff_mat());

                  for(int curLocus=0;curLocus<num_linked_loci[0];++curLocus) {
                     int k=num_ind_loci*num_linked_loci[0];
                     mean_T_MRCA[k + curLocus]+=T_MRCA[k + curLocus];
                     //mean_within_loci_mean_time[curLocus]+=within_loci_mean_time[curLocus];
                     //mean_within_loci_sd_time[curLocus]+=within_loci_sd_time[curLocus];
                     //mean_within_loci_sum_time[curLocus]+=sumTime[curLocus];
                     //covTime_first_others[curLocus]+=ProdSumTime[curLocus];
                     //Between_loci_mean_T[curLocus]+=Between_loci_mean_T[curLocus];
                     //Between_loci_sd_T[curLocus]+=Between_loci_sd_T[curLocus];
                     //THudson[curLocus]+=THudson[curLocus];
                  }

                  for(int curLocus=0;curLocus<(num_linked_loci[0]-1);++curLocus) {
                     int k=num_ind_loci*(num_linked_loci[0]-1);
                     mean_Num_recomb_hits[k + curLocus]+=Num_recomb_hits[k + curLocus];
                     mean_Num_effective_recomb_hits[k + curLocus]+=Num_effective_recomb_hits[k + curLocus];
                  }

                  ofstream ofcoaltime(CoalTimes_out[num_ind_loci].c_str(),ios::app);
                  int num_stat_displayed=1;
                  list_name_stat[0]="T_MRCA   ";
                  //list_name_stat[1]=" wthLociMeanTime";
                  //list_name_stat[2]=" wthLociStdvTime";
                  write_output_coalTimes(ofcoaltime, num_stat_displayed, list_name_stat,
                                         list_stat_time, num_ind_loci, num_linked_loci[0]);

                  num_stat_displayed=2;
                  list_name_stat[1]="Rho      ";
                  list_name_stat[2]="Rho'     ";
                  write_output_RecomStats(ofcoaltime, num_stat_displayed, list_name_stat,
                                         list_stat_time, num_ind_loci, num_linked_loci[0]);

                   //if(num_linked_loci == 2) {
                   for(int j=0; j<(num_linked_loci[0]-1); ++j) {
                     int k=(int) Num_recomb_hits[num_ind_loci*(num_linked_loci[0]-1) + j];
                     if(k>=upBoundary) {
                        ++recomb_distrid[upBoundary-1][num_ind_loci*(num_linked_loci[0]-1) + j];
                     }
                     else {
                        ++recomb_distrid[k][num_ind_loci*(num_linked_loci[0]-1) + j];
                     }
                     k=(int) Num_effective_recomb_hits[num_ind_loci*(num_linked_loci[0]-1) + j];
                     if(k>=upBoundary) {
                        ++eff_recomb_distrid[upBoundary-1][num_ind_loci*(num_linked_loci[0]-1) + j];
                     }
                     else {
                        ++eff_recomb_distrid[k][num_ind_loci*(num_linked_loci[0]-1) + j];
                     }
                   }
                }
                #endif //_COMPUTE_MOMENTS_

      		}
      		myDemes.reset(myDemeSizes,myGrowthRates);
         } //end of the loop for number of independent loci to simulate


         if (!fmod(1.0*(tot_cases+1),10.0)) {
            #ifdef _VERBOSE_
            cout << "\nGenealogy # " << (tot_cases+1) << "/" << num_rand_samples;
            #endif
         }

         #ifdef _PRINT_ARLEQUIN_OUTPUT_
                ArlInFile << "\n";
                myDemes.write_loci_to_Arlequin_file(ArlInFile, data_type, genotypic_data);
                myDemes.write_group_section_to_Arlequin_file(ArlInFile);
                ArlInFile.close();
         #endif //_PRINT_ARLEQUIN_OUTPUT_

         #ifdef _PRINT_PHASE_OUPUT_
      	for (int i=0, j=0; i<num_pop; ++i) if (myDemes[i].sample_size()) {
            //Build PHASE file name
            int deca, signa, ndiga = 0;
            my_string phase_name = 	paup_generic_name + "_" +
            								fcvt(tot_cases, ndiga, &deca, &signa) + "_";
            phase_name+= (my_string) (fcvt(i, ndiga, &deca, &signa)) + ".pha";
            ++j;

        		ofstream	PhaseFile(phase_name.c_str());

            //write_PHASE_Header(PhaseFile, num_linked_loci, mySampleSizes[i]/2, data_type);
            //myDemes[i].print_loci_for_Phase(PhaseFile, i, data_type);

            //Recombination: mix of different loci
            myDemes[i].print_loci_for_Phase(PhaseFile, i, num_linked_loci, &data_type_perLoc);

            PhaseFile.close();
         }
         #endif //_PRINT_PHASE_OUPUT_

         #ifdef _PRINT_HAPLOTYPER_OUTPUT_
      	if ((data_type==DNA  || data_type==SNP) && genotypic_data) //Output only for DNA sequences and genotypic data
      	for (int i=0, j=0; i<num_pop; ++i) if (myDemes[i].sample_size()) {
            //Build HAPLOTYPER file name
  		  		int deca, signa, ndiga = 0;
            my_string haplotyper_name = 	paup_generic_name + "_" +
            										fcvt(tot_cases, ndiga, &deca, &signa) + "_";
            haplotyper_name+= (my_string) (fcvt(i, ndiga, &deca, &signa)) + ".hap";
            ++j;

        		ofstream	HaplotyperFile(haplotyper_name.c_str());

            myDemes[i].print_loci_for_Haplotyper(HaplotyperFile, i);

            HaplotyperFile.close();
         }
         #endif //_PRINT_HAPLOTYPER_OUTPUT_

         myDemes.flushLoci();  //Get rid of previously simulated loci
         #ifdef _COMPUTE_MOMENTS_
         if (num_chrom_structs==1) {
            myDemes.resetRecCounts();
         }
         #endif
   	}//End of loop for simulations

      #ifdef _COMPUTE_MOMENTS_
      if (num_chrom_structs==1) {
         double* RHudson= new double[num_pop];
         double meanRHudson=0;
         for (int i=0; i<num_pop; ++i) {
            RHudson[i]=2*myDemeSizes[i]*rec_rate[0];
            meanRHudson+=RHudson[i];
         }
         meanRHudson=meanRHudson/num_pop;

         for (int num_ind_loci=0; num_ind_loci<num_indep_loci; ++num_ind_loci ) {
            ofstream ofcoaltime(CoalTimes_out[num_ind_loci].c_str(),ios::app);
            double meanDemeSize=0;
            for (int i=0; i<num_pop; ++i) {
               ofcoaltime << "\nDeme " << (i+1) << ";";
               ofcoaltime << "\n\tR=4*"<<(myDemeSizes[i]/2.)<<"*"<<rec_rate[0]<<"=" << RHudson[i] << "; ";
               ofcoaltime << "\n\tTheta=4*"<<(myDemeSizes[i]/2.)<<"*"<<mut_rate[0]<<"=" << 4*(myDemeSizes[i]/2.)*mut_rate[0] << "; ";
               meanDemeSize+=myDemeSizes[i];
            }
            meanDemeSize=meanDemeSize/num_pop;
            if (num_pop>02){
               ofcoaltime << "\nmean betweeen demes;";
               ofcoaltime << "\n\tR=4*"<<(meanDemeSize/2.)<<"*"<<rec_rate[0]<<"=" << meanRHudson << "; ";
               ofcoaltime << "\n\tTheta=4*"<<(meanDemeSize/2.)<<"*"<<mut_rate[0]<<"=" << 4*(meanDemeSize/2.)*mut_rate[0] << "; ";
            }

            if(num_rand_samples>1) {
               for(int curLocus=0;curLocus<num_linked_loci[0];++curLocus) {
                  int k=num_ind_loci*num_linked_loci[0];
                  mean_T_MRCA[k + curLocus]=mean_T_MRCA[k + curLocus]/num_rand_samples;
                  //mean_within_loci_mean_time[curLocus]=mean_within_loci_mean_time[curLocus]/num_rand_samples;
                  //mean_within_loci_sd_time[curLocus]=mean_within_loci_sd_time[curLocus]/num_rand_samples;
                  //mean_within_loci_sum_time[curLocus]=mean_within_loci_sum_time[curLocus]/num_rand_samples;
                  //covTime_first_others[curLocus]=covTime_first_others[curLocus]/num_rand_samples;
                  //covTime_first_others[curLocus]=covTime_first_others[curLocus] -
                  //(mean_within_loci_sum_time[0]*mean_within_loci_sum_time[curLocus]);
                  //if(sqrt( (covTime_first_others[curLocus])*(covTime_first_others[curLocus]) )<0.000000001) {
                  //   covTime_first_others[curLocus]=0;
                  //}//Between_loci_mean_T[curLocus]/=num_rand_samples;
                  //Between_loci_sd_T[curLocus]/=num_rand_samples;
                  //THudson[curLocus]/=num_rand_samples;
               }
               for(int curLocus=0;curLocus<(num_linked_loci[0]-1);++curLocus) {
                  int k=num_ind_loci*(num_linked_loci[0]-1);
                  mean_Num_recomb_hits[k + curLocus]=mean_Num_recomb_hits[k + curLocus]/num_rand_samples;
                  mean_Num_effective_recomb_hits[k + curLocus]=mean_Num_effective_recomb_hits[k + curLocus]/num_rand_samples;
               }
               ofcoaltime << "\n\nBetween " << num_rand_samples << " simulations;";
               ofcoaltime << "\nT_MRCA   ;";
               for (int cpt=0;cpt<num_linked_loci[0];cpt++) {
               ofcoaltime << fixed << setw(11) << setprecision(1) << showpoint << mean_T_MRCA[num_ind_loci*num_linked_loci[0] + cpt]       << ";" ;
               }
               ofcoaltime << "\nT_MRCA/2N;";
               for (int cpt=0;cpt<num_linked_loci[0];cpt++) {
               ofcoaltime << fixed << setw(11) << setprecision(3) << showpoint << mean_T_MRCA[num_ind_loci*num_linked_loci[0] + cpt]/meanDemeSize<< ";" ;
               }

               ofcoaltime << "\nr        ;       ";
               for (int cpt=0;cpt<(num_linked_loci[0]-1);cpt++) {
                double temp=*rec_rate_perLoc[0].elem(cpt);
                  ofcoaltime << fixed <<  setw(11)  << setprecision(7) << showpoint
                           <<  temp << ";" ;
               }

               ofcoaltime << "\n4Nr      ;       ";
               for (int cpt=0;cpt<(num_linked_loci[0]-1);cpt++) {  
                double temp=*rec_rate_perLoc[0].elem(cpt);
                  ofcoaltime << fixed <<  setw(11)  << setprecision(2) << showpoint
                             << ( 2*meanDemeSize*temp ) << ";" ;
               }

               ofcoaltime << "\npred Rho ;       ";
               for (int cpt=0;cpt<(num_linked_loci[0]-1);cpt++) {    
                double temp=*rec_rate_perLoc[0].elem(cpt);
                  ofcoaltime << fixed << setw(11) << setprecision(2) << showpoint
                           << (   6*(2*meanDemeSize*temp)/(6 +
                                       2*meanDemeSize*temp)   ) << ";" ;
               }
               ofcoaltime << "      => 6*R / (6+R) " ;
               ofcoaltime << "\nRho      ;       " ;
               for (int cpt=0;cpt<(num_linked_loci[0]-1);cpt++) {
                  ofcoaltime << fixed << setw(11) << setprecision(2) << showpoint
                           << mean_Num_recomb_hits[num_ind_loci*(num_linked_loci[0]-1) + cpt]  << ";";
               }
               ofcoaltime << "\nRho'     ;       ";
               for (int cpt=0;cpt<(num_linked_loci[0]-1);cpt++) {
                  ofcoaltime << fixed << setw(11) << setprecision(2) << showpoint
                             << mean_Num_effective_recomb_hits[num_ind_loci*(num_linked_loci[0]-1) + cpt]  << ";";
               }

               int downBoundary_displayed=0;
               int UpBoundary_displayed=100;
               ofcoaltime << "\n\nRho and Rho' distribution for every recombination position (position i between locus i and locus i+1):\n";
               for(int j=0; j<(num_linked_loci[0]-1); ++j) {
                  ofcoaltime << "\n--Position " << (j+1) << "--\n";
                  for(int i=downBoundary_displayed; i<UpBoundary_displayed; ++i) {
                     ofcoaltime << fixed << setw(5) << i << ";";
                  }
                  ofcoaltime << "\n" ;
                  for(int i=downBoundary_displayed; i<UpBoundary_displayed; ++i) {
                     ofcoaltime << fixed << setw(5)
                              << recomb_distrid[i][num_ind_loci*(num_linked_loci[0]-1) + j] << ";";
                  }
                  ofcoaltime << "\n" ;
                  for(int i=downBoundary_displayed; i<UpBoundary_displayed; ++i) {
                     ofcoaltime << fixed << setw(5)
                              << eff_recomb_distrid[i][num_ind_loci*(num_linked_loci[0]-1) + j] << ";";
                  }
               }
            }
            ofcoaltime.close();
         }//end write in every chromosome file

         if(RHudson)  delete[] RHudson;
      }
      #endif //_COMPUTE_MOMENTS_

      if (num_linked_loci) {
         tree_exp_length/=num_linked_loci[0];
         tree_exp_length2/=num_linked_loci[0]*num_linked_loci[0];
      }
      else {
         tree_exp_length=0;
         tree_exp_length2=0;
      }

      #ifdef _PRINT_SEQ_PAUP_
      if (num_chrom_structs==1) {
		   write_PAUP_footer(PaupFile);
      }
  	   #endif


      #ifdef _PRINT_TREE_PAUP_
      if (num_chrom_structs==1) {
         double tempsqrt=0.0;
         write_PAUP_end(PaupTrueTreesFile);

         if(data_type != SNP){
            write_PAUP_end(PaupMutTreesFile);
         }

         if (tot_num_rand_samples>1) {
            tempsqrt=((tree_exp_length2 - tree_exp_length*tree_exp_length/tot_num_rand_samples)/
                  			 (tot_num_rand_samples-1));
            if(fabs(tempsqrt)<= 0.000000000000001) {   // In some case tempsqrt can be equal to 0 but
               tempsqrt=0.0;                             // the difference between 2 floating numbers
            }
         }                                            // is not exactly equal to 0 (less than 10-17)
   		PaupTrueTreesFile
         	<< "\nMean length of the trees : "
            << (tree_exp_length/tot_num_rand_samples);
         if (tot_num_rand_samples>1)
         	PaupTrueTreesFile
               << " +- "
            	<< setiosflags(ios::showpoint | ios::fixed | ios::right)
            	<< setw(10)
            	<< setprecision(5)
               //<< ( sqrt((tree_exp_length2 - tree_exp_length*tree_exp_length/tot_num_rand_samples)/
               //			 (tot_num_rand_samples-1)) )
               << sqrt(tempsqrt)
               << endl;

         if(data_type != SNP){
            PaupMutTreesFile
         	   << "\nMean length of the trees : "
               << (tree_mut_length/tot_num_rand_samples);
            if (tot_num_rand_samples>1)
         	   PaupMutTreesFile
                  << " +- "
            	   << setiosflags(ios::showpoint | ios::fixed | ios::right)
            	   << setw(10)
            	   << setprecision(5)
                  << ( sqrt((tree_mut_length2 - tree_mut_length*tree_mut_length/tot_num_rand_samples)/
               			   (tot_num_rand_samples-1)) )
                  << endl;
         }
      }
      #endif //_PRINT_TREE_PAUP_

   	start=clock()-start;
   	//cout << "\n\ndone in " <<  start << "ms\n" << endl;
        double donein=start/ (double) CLOCKS_PER_SEC;
      #ifdef _VERBOSE_
      cout << "\n\ndone in " <<  start/ (double) CLOCKS_PER_SEC << "sec\n" << endl;
      cout << "\n\ndone in " <<  donein << "sec\n" << endl;
      #endif

      #ifdef _PRINT_ARLEQUIN_OUTPUT_
      ArlSimulConditions << endl ;
      ArlSimulConditions << "\n\ndone in " <<  donein << "sec\n" << endl;
      #endif

      #ifdef _COMPUTE_MOMENTS_
      if (num_chrom_structs==1) {
   	   MeanCoalMat.compute_mean(tot_num_rand_samples);
   		sdCoalMat.compute_sd(tot_num_rand_samples,MeanCoalMat);

   	   MeanPairDiff.compute_mean(tot_num_rand_samples);
   		sdPairDiff.compute_sd(tot_num_rand_samples,MeanPairDiff);
      }
      #endif

      //Indices for computing pairwise Fst's (Weir's theta)
      double 	theta, t1, t0, t0_1, t0_2, time;

      #ifdef _VERBOSE_
   	cout << "\n";
      #endif

      //ofs	<< "\nFor every region, total mutation rate per generation          : "
      //   	<< mut_rate;
      //ofs	<< "\nMutation rate per locus per generation                        : "
      //   	<< mut_rate/num_linked_loci;
      //ofs   << "\n\n";


      #ifdef _COMPUTE_MOMENTS_

   		//ofs 	<< "\nMean coalescence times within and among demes over " << tot_cases <<" random samples\n"
   		//		<< MeanCoalMat;
   		//ofs 	<< "\nS.D. coalescence times within and among demes \n"
   		//		<< sdCoalMat;
   		//ofs 	<< "\nMean number of pairwise differences within and among demes over " << tot_cases <<" random samples\n"
   		//		<< MeanPairDiff;
   		//ofs 	<< "\nS.D. number of pairwise differences within and among demes \n"
   		//		<< sdPairDiff;

         if (num_chrom_structs==1) {
                for(int i=0;i<num_indep_loci;++i) {
                        CoalTimesFile[i].close();
                }
         }
         if(CoalTimesFile)                delete []   CoalTimesFile;
         if(CoalTimes_out)                delete []   CoalTimes_out;

         if(list_name_stat)               delete []   list_name_stat;
         if(list_stat_time)               delete []   list_stat_time;
         if(T_MRCA)                       delete []   T_MRCA;
         //if(within_loci_mean_time)        delete []   within_loci_mean_time;
         //if(within_loci_sd_time)          delete []   within_loci_sd_time;
         if(Num_recomb_hits)              delete []   Num_recomb_hits;
         if(Num_effective_recomb_hits)    delete []   Num_effective_recomb_hits;
         //if(ProdSumTime)                  delete []   ProdSumTime;
         //if(sumTime)                      delete []   sumTime;
         //if (THudson)                 delete []   THudson;
         //if (Between_loci_mean_T)     delete []   Between_loci_mean_T;
         //if (Between_loci_sd_T)       delete []   Between_loci_sd_T;


         if(mean_T_MRCA)                  delete []   mean_T_MRCA;
         //if(mean_within_loci_mean_time)   delete []   mean_within_loci_mean_time;
         //if(mean_within_loci_sd_time)     delete []   mean_within_loci_sd_time;
         //if(mean_within_loci_sum_time)    delete []   mean_within_loci_sum_time;
         //if(covTime_first_others)         delete []   covTime_first_others;
         if(mean_Num_recomb_hits)         delete []   mean_Num_recomb_hits;
         if(mean_Num_effective_recomb_hits) delete []   mean_Num_effective_recomb_hits;

         for(int i=0;i<upBoundary;++i) {
            if(recomb_distrid[i])      delete [] recomb_distrid[i];
            if(eff_recomb_distrid[i])  delete [] eff_recomb_distrid[i];
         }
         if(recomb_distrid)               delete []   recomb_distrid;
         if(eff_recomb_distrid)           delete []   eff_recomb_distrid;

      #endif //_COMPUTE_MOMENTS_

      //Do not forget to reset the static member of TNode
      TNode resetNode;
      resetNode.reset_node_count();

	  //removed so that directory doesnt change mjj 5/15/06
	  /*
      #ifdef  _LINUX_
      chdir("..");
      #else
      chdir("..");
      #endif
	   */
      if (mut_ratio)                delete[] mut_ratio;
      if (prop_sites)               delete[] prop_sites;
      if (num_rate_categories)      delete[] num_rate_categories;
      if (num_linked_loci)          delete[] num_linked_loci;
      if (rec_rate)                 delete[] rec_rate ;
  		if (mut_rate)                 delete[] mut_rate;
      if (num_linked_loci_perBlock) delete[] num_linked_loci_perBlock;
      if (transRatePerLoc)          delete[] transRatePerLoc;
      if (rec_rate_perLoc)          delete[] rec_rate_perLoc;
      if (mut_rate_perLoc)          delete[] mut_rate_perLoc;
      if (geom_param_perLoc)        delete[] geom_param_perLoc;
      if (rangeConstraintPerLoc)    delete[] rangeConstraintPerLoc;
      if (freq_SNP_min) {
         for (int p=0; p<num_chrom_structs; ++p) {
            if (freq_SNP_min[p]) delete []  freq_SNP_min[p];
         }
         delete [] freq_SNP_min;
      }
      if (data_type_perLoc)         delete[] data_type_perLoc;
      if (Data_mix)                 delete[] Data_mix ;
      if (num_linkage_blocks)       delete[] num_linkage_blocks;
	}
   return 1;
};









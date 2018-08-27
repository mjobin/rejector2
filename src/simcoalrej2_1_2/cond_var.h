#ifndef _COND_VAR_
   //#define _COND_VAR_
   #define _LINUX_              		 //Compile on an Unix plateform
                                        //By default Windows
   
   #define _PRINT_ARLEQUIN_OUTPUT_      //Output simulazed data in Arlequin format
   //#define _PRINT_SNP_ONLY_             //Output only polymorphic positions (for DNA data type only)

   
   //#define _PRINT_TREE_PAUP_            //Output generated trees in PAUP format (one tree per locus)
   
   //#define _COMPUTE_MOMENTS_            //Compute moments of coalescence times
                                        // compute distribution of coalescence time 
                                        // for ervery partially linked loci
                                        // store them in a file name "CoalTimesFile.txt"

                                       //by default one coalescent event   
   //#define _MCE_                     //Multiple coalescent event
   #define _MIXT_                      //Automatic check if one or multiple coalescent event
                                       //(MIXT and MRCE cannot be simultaneously uncommented)
                                       
                                       //by default one Recombination event
   #define _MRE_                       //Multiple Recombination event per generation

   //#define _VERBOSE_                                          
                                       
#endif

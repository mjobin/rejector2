#!/usr/bin/perl
# use strict

# A script for making rejector2 input files.

# Written by Matthew Jobin, Department of Anthropological Sciences, Stanford University

$numArgs = $#ARGV + 1;


if($numArgs != 2) {
	print "usage: perl makerejin.pl <priors and tests> <outfile>\n";
	exit;
	}

$priorsandtests = $ARGV[0];
$outfile = $ARGV[1];

open(INF, $priorsandtests) || die "couldn't open " . $priorsandtests;

# and output
open(OUTF, ">$outfile");
print "Writing to " . $outfile . "\n";

print OUTF "REJECTOR\n\n";




#
my $linked = 0;


my $recrate = 0.0;
my $avphysdist = 100000;

print "All same?(Y/N)\n";
$a = <STDIN>;
chop $a;
do
{
	if($a eq "Y") 
	{
		$linked = 1;
		

		
	}
	elsif($a eq "N")
	{
		$linked = 0;
	}
	else
	{
		print "Answer Y or N\n";
		$a = <STDIN>;
		chop $a;
	}
}
while($a ne "Y" && $a ne "N");


		print "Recombination rate?\n";
		$a = <STDIN>;
		chop $a;
		$recrate = $a;
		
		
		
		print "Average Phyisical Distance?\n";
		$a = <STDIN>;
		chop $a;
		$avphysdist = $a;


my $locitype = 0;
print "Type 1 for SNPs, 2 for STRs, 3 for UEPs, 4 for DNA sequence.\n";
$a = <STDIN>;
chop $a;
$locitype = $a;

do
{
	$locitype = $a;
	if($locitype < 1 || $locitype > 7)
	{
		print "Choose 1 to 5.\n";
		$a = <STDIN>;
		chop $a;
	}
}
while ($locitype < 1 && $locitype > 7);


my $loci = 0;
print "How many loci?\n";
$a = <STDIN>;
chop $a;
$loci = $a;
do
{
	$loci = $a;
	if($locitype < 1)
	{
		print "Positive integer please.\n";
		$a = <STDIN>;
		chop $a;
	}
}
while ($loci < 1 );


my $snprate = 1.0e-9;
my $strrate = 2.6e-4;
my $ueprate = 1.0e-9;
my $dnarate = 3.0e-08;
my $dnasites = 1000;

my $snpanc = "A";
my $stranc = "NA";
my $uepanc = "S";
my $dnaanc = "NA";


if ($locitype == 1)
{
	print "SNP mutation rate?\n";
	$a = <STDIN>;
	chop $a;
	$snprate = $a;
}

elsif ($locitype == 2)
{
	print "STR mutation rate?\n";
	$a = <STDIN>;
	chop $a;
	$strrate = $a;
}

elsif ($locitype == 3)
{
	print "UEP mutation rate? \n";
	$a = <STDIN>;
	chop $a;
	$ueprate = $a;
}

elsif ($locitype == 4)
{
	print "num sites per locus?\n";
	$a = <STDIN>;
	chop $a;
	$dnasites = $a;
	

	
}


my $pops = 0;
print "How many populations?\n";
$a = <STDIN>;
chop $a;
do
{
	$pops = $a;
	if($pops < 1)
	{
		print "Positive integer please.\n";
		$a = <STDIN>;
		chop $a;
	}
}
while ($pops < 1);

#


my @samplist;
for $i (1 .. $pops)
{
	print "How many samples in pop $i?\n";
	$a = <STDIN>;
	chop $a;
	push(@samplist, $a);
	
	

}



print OUTF "Number of Populations\n";
print OUTF $pops . "\n";




my $infline;
foreach $infline (<INF>){
	chomp $infline;
	print OUTF $infline . "\n";
}
close (INF);

print OUTF "\n";


print OUTF "/--Data\nVnaught\t0\n\n";

my $i;

if($linked)
{
print "linked\n";
	print OUTF "Loci\t";

		if ($locitype == 1)
		{
			print OUTF "SNP\t";
		}
		elsif ($locitype == 2)
		{
			print OUTF "STR\t";
		}
		elsif ($locitype == 3) 
		{
			print OUTF "UEP\t";
		}

		elsif ($locitype == 4) 
		{
			print OUTF "DNA\t";
		}
	print OUTF "\n";
	
	print OUTF "MutRates\t";

		if ($locitype == 1)
			{
				print OUTF $snprate . "\t";
			}
			elsif ($locitype == 2)
			{
				print OUTF $strrate . "\t";
			}
			elsif ($locitype == 3) 
			{
				print OUTF $ueprate . "\t";
			}
			
			
			elsif ($locitype == 4) 
			{
				print OUTF $dnarate . "\t";
			}
	
	print OUTF "\n";
	
	print OUTF "Ancestral\t";

		if ($locitype == 1)
		{
				print OUTF $snpanc . "\t";
			}
			elsif ($locitype == 2)
			{
				print OUTF $stranc . "\t";
			}
			elsif ($locitype == 3) 
			{
				print OUTF $uepanc . "\t";
			}

			elsif ($locitype == 4) 
			{
				print OUTF $dnaanc . "\t";
			}
	


	print OUTF "\n";
	


print OUTF "RecombRt\t" . $recrate . "\n";

print OUTF "PhysDist\t" . $avphysdist . "\n";






print OUTF "NumLoci\t" . $loci . "\n";

print OUTF "Length\t";


		if ($locitype == 1)
		{
				print OUTF "1\t";
			}
			elsif ($locitype == 2)
			{
				print OUTF "2\t";
			}
			elsif ($locitype == 3) 
			{
				print OUTF "1\t";
			}

			elsif ($locitype == 4) 
			{
				print OUTF $dnasites. "\t";
			}

}

else
{

	
	print OUTF "Loci\t";
	for $i (1 .. $loci)
	{
		
		if ($locitype == 1)
		{
			print OUTF "SNP\t";
		}
		elsif ($locitype == 2)
		{
			print OUTF "STR\t";
		}
		elsif ($locitype == 3) 
		{
			print OUTF "UEP\t";
		}
		elsif ($locitype == 4) 
		{
			print OUTF "DNA\t";
		}

	}
	print OUTF "\n";
		
		
	print OUTF "MutRates\t";
	for $i (1 .. $loci){
		
		if ($locitype == 1)
		{
			print OUTF $snprate . "\t";
		}
		elsif ($locitype == 2)
		{
			print OUTF $strrate . "\t";
		}
		elsif ($locitype == 3) 
		{
			print OUTF $ueprate . "\t";
		}
		
		
		elsif ($locitype == 4) 
		{
			print OUTF $dnarate . "\t";
		}
		
	}
	print OUTF "\n";
		
	
	print OUTF "Ancestral\t";
	for $i (1 .. $loci){
		if ($locitype == 1)
		{
			print OUTF $snpanc . "\t";
		}
		elsif ($locitype == 2)
		{
			print OUTF $stranc . "\t";
		}
		elsif ($locitype == 3) 
		{
			print OUTF $uepanc . "\t";
		}
		
		elsif ($locitype == 4) 
		{
			print OUTF $dnaanc . "\t";
		}
		
		
	}
	print OUTF "\n";
	
	print OUTF "RecombRt\t";
	for $i (1 .. $loci){
		print OUTF $recrate . "\t";
	}
	print OUTF "\n";
	
	
	print OUTF "PhysDist\t";
	for $i (1 .. $loci){
		print OUTF (rand() * $avphysdist) . "\t";
	}
	print OUTF "\n";
	
	print OUTF "NumLoci\t";
	for $i (1 .. $loci){
		print OUTF "1\t";
	}
	print OUTF "\n";
	
	print OUTF "Length\t";
	for $i (1 .. $loci){
		
		if ($locitype == 1)
		{
			print OUTF "1\t";
		}
		elsif ($locitype == 2)
		{
			print OUTF "2\t";
		}
		elsif ($locitype == 3) 
		{
			print OUTF "1\t";
		}
		
		elsif ($locitype == 4) 
		{
			print OUTF $dnasites. "\t";
		}
	}
	print OUTF "\n";
		
		
	print OUTF "\n\n";
}


print OUTF "\n\n";
#
print OUTF "Tag\tPopulation\n";

if($linked)
{

	for $i (1 .. $pops)
	
{
	
	my $j;
	
	for $j (1 .. $samplist[($i-1)])
	{
		print OUTF ($i-1) . "_" . $j. "\t" . $i . "\t";
		my $k;
		
		for $k (1 .. $loci)
		{
			
			
			
		if ($locitype == 1)
		{
			if(rand() > 0.5)
			{
				print OUTF "A";
			}
			else
			{
				print OUTF "C";
			}
			
			if($k < $loci)
			{
				print OUTF "\t";
			}
			
		}
		elsif ($locitype == 2)
		{
			print OUTF int(rand(15)) + 1;
			
			if($k < $loci)
			{
				print OUTF "\t";
			}
		}
		elsif ($locitype == 3) 
		{
			if(rand() > 0.5)
			{
				print OUTF "S";
			}
			else
			{
				print OUTF "L";
			}
		if($k < $loci)
			{
				print OUTF "\t";
			}
		}

			
		elsif ($locitype == 4) 
		{
			for $ds (1 .. $dnasites)
			{
				if(rand() > 0.25)
				{
					print OUTF "A";
				}
				elsif(rand() > 0.25 && rand() < 0.5)
				{
					print OUTF "C";
				}
				elsif(rand() > 0.5 && rand() < 0.75)
				{
					print OUTF "G";
				}
				else
				{
					print OUTF "T";
				}
				
			if($k < $loci)
			{
				print OUTF "\t";
			}
			
			}
		
		
		
		}

			
	
		}
	print OUTF "\n";
	
	}
	
	}
	
	
}



else
{
	
	for $i (1 .. $pops)
		
	{
		
		my $j;
		
		for $j (1 .. $samplist[($i-1)])
		{
			print OUTF ($i-1) . "_" . $j. "\t" . $i . "\t";
			my $k;
			
			for $k (1 .. $loci)
			{
				
				
				
				if ($locitype == 1)
				{
					if(rand() > 0.5)
					{
						print OUTF "A";
					}
					else
					{
						print OUTF "C";
					}
					
					if($k < $loci)
					{
						print OUTF "\t";
					}
					
				}
				elsif ($locitype == 2)
				{
					print OUTF int(rand(15)) + 1;
					
					if($k < $loci)
					{
						print OUTF "\t";
					}
				}
				elsif ($locitype == 3) 
				{
					if(rand() > 0.5)
					{
						print OUTF "S";
					}
					else
					{
						print OUTF "L";
					}
					if($k < $loci)
					{
						print OUTF "\t";
					}
				}
				
				
				elsif ($locitype == 4) 
				{
					for $ds (1 .. $dnasites)
					{
						if(rand() > 0.25)
						{
							print OUTF "A";
						}
						elsif(rand() > 0.25 && rand() < 0.5)
						{
							print OUTF "C";
						}
						elsif(rand() > 0.5 && rand() < 0.75)
						{
							print OUTF "G";
						}
						else
						{
							print OUTF "T";
						}
						
						if($k < $loci)
						{
							print OUTF "\t";
						}
						
					}
					
					
					
				}
				
				
				
			}
			print OUTF "\n";
			
		}
		
	}
	
	
}






close(OUTF);
		
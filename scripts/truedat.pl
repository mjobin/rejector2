#!/usr/bin/perl
# use strict;
# Processes all .txt outfiles in folder. goes from <file>.rej0.txt up to the first missing number in a series
# , so if you have <file>.rej0.txt, <file>.rej1.txt and <file>.rej3.txt it won't read the last one.
# 
# usage - perl gluedat.pl <basefile>, where basefil is the base of the common filename 
# .rej0.txt. In ither words, if you have x.rej0.txt and x.rej1.txt, type 
# perl gluedat.pl x

# Written by Matthew Jobin, Department of Anthropological Sciences, Stanford University

$numArgs = $#ARGV + 1;


if($numArgs != 1) {
	print "usage:  perl gluedat.pl <basefile>, where basefil is the base of the common filename .rej0.txt. In other words, if you have x.rej0.txt and x.rej1.txt, type perl gluedat.pl x";
	exit;
	}
	
$filebase = $ARGV[0];



$firstfile = $filebase . ".rej0.txt";



if (!-e $firstfile) {
	print "can't find first file: " . $firstfile . "\n";
	exit;
	}
	
$outfile = $filebase . "-out.txt";
open(OUTF, ">$outfile");
print "Writing to " . $outfile . "\n";
	
$i = 0;
my $headersize = -1;
my @data;
my $lncount = 0;
$filesin = 1;

while($filesin){
$j = 0;
	$infile = $filebase . ".rej" . $i . ".txt";

	if (-e $infile) {
		# convert accursed mac line breaks!
		$cmd = "perl -pi -e 's/\r/\n/g' ". $infile;
		#print $cmd."\n";
		if(system($cmd)) { print "convert failed\n"; }
		
		
		open(INF, $infile) || die "couldn't open " . $infile;
		print "Reading " . $infile . "\n";


		foreach $infline (<INF>){
			chomp($infline);
			
			if($infline =~ m/^Rejector Output/)
			{
				#do nothing
			}
			
			elsif($infline =~ m/^Rejector2 Output/)
			{
				#do nothing
			}
			
			elsif($infline =~ m/^Infile/)
			{
				#do nothing
			}
			
			elsif($infline =~ m/^Run Limit/)
			{
				#do nothing
			}
			
			elsif($infline =~ m/^Printing started/)
			{
				#do nothing
			}
			elsif($infline =~ m/^Deme/)
			{

				if($i == 0){ # only deal with header the first time
					@headerline = split(/\t/, $infline);
					$headersize = @headerline;
					print "size of header: " . $headersize . "\n";
					print OUTF "$infline\n";
				}
			}
			elsif($infline =~ m/^\d/ || $infline =~ m/^-/) #Match first part of line to either a digit or - sign
			{
				my $trueline = 1;
				
				@dataline = split(/\t/, $infline);

				if(@dataline != $headersize) {
					print "Error: File " . $infile . " line " . $j . ". Line has " . @dataline . " elements, but the first file's header has " . $headersize . ". Exiting.\n";
					exit;
				}
				
				
				foreach (@dataline) {
					if($_ eq "F"){
						$trueline = 0;
						}

					}
				
				if($trueline){
				print OUTF "$infline\n";
				}
				$lncount++;

				
			}
			elsif($infline =~ m /^$/)
			{
				#empty line, do nothing
			
			}
			else
			#elsif($infline =~ m/^\D/)
			{
				print "Error. Found line starting with non-digit that is not in the header. Skipping." . "\n";
				print $infline . "\n";
			}


		
			$j++;
		}
		
		
		close(INF);
		
	
		
	
	}
	else{
		#print "No match. Exiting read-in.\n";
		$filesin = 0;
	}
	$i++;
}

print $lncount . " lines read.\n";

close(OUTF);

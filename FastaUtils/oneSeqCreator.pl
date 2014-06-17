#!/usr/bin/perl
#Create a one sequence fasta file from a multi-fasta file.
#Usage : oneSeqCreator.pl multi-fasta-file one-fasta-file
# Clement DELESTRE 5/05/2014
$inputfile=@ARGV[0];  #Multi fasta file
$outputFile=@ARGV[1];  #One fasta file
$count=0;
open(FILE,$inputfile);
open(OUT,">".$outputFile) or die ("Error with output file\n");
while(<FILE>){
	if ($_=~/^>(.*)/){ 
		if ($count!=0){
			break;
		}
		else {
			print OUT $_;
			$count++;
		}
	}
	else {
		print OUT $_;
	}
}
close(OUT);
close(FILE);
#!/usr/bin/perl
# Extract all sequences from a file of sequences (fasta format). One file by sequence will be created.
# Usage : extractAllSeq.pl file
# Clement DELESTRE 5/05/2014
$file=@ARGV[0];
open(FILE,$file) || die ("File error\n") ;
$write=0;
$first=1;
while(<FILE>){
	chomp();
	if ($_=~/^>(.*)$/){
		if ($first==0){
			close(RESULT);
		}
		open(RESULT,">>$1.fa") || die ("Result file error \n") ;
		$write=1;
		$first=0;
	}
	elsif ($write==1){
		print RESULT"$_\n";
	}
}
close(RESULT);
close(FILE);
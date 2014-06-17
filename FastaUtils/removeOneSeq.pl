#!/usr/bin/perl
# Remove one sequence from a file of sequences (FASTA format)
# Usage : removeOneSeq.pl query multi-fasta-file > out
# Clement DELESTRE 5/05/2014
$IN=@ARGV[0];  # sequence name to look for in the db
$DB=@ARGV[1];  # db name
$SEQ=">"."$IN";
$found=$iflag=0;
open(DB,"$DB") || die "can't open $DB\n";
while(<DB> ){
     chomp();
     if ( m />/ && $iflag == 1 ) {
		$iflag=0;
	}
     if ( $_ eq $SEQ ){ 
		$iflag=1; 
	}
     if ( $iflag == 0 ){
       print $_."\n";
   }
}
close(DB);

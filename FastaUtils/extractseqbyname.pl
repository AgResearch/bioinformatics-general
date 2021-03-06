#!/usr/bin/perl
# Extract one sequence from a file of sequences
# Usage : extractseqbyname.pl sequence database > outputfile
# ET 25/10/2001
# Update CD 5/05/2014
$IN=@ARGV[0];  # sequence name to look for in the db
$DB=@ARGV[1];  # db name
$SEQ=">"."$IN";
$found=$iflag=0;
open(DB,"$DB") || die "can't open $DB\n";
while(<DB> ){
     chomp();
     if ( m />/ && $found == 1 ) {
		exit 0;
	}
     if ( $_ eq $SEQ ){ 
		$iflag=1; 
	}
     if ( $iflag == 1 ){
       print $_."\n";
       $found=1;
   }
}
close(DB);

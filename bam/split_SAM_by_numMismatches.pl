#!/usr/bin/env perl

$counter = 0 ;
while( $line = <STDIN> ) {
	# SAM/BAM header
	if ( $line =~/^@/ ) {
		print $line;
		next;
	}
	# Reads
	if ($line =~/MD\:Z\:(.*?)\s/ ) {	# Only present for mapped reads
		$a = $1 ;						# Get mismatches description tag field
		$a =~s/\d+//g ;					# Lose the runs of matching positions.
		$a =~s/\^[A-Z]+//g ;			# Don't count deletions for the stratification.
		$num = length $a ;				# How many letters are left?
		if ($num >= $ARGV[0] && $num <= $ARGV[1]) {
			$counter++ ;
			print STDOUT $line ;
#			print STDERR "$num - @mis\n"
		}
	}
}
print STDERR "$counter entries found with $ARGV[0] - $ARGV[1] mismatches.\n" ;

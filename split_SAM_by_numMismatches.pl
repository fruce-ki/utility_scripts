#!/users/kimon.froussios/miniconda3/envs/mybasics/bin/perl

$counter = 0 ;
while( $line = <STDIN> ) {
	# SAM/BAM header
	if ( $line =~/^@/ ) {
		print $line;
		next;
	}
	# Reads
	$line =~/MD\:Z\:(.*?)\s/ ; 
	$a = $1 ; 
	@mis = split /\d+/, $a ; 
	$num = $#mis ;
	if ($num >= $ARGV[0] && $num <= $ARGV[1]) {
		$counter++ ;
		print STDOUT $line ;
#		print STDERR "$num - @mis\n"
	}
}
print STDERR "$counter entries found with $ARGV[0] - $ARGV[1] mismatches.\n" ;

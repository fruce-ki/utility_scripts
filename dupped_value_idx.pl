#!/usr/bin/env perl

$header = <STDIN>;
chomp $header;
@fields = split "\t", $header;
$i=0;
foreach $field (@fields) {
	if ($field eq $ARGV[0]) {
		print "$i " ;
	}
	$i++;
}
print "\n";

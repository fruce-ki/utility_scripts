#!/usr/bin/env sh

comparisondir=$1 ; shift
outputdir=$1 ; shift

if [[ ! -d ${outputdir} ]]; then
	 mkdir -p ${outputdir}
fi

treatment=$(perl -e '$ARGV[0]=~/\/([a-zA-Z0-9_\-.]+)_vs_(\S+)$/;print $1;' $comparisondir)
reference=$(perl -e '$ARGV[0]=~/\/([a-zA-Z0-9_\-.]+)_vs_(\S+)$/;print $2;' $comparisondir)

echo "Treat: ${treatment}   Ref: ${reference}"

cat ${comparisondir}/guides_stats.txt \
	| perl -e '$header = <STDIN>;
		   chomp $header ;
			 
		   @fields = split /\t/, $header;
		   foreach $name (@fields) {
		   	if ($name ne "id" && $name ne "group") {
			   	$header =~s/(?<![\w.])$name(?![\w.])/M.sg.$name.$ARGV[0].vs.$ARGV[1]/;
			}
		   }
		   #$header =~s/_/./g;
		   print STDOUT "$header\n";
		   while($line = <STDIN>) {
			print STDOUT $line;
		   }' $treatment $reference \
	> ${outputdir}/guides_stats.txt


cat ${comparisondir}/genes_pos_stats.txt \
	| perl -e '$header = <STDIN>;
		   chomp $header ;
		   @fields = split /\t/, $header;
		   foreach $name (@fields) {
		   	if ($name ne "group") {
			   	$header =~s/(?<![\w.])$name(?![\w.])/M.g+.$name.$ARGV[0].vs.$ARGV[1]/;
			}
		   }
		   #$header =~s/_/./g;
		   print STDOUT "$header\n";
		   while($line = <STDIN>) {
			print STDOUT $line;
		   }' $treatment $reference \
	> ${outputdir}/genes_pos_stats.txt


cat ${comparisondir}/genes_neg_stats.txt \
	| perl -e '$header = <STDIN>;
		   chomp $header ;
		   @fields = split /\t/, $header;
		   foreach $name (@fields) {
		   	if ($name ne "group") {
			   	$header =~s/(?<![\w.])$name(?![\w.])/M.g-.$name.$ARGV[0].vs.$ARGV[1]/;
			}
		   }
		   #$header =~s/_/./g;
		   print STDOUT "$header\n";
		   while($line = <STDIN>) {
			print STDOUT $line;
		   }' $treatment $reference \
	> ${outputdir}/genes_neg_stats.txt

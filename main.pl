#!usr/bin/perl
use warnings;
use strict;
use lib './';
use HydrophobicMoment;
use Shuffe;

for(my $w=0;$w<5;$w++){
	unless(open(FILE,"helixnr70curated.fas")){
		die ("cannot open the input file");
	}

	my $out="2021HM". 10**$w;
	unless(open(OUT,">$out.csv")){
		die ("cannot create the output file");
	}

	while(my $line=<FILE>){
		chomp $line;
		if($line!~/>/){
			for(my $i=0;$i<100;$i++){
				printf OUT "%.3f;",Shuffe::averageHM($line,10**$w);
			}
			print OUT "\n";
		}
	}
	close FILE;
	close OUT;
}

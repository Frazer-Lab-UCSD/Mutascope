#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);



open(OUT, ">$ARGV[1]");
open(IN, $ARGV[0]);
while(my $line = <IN>){
	chomp $line;
	next if($line =~ m/^#/);
	my ($chr, $pos, $id, $ref, $alt, $qual, $filt, $info, $format, $normal, $tumor) = split(/\t/, $line);
	next if($filt ne "PASS");

	my @infos = split(/;/, $info);
	my $vt;

	foreach (@infos){
		if($_ =~ m/VT=/){
			$vt = $_;
			$vt =~ s/VT=//;
		}
	}
	
	my $af;
	my $ss;
	my @formats = split(/:/, $format);
	for(my $i = 0; $i <= $#formats; $i++){
		if($formats[$i] =~ m/AF/){
			$af = $i;
		}
		if($formats[$i] =~ m/SS/){
			$ss = $i;
		}
	}	
	
	my @norms = split(/:/, $normal);
	my @tums = split(/:/, $tumor);
	
	my $tss = $tums[$ss];
	my $nfreq = $norms[$af];
	my $tfreq = $tums[$af];
	
	print OUT "$chr\t$pos\t$nfreq\t$tfreq\t$tss\t$vt\n";
	
}
close(OUT);
close(IN);



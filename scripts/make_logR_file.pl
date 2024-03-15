#! /usr/bin/perl -w

use warnings;
use strict;

sub log2($);

my $tot_normal = 0;
my $tot_tumor = 0;
my %norm_amps;
my %tum_amps;
my %genes;

open(FILE, $ARGV[0]);
while(my $line = <FILE>){
	chomp $line;
	my ($chr, $start, $stop, $gene, $cov) = split(/\t/, $line);
	$norm_amps{$chr}{$start}{$stop} = $cov;
	$genes{$chr}{$start}{$stop} = $gene;
	$tot_normal += $cov;
}
close(FILE);

open(FILE, $ARGV[1]);
while(my $line = <FILE>){
	chomp $line;
	my ($chr, $start, $stop, $gene, $cov) = split(/\t/, $line);
	$tum_amps{$chr}{$start}{$stop} = $cov;
	$tot_tumor += $cov;
}
close(FILE);

print "#CHR\tSTART\tSTOP\tAMPLICON\tNORMAL_RC\tTUMOR_RC\tLOGR\n";

foreach my $chr (sort keys %genes){
	foreach my $start (sort {$a <=> $b} keys %{$genes{$chr}}){
		foreach my $stop (sort {$a <=> $b} keys %{$genes{$chr}{$start}}){
			my $g = $genes{$chr}{$start}{$stop};
			my $tcov = $tum_amps{$chr}{$start}{$stop};
			my $ncov = $norm_amps{$chr}{$start}{$stop};
			
			$tcov = 1 if($tcov == 0);
			$ncov = 1 if($ncov == 0);
			
			my $logr = log2(($tcov/$tot_tumor)/($ncov/$tot_normal));
			print "$chr\t$start\t$stop\t$g\t$ncov\t$tcov\t$logr\n";
		}
	}
}


sub log2($){
	my ($n) = @_;
	return (log($n)/log(2));
}

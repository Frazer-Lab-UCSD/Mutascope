#! /usr/bin/perl -w

use warnings;
use strict;

my %beds;
my %genes;

open(BED, $ARGV[0]);
while(my $line = <BED>){
	chomp $line;
	my ($chr, $start, $stop, $gene, @other) = split(/\t/, $line);
	$beds{$chr}{$start}{$stop} = 0;
	$genes{$chr}{$start}{$stop} = $gene;
}

open(FILE, $ARGV[1]);
while(my $line = <FILE>){
	chomp $line;
	next if($line =~ m/^@/);
	my (@fields) = split(/\t/, $line);
	my $rg_string;
	foreach (@fields){
		if($_ =~ m/^RG:Z:/){
			$rg_string = $_;
		}
	}
	
	my (@rgs) = split(/_/, $rg_string);
	my $location = $rgs[$#rgs];
	my ($chr, $pos) = split(/:/, $location);
	my ($start, $stop) = split(/-/, $pos);
	
	$beds{$chr}{$start}{$stop}++;
}

foreach my $c (sort keys %beds){
	foreach my $st (sort {$a <=> $b} keys %{$beds{$c}}){
		foreach my $se (sort {$a <=> $b} keys %{$beds{$c}{$st}}){
			print "$c\t$st\t$se\t$genes{$c}{$st}{$se}\t$beds{$c}{$st}{$se}\n";
		}
	}
}



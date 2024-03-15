#! /usr/bin/perl -w

use warnings;
use strict;

my @chroms = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"); 

my %covs;
my %refs;
my %rgs;
my %rp;
my %indls;


foreach my $file (@ARGV){
	open(FILE, $file);
	while(my $line = <FILE>){
		chomp $line;
		my ($chr, $pos, $ref, $cov, $rg_info, $indel_info) = split(/\t/, $line);
		if(exists($covs{$chr}{$pos})){
			$covs{$chr}{$pos} += $cov;
			push(@{$rgs{$chr}{$pos}}, $rg_info);
			push(@{$indls{$chr}{$pos}}, $indel_info) if($indel_info);
		}else{
			$covs{$chr}{$pos} = $cov;
			$refs{$chr}{$pos} = $ref;
			@{$rgs{$chr}{$pos}} = ();
			push(@{$rgs{$chr}{$pos}}, $rg_info);
			@{$indls{$chr}{$pos}} = ();
			push(@{$indls{$chr}{$pos}}, $indel_info) if($indel_info);
		}
	}
	close(FILE);
}

foreach my $chr (@chroms){
	if(exists($refs{$chr})){
		foreach my $pos (sort {$a <=> $b} keys %{$refs{$chr}}){
			print "$chr\t$pos\t$refs{$chr}{$pos}\t$covs{$chr}{$pos}\t";
			
			if(exists($indls{$chr}{$pos})){
				my $indf = 0;
				foreach (@{$indls{$chr}{$pos}}){
					print "$_" if($indf == 0);
					print " $_" if($indf == 1);
					$indf = 1;
				}
			}
			foreach (@{$rgs{$chr}{$pos}}){
				print "\t$_";
			}
			print "\n";
		}
	}
}

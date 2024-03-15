#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;

my %dbsnp;
my %freqs;

my $n = 3;
my $c = 5;
my $g = 1;
my $e = 6;

my $matchScore = 1;
my $mismatchScore = -3;
my $clipScore = 0;
my $gapOpenScore = -5;
my $gapExtendScore = -2;
my $read;

GetOptions("n=i" => \$n,
						"c=i" => \$c,
						"g=i" => \$g,
						"e=i" => \$e,
						"read=s" => \$read);

my $count = 0;
my $total = 0;

open(FILE, $ARGV[0]);
while(my $line = <FILE>){
	chomp $line;
	print "$line\n" if($line =~ m/^@/);
	next if($line =~ m/^@/);
	my ($rname, $flag, $chr, $pos, $mq, $cigar, @fields) = split(/\t/, $line);
	my @cigar_counts = split(/\D+/, $cigar);
	my @cigar_types = split(/\d{1,3}/, $cigar);
	shift @cigar_types;
	my $length = 0;

	for(my $i = 0; $i <= $#cigar_types; $i++){
		$length += $cigar_counts[$i] if($cigar_types[$i] ne "D");
	}
		
	my $m = $length - $n - $c - $g*$e;

	my $min_score = $matchScore*$m + $mismatchScore*$n + $clipScore*$c+$g*($gapOpenScore + $gapExtendScore * ($e-1));
	
	
	my $as_tmp;
	foreach (@fields){
		$as_tmp = $_ if($_ =~ m/^AS:i:/);
	}
	my ($tmp1, $tmp2, $as_score) = split(/:/, $as_tmp);
	if($as_score >= $min_score){
		print "$line\n";
		$count++;
	}
	$total++;
}

open(OUT, ">>$ARGV[1]");
print OUT "$read\_reads_with_high_SW_score\t$count\n";
close(OUT);


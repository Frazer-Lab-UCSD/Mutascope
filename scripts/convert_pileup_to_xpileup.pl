#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

my %SNPs = ('A',0,'T',0,'C',0,'G',0);

my $sample;
my $bed;
my $regions;
my $length=151;
my $tmp_counter=1;

GetOptions("sample=s" => \$sample,
						"bed=s" => \$bed,
						"regions=s" => \$regions,
						"tmp_counter=i" => \$tmp_counter,
						"length=i" => \$length);

my %beds;
my %regions;

open(BED, $bed);
while(my $line = <BED>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);	
	my ($chr, $bstart, $bend, $name, $score, $strand, $pstart, $pend, @other) = split(/\t/, $line);
	$beds{$chr}{$bstart}{$bend} = $pstart."\t".$pend;
}
close(BED);

open(REGIONS, $regions);
while(my $line = <REGIONS>){
	chomp $line;
	my ($chr, $bstart, $bend, $read) = split(/\t/, $line);
	$regions{$chr}{$bstart}{$bend} = $read;
}
close(REGIONS);


open(IN, $ARGV[0]);
my $out = $ARGV[0];
$out =~ s/.txt//;
open(OUT, ">$out\_$tmp_counter.txt");


READLINE: while(my $line = <IN>){
	chomp $line;
	my ($chr, $pos, $ref, $cov, $read_bases, $base_quality, $map_quality, $read_position) = split(/\t/, $line);
	my $rg;
	my $primer_one_end;
	my $primer_two_start;
LOOP:	foreach my $st (sort {$a <=> $b} keys %{$regions{$chr}}){
		foreach my $sp (sort {$a <=> $b} keys %{$regions{$chr}{$st}}){
			if($pos >= $st and $pos <= $sp){
				my $read = $regions{$chr}{$st}{$sp};
				($primer_one_end, $primer_two_start) = split(/\t/, $beds{$chr}{$st}{$sp});
				$rg = "$sample\_$read\_$chr:$st-$sp";
				last LOOP;
			}
		}
	}
	
	my %indels;
	my $indel_flag = 0;
	
	next READLINE if($pos < $primer_one_end or $pos > $primer_two_start);

	if($read_bases =~ m/[\$\^\+-]/){
		$read_bases =~ s/\^.//g; #removing the start of the read segement mark
		$read_bases =~ s/\$//g; #removing end of the read segment mark
		
		while($read_bases =~ m/[-]{1}(\d+)/) {
			$indel_flag = 1;
			my $indel_len = $1;
			my $count = 0;
			my $pat;
			while($read_bases =~ m/([-]{1}$indel_len\D{$indel_len})/g){
				$pat = $1;
				$count++;
			}
			$indels{$pat} = $count; 
			$read_bases =~ s/$pat//gi; # remove indel info from read base field
		}
		
		while ($read_bases =~ m/[\+]{1}(\d+)/) {
			$indel_flag = 1;
			my $indel_len = $1;
			my $pat;
			my $count = 0;
			while($read_bases =~ m/[\+]{1}($indel_len\D{$indel_len})/g){
				$pat = $1;
				$count++;
			}
			$indels{"+$pat"} = $count; 
			$read_bases =~ s/\+$pat//gi; # remove indel info from read base field
		}
		
	}

	my @bases = split(//, $read_bases);
	my @rp = split(/,/, $read_position);
	
	for(my $i = 0; $i <= $#rp; $i++){
		$rp[$i] = $length + 1 - $rp[$i] if($bases[$i] eq "," or $bases[$i] =~ m/(a|t|c|g)/);
	}

	
	$read_bases = join("", @bases);
	$read_position = join(",", @rp);
	
	print OUT "$chr\t$pos\t$ref\t$cov\t$rg $read_bases $base_quality $map_quality $read_position";
	if($indel_flag == 1){
		print OUT "\t$rg";
		foreach my $ind (keys %indels){
			print OUT "_$ind:$indels{$ind}";
		}
	}
	print OUT "\n";

}

close(IN);

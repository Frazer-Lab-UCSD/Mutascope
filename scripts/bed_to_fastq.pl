#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);


my $length = 151;
my $bed;
my $twoBit;
my $sub_name;
my %fastas;

GetOptions("length=i" => \$length,
						"bed=s" => \$bed,
						"twoBit=s" => \$twoBit,
						"sub_name=s" => \$sub_name);


open(BED, $bed);
open(OUT1, ">$sub_name\_tmp_file.txt");
while(my $line = <BED>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);
	my ($chr, $start, $stop, @other) = split(/\t/, $line);
	$start--;
	print OUT1 "$chr:$start-$stop\n";
}
close(OUT1);
close(BED);


`$abs_path/twoBitToFa -noMask -seqList=$sub_name\_tmp_file.txt $twoBit $sub_name.fa`;

`rm $sub_name\_tmp_file.txt`;

my $rname = "tmp";
my $fasta = "";

open(TMP, "$sub_name.fa");
open(FOR, ">$sub_name\_Forward_reads.fastq");
open(REV, ">$sub_name\_Reverse_reads.fastq");

while(my $line = <TMP>){
	chomp $line;
	if($line =~ m/^>/){
		if($rname ne "tmp"){
			my $for;
			my $rev;
			
			if(length($fasta) < $length){
				$for = $fasta;
			}else{
				$for = substr($fasta, 0, $length);
			}
			
			$fasta = reverse($fasta);
			$fasta =~ tr/ATCG/TAGC/;
			
			if(length($fasta) < $length){
				$rev = $fasta;
			}else{
				$rev = substr($fasta, 0, $length);
			}

			print FOR "\@$rname\_FOR\n";
			print FOR "$for\n";
			print FOR "+\n";
			for(my $i = 1; $i <= length($for); $i++){print FOR "B";}
			print FOR "\n";
			
			print REV "\@$rname\_REV\n";
			print REV "$rev\n";
			print REV "+\n";
			for(my $i = 1; $i <= length($rev); $i++){print REV "B";}
			print REV "\n";
			
		}
		
		$line =~ s/>//;
		my ($chr, $loc) = split(/:/, $line);
		my ($start, $stop) = split(/-/, $loc);
		$start++;
		$line = "$chr:$start-$stop";
		$rname = $line;
		$fasta = "";
		
	}else{
		$fasta .= "$line";
	}
}
my $for;
my $rev;
			
if(length($fasta) < $length){
	$for = $fasta;
}else{
	$for = substr($fasta, 0, $length);
}
$fasta = reverse($fasta);
$fasta =~ tr/ATCG/TAGC/;
			
if(length($fasta) < $length){
	$rev = $fasta;
}else{
	$rev = substr($fasta, 0, $length);
}

print FOR "\@$rname\_FOR\n";
print FOR "$for\n";
print FOR "+\n";
for(my $i = 1; $i <= length($for); $i++){print FOR "B";}
print FOR "\n";
			
print REV "\@$rname\_REV\n";
print REV "$rev\n";
print REV "+\n";
for(my $i = 1; $i <= length($rev); $i++){print REV "B";}
print REV "\n";



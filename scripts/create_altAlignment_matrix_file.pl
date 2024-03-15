#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);

sub needlemanWinsch($$);

my $twoBit;

GetOptions("twoBit=s" => \$twoBit);

my %printed;
my %mms;

open(SAM, $ARGV[0]);
while(my $line = <SAM>){
	chomp $line;
	next if($line =~ m/^@/);
	my @fields = split(/\t/, $line);
	my $seq1 = $fields[9];
	$seq1 = reverse($seq1) if($fields[1] eq 16);
	$seq1 =~ tr/ATCG/TAGC/ if($fields[1] eq 16);
	my $chr = $fields[2];
	my $pos = $fields[3];
	my ($tmp1, $tmp2) = split(/-/, $fields[0]);
	my ($chr2, $pos2) = split(/:/, $tmp1);
	
	my $xa_tmp;
	foreach (@fields){
		if($_ =~ m/^XA:Z:/){
			$xa_tmp = $_;
		}
	}
	
	my %done;
	
	if($xa_tmp){
		my ($tmp1, $tmp2, $tmp_xa) = split(/:/, $xa_tmp);
		my @xas = split(/;/, $tmp_xa);
		LOOP: foreach my $xa (@xas){
			my ($cur_chr, $cur_pos, $cur_cigar, $cur_mm) = split(/,/, $xa);
			my $strand;
			$strand = "+" if($cur_pos =~ m/^\+/);
			$strand = "-" if($cur_pos =~ m/^\-/);
			$cur_pos =~ s/\+|-//;
			next LOOP if($cur_pos >= ($pos - 2) and $cur_pos <= ($pos + 2) and $cur_chr eq $chr);
			next LOOP if($cur_pos >= ($pos2 -2) and $cur_pos <= ($pos2 + 2) and $cur_chr eq $chr2);
				
				my $seqfile = "$cur_chr:".($cur_pos - 1)."-".($cur_pos + 150);
				`$abs_path/twoBitToFa -seq=$seqfile -noMask $twoBit $ARGV[1]\/tmp_fastq_out.txt`;
				open(TMP, "$ARGV[1]/tmp_fastq_out.txt");
				my $seq2 = "";

				SUBLOOP: while(my $ln = <TMP>){
					chomp $ln;
					next SUBLOOP if($ln =~ m/^>/);
					$seq2 .= $ln;
				}
				close(TMP);
				`rm $ARGV[1]\/tmp_fastq_out.txt`;
				if($strand eq "-"){
					my $tmp = reverse($seq2);
					$tmp =~ tr/[a-z]/[A-Z]/;
					$tmp =~ tr/ACGT/TGCA/;
					$seq2 = $tmp;
				}
				my ($s1, $s2) = needlemanWinsch($seq1, $seq2);
				
				if($fields[1] == 16){
					$s1 = reverse($s1);
					$s2 = reverse($s2);
					$s1 =~ tr/ACGT/TGCA/;
					$s2 =~ tr/ACGT/TGCA/;
				}
				
				my @fasta1 = split(//, $s1);
				my @fasta2 = split(//, $s2);
				my $tmp_pos = $pos;
				for(my $i = 0; $i <= $#fasta1; $i++){
					if($fasta1[$i] eq $fasta2[$i]){
						$fasta1[$i] = ".";
						$fasta2[$i] = ".";
					}elsif($fasta1[$i] eq "-"){
						$mms{$chr}{$tmp_pos}{$fasta1[$i]}{$fasta2[$i]} = 1;
					}elsif($fasta2[$i] eq "-"){
						$mms{$chr}{$tmp_pos}{$fasta1[$i]}{$fasta2[$i]} = 1;
					}else{
						$mms{$chr}{$tmp_pos}{$fasta1[$i]}{$fasta2[$i]} = 1;
					}
					$tmp_pos++ if($fasta1[$i] ne "-");
				}
				
				$s2 = join("", @fasta2);
				$done{$seq1}{$cur_chr}{$cur_pos}=1;
				
		}
	}
}

foreach my $chr (sort keys %mms){
	foreach my $pos (sort {$a <=> $b} keys %{$mms{$chr}}){
		foreach my $ref (sort keys %{$mms{$chr}{$pos}}){
			foreach my $alt (sort keys %{$mms{$chr}{$pos}{$ref}}){
				print "$chr\t$pos\t$ref\t$alt\n";
			}
		}
	}
}





sub needlemanWinsch($$){

	my ($seq1, $seq2) = @_;
	# scoring scheme
	my $MATCH     =  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP       = -3; # -1 for any gap

	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	for(my $j = 1; $j <= length($seq1); $j++) {
		$matrix[0][$j]{score}   = $GAP*$j;
		$matrix[0][$j]{pointer} = "left";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
		$matrix[$i][0]{score}   = $GAP*$i;
		$matrix[$i][0]{pointer} = "up";
	}

	# fill
	my $max_i     = 0;
	my $max_j     = 0;
	my $max_score = 0;

	for(my $i = 1; $i <= length($seq2); $i++) {
		for(my $j = 1; $j <= length($seq1); $j++) {
			my ($diagonal_score, $left_score, $up_score);

			# calculate match score
			my $letter1 = substr($seq1, $j-1, 1);
			my $letter2 = substr($seq2, $i-1, 1);      
			if ($letter1 eq $letter2) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
			}
			else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
			}

			# calculate gap scores
			$up_score   = $matrix[$i-1][$j]{score} + $GAP;
			$left_score = $matrix[$i][$j-1]{score} + $GAP;

			# choose best score
			if ($diagonal_score >= $up_score) {
				if ($diagonal_score >= $left_score) {
					$matrix[$i][$j]{score}   = $diagonal_score;
					$matrix[$i][$j]{pointer} = "diagonal";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			} else {
				if ($up_score >= $left_score) {
					$matrix[$i][$j]{score}   = $up_score;
					$matrix[$i][$j]{pointer} = "up";
				}
				else {
					$matrix[$i][$j]{score}   = $left_score;
					$matrix[$i][$j]{pointer} = "left";
				}
			}
		}
	}

	 # trace-back

	my $align1 = "";
	my $align2 = "";

	my $j = length($seq1);
	my $i = length($seq2);

	while (1) {
		last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix
			 
		if($matrix[$i][$j]{pointer} eq "diagonal") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= substr($seq2, $i-1, 1);
			$i--; 
			$j--;
		}elsif ($matrix[$i][$j]{pointer} eq "left") {
			$align1 .= substr($seq1, $j-1, 1);
			$align2 .= "-";
			$j--;
		}elsif ($matrix[$i][$j]{pointer} eq "up") {
			$align1 .= "-";
			$align2 .= substr($seq2, $i-1, 1);
			$i--;
		}  
	}
	
	$align1 = reverse $align1;
	$align2 = reverse $align2;
	return ($align1, $align2); 

}

















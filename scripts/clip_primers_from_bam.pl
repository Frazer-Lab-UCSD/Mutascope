#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;

sub get_length($);

my %primers;

open(FILE, $ARGV[0]);
while(my $line = <FILE>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);	
	my ($chr, $bstart, $bend, $name, $score, $strand, $pstart, $pend, @other) = split(/\t/, $line);
	$pend++;
	$primers{"$chr:$bstart-$bend"} = "$pstart\t$pend";
}
close(FILE);

open(SAM, $ARGV[1]);
while(my $line = <SAM>){
	chomp $line;
	if($line =~ m/^@/){
		print "$line\n";
		next;
	}
	my @fields = split(/\t/, $line);
	my $rname = $fields[0];
	my $flag = $fields[1];
	my $chr = $fields[2];
	my $pos = $fields[3];
	my $cigar = $fields[5];
	my $length = get_length($cigar);

	my $rgroup;
	foreach (@fields){
		if($_ =~ m/^RG:Z:/){
			$rgroup = $_;
			$rgroup =~ s/RG:Z://;
		}
	}
	my ($samp, $strand, $loc) = split(/_/, $rgroup);
	my ($c1, $tmp) = split(/:/, $loc);
	my ($start, $end) = split(/-/, $tmp);
	my ($pstart, $pend) = split(/\t/, $primers{$loc});
	
	my @cigar_counts = split(/\D+/, $cigar);
	my @cigar_types = split(/\d{1,3}/, $cigar);
	shift @cigar_types;
	my @cigar_array;
	my @original_array;
	my $del_count = 0;
	my $ins_count = 0;
	
	for(my $i = 0; $i <= $#cigar_counts; $i++){
		$del_count += $cigar_counts[$i] if($cigar_types[$i] eq "D");
		$ins_count += $cigar_counts[$i] if($cigar_types[$i] eq "I");
		for(my $j = 1; $j <= $cigar_counts[$i]; $j++){
			push(@cigar_array, $cigar_types[$i]);
			push(@original_array, $cigar_types[$i]);
		}
	}

####  CHECK FOR DELETIONS!!!!
	
#print "@cigar_array\n";

#print "$rgroup\t$pos\t$cigar\t$pstart\t$pend\n";

	if($cigar =~ m/^\d+S/){
		my ($len, @tmps) = split(/S/, $cigar);
		$pos -= $len;
	}

	my $read_end = $pos+$length+$del_count-$ins_count;
	my $read_len = $length+$del_count;
	
	my $del_in_primer = 0;
	my $tmp_count = 0;
	if($pos < $pstart){
		for(my $i = $pos; $i < $pstart; $i++){
			if($cigar_array[$tmp_count] eq "I"){
				$i--;
				$pos--;
			}
			$del_in_primer++ if($cigar_array[$tmp_count] eq "D");
			$cigar_array[$tmp_count] = "S";
			$tmp_count++;
		}
	}
	
	if($pend < $read_end){
		$tmp_count = 1;
		for(my $j = 1; $j+$pend <= $read_end; $j++){
			$j-- if($cigar_array[$read_len-$tmp_count] eq "I");
			$cigar_array[$read_len-$tmp_count] = "S";
			$tmp_count++;
		}
	}
	
	my $new_cigar = "";
	my $count = 0;
	my $cur_type = $cigar_array[0];
	for(my $i = 0; $i<=$#cigar_array;$i++){
		if($cigar_array[$i] eq "S" and $original_array[$i] eq "D"){
			
		}else{
			if($cur_type eq $cigar_array[$i]){
				$count++;
			}else{
				$new_cigar .= ($count.$cur_type);
				$cur_type = $cigar_array[$i];
				$count = 1;
			}
		}
	}
	$new_cigar .= ($count.$cur_type);
#print "@cigar_array\n";
	if($new_cigar =~ m/^\d+S/){
		my ($len, @tmps) = split(/S/, $new_cigar);
		$pos += $len;
	}
	
	$fields[3] = $pos + $del_in_primer;
	$fields[5] = $new_cigar;
	$line = join("\t", @fields);
	print "$line\n";
	
	
}
close(SAM);


sub get_length($){
	my ($cigar) = @_;
	my @cigar_counts = split(/\D+/, $cigar);
	my @cigar_types = split(/\d{1,3}/, $cigar);
	shift @cigar_types;
	my $length = 0;

	for(my $i = 0; $i <= $#cigar_types; $i++){
		$length += $cigar_counts[$i] if($cigar_types[$i] ne "D");
	}
	return($length);
}




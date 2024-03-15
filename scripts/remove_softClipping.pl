#! /usr/bin/perl -w

use warnings;
use strict;

open(SAM, $ARGV[0]);
while(my $line = <SAM>){
	chomp $line;
	if($line =~ m/^@/){
		print "$line\n";
		next;
	}
	my @fields = split(/\t/, $line);
	my $cigar = $fields[5];

	if($cigar =~ m/S/){		
		if($fields[1] == 0){
			if($cigar =~ m/\d{1,3}S$/){

				my @cigar_counts = split(/\D+/, $cigar);
				my @cigar_types = split(/\d{1,3}/, $cigar);
				shift @cigar_types;
				my @cigar_array;
				for(my $i = 0; $i <= $#cigar_counts; $i++){
					for(my $j = 1; $j <= $cigar_counts[$i]; $j++){
						push(@cigar_array, $cigar_types[$i]);
					}
				}		
				my $cigar_flag = 0;
				for(my $i = 0; $i <= $#cigar_array; $i++){
					if($cigar_array[$i] eq "S" and $cigar_flag == 1){
						$cigar_array[$i] = "M";
					}
					$cigar_flag = 1 if($cigar_array[$i] eq "M");
				}


				my $new_cigar = "";
				my $count = 0;
				my $cur_type = $cigar_array[0];
				for(my $i = 0; $i<=$#cigar_array;$i++){
					if($cur_type eq $cigar_array[$i]){
						$count++;
					}else{
						$new_cigar .= ($count.$cur_type);
						$cur_type = $cigar_array[$i];
						$count = 1;
					}
				}
				$new_cigar .= ($count.$cur_type);
				
				$fields[5] = $new_cigar;
				$line = join("\t", @fields);
				print "$line\n";
				next;
			}else{
				print "$line\n";
				next;
			}
		}else{
			if($cigar =~ m/^\d{1,3}S/){				
				my ($len, @other) = split(/S/, $cigar);

				$fields[3]-= $len;

				my @cigar_counts = split(/\D+/, $cigar);
				my @cigar_types = split(/\d{1,3}/, $cigar);
				shift @cigar_types;
				my @cigar_array;
				for(my $i = 0; $i <= $#cigar_counts; $i++){
					for(my $j = 1; $j <= $cigar_counts[$i]; $j++){
						push(@cigar_array, $cigar_types[$i]);
					}
				}
				my $cigar_flag = 1;
				for(my $i = 0; $i <= $#cigar_array; $i++){
					$cigar_flag = 0 if($cigar_array[$i] eq "M");
					if($cigar_array[$i] eq "S" and $cigar_flag == 1){
						$cigar_array[$i] = "M";
					}
				}


				my $new_cigar = "";
				my $count = 0;
				my $cur_type = $cigar_array[0];
				for(my $i = 0; $i<=$#cigar_array;$i++){
					if($cur_type eq $cigar_array[$i]){
						$count++;
					}else{
						$new_cigar .= ($count.$cur_type);
						$cur_type = $cigar_array[$i];
						$count = 1;
					}
				}
				$new_cigar .= ($count.$cur_type);
				
				$fields[5] = $new_cigar;
				$line = join("\t", @fields);
				print "$line\n";
				next;
			}else{
				print "$line\n";
				next;
			}
		}
	}else{
		print "$line\n";
		next;
	}
}

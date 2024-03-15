#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);


sub snp_indel_cluster($$);
sub indel_cluster($$);
sub add_filter($$);

my %vcf_info;
my %indels;
my %snps;
my %vcf_order;

my $i = 1;
my $header = "";

my $window_indel_size = 10;
my $snp_cluster_window = 10;
my $snp_cluster_count = 3;
my $indel_cluster_window = 10;
my $indel_cluster_count = 2;

my $vcf;

GetOptions("vcf=s" => \$vcf);


open(OUT, ">$vcf\_tmp_bq.txt");

my $header_flag = 0;

open(VCF, $vcf);
while(my $line = <VCF>){
	chomp $line;
	if($line =~ m/^#/){
		#ADD FILTERS
		if($line =~ m/^##FILTER/ and $header_flag == 0){
			$header .= "##FILTER=<ID=SBFILT,Description=\"Alternate allele seen only on 1 strand while reference allele is seen on both\">\n";
			$header .= "##FILTER=<ID=RGBF,Description=\"The alternate allele has a Read-Group Bias\">\n";
			$header .= "##FILTER=<ID=DSF,Description=\"The First and Second allele frequencies are too close\">\n";
			$header .= "##FILTER=<ID=MAMQ,Description=\"Low average Mapping Quality score\">\n";
			$header .= "##FILTER=<ID=MFEP,Description=\"Low Fisher-exact P-value for Somatic variant\">\n";
			$header .= "##FILTER=<ID=MBPF,Description=\"Low Binomial P-value to call the position a variant\">\n";
			$header .= "##FILTER=<ID=MAABQ,Description=\"Low average alternate allele Base Quality score\">\n";
			$header .= "##FILTER=<ID=MAAMQ,Description=\"Low average alternate allele Mapping Quality score\">\n";
			$header .= "##FILTER=<ID=CNIF,Description=\"The SNP is within $window_indel_size"."bp of an Indel\">\n";
			$header .= "##FILTER=<ID=SNPCF,Description=\"The SNP is within $snp_cluster_window"."bp of ".($snp_cluster_count-1)." other SNPs\">\n";
			$header .= "##FILTER=<ID=ICF,Description=\"The Indel is within $indel_cluster_window"."bp of ".($indel_cluster_count-1)." other Indel\">\n";

			$header_flag = 1;
		}
		$header .= "$line\n";
		next;
	}
	

	my ($chr, $pos, $rsid, $ref, $alt, $qual, $filt, $info, $format, @samples) = split(/\t/, $line);
	my @infos = split(/;/, $info);
	my $vt;

	foreach (@infos){
		if($_ =~ m/VT=/){
			$vt = $_;
			$vt =~ s/VT=//;
		}
	}

	my @nother = split(/:/, $samples[0]);
	my @tother = split(/:/, $samples[1]);
		
	my $bq2;
	my $ad;
	my $af;
	my $ss;
	my $rgb;
	my @formats = split(/:/, $format);
	for(my $i = 0; $i <= $#formats; $i++){
		if($formats[$i] =~ m/BQ2/){
			$bq2 = $i;
		}
		if($formats[$i] =~ m/AD/){
			$ad = $i;
		}
		if($formats[$i] =~ m/AF/){
			$af = $i;
		}
		if($formats[$i] =~ m/SS/){
			$ss = $i;
		}
		if($formats[$i] =~ m/RGB/){
			$rgb = $i;
		}		
	}
	
	my $tss = $tother[$ss];
	
	if($vt ne "SNP"){
		$indels{$chr}{$pos} = $filt;
	}else{
		$snps{$chr}{$pos} = $filt;
	}
	
	$vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} = $vt;
	$vcf_info{$chr}{$pos}{$ref}{$alt}{"LINE"} = $line;
	$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = $filt;
	
	$vcf_order{$i} = "$chr\t$pos\t$ref\t$alt";
	$i++;
	
	if($vt eq "SNP"){
		my ($nrefbq, $naltbq) = split(/,/, $nother[$bq2]);
		my ($trefbq, $taltbq) = split(/,/, $tother[$bq2]);	
		
		$vcf_info{$chr}{$pos}{$ref}{$alt}{"ABQS"} = $taltbq;
		print OUT "$taltbq\n";
	}	

}

`R --slave --args $vcf\_tmp_bq.txt $vcf\_tmp_minbq.txt < $abs_path/find_minaltBQ_with_mclust.r`;

my $min_aabq;
open(IN, "$vcf\_tmp_minbq.txt");
while(my $line = <IN>){
	chomp $line;
	$min_aabq = $line;
}


if($min_aabq > 25){
	$min_aabq = 25;
}


close(IN);
close(OUT);
close(VCF);
`rm $vcf\_tmp_minbq.txt`;
`rm $vcf\_tmp_bq.txt`;


print "$header";

### Fist filter based on Binomial P-value and Fisher-test
foreach my $l (sort {$a <=> $b} keys %vcf_order){
	my ($chr, $pos, $ref, $alt) = split(/\t/, $vcf_order{$l});
	my $vcf_line = $vcf_info{$chr}{$pos}{$ref}{$alt}{"LINE"};
	my @fields = split(/\t/, $vcf_line);
	my ($c, $p, $id, $rf, $at, $qual, $filt, $info, $format, $sample1, $sample2) = split(/\t/, $vcf_line);
	$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = "PASS" if ($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} eq "MBPF"); # need to ease the old MBPF flag in case it exists
	$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "MBPF") if($qual < 20);
	
	my @infos = split(/;/, $info);
	my $fep;
	foreach (@infos){
		if($_ =~ m/FEP=/){
			$fep = $_;
			$fep =~ s/FEP=//;
		}
	}
	
	
	my $ss;
	my $rgb;
	my @formats = split(/:/, $format);
	for(my $i = 0; $i <= $#formats; $i++){
		if($formats[$i] =~ m/SS/){
			$ss = $i;
		}
		if($formats[$i] =~ m/RGB/){
			$rgb = $i;
		}	
	}

	my @tother = split(/:/, $sample2);
	my @nother = split(/:/, $sample1);
	my $tss = $tother[$ss];
	my $nss = $tother[$ss];
	my $trgb = $tother[$rgb] if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP");
	my $nrgb = $nother[$rgb] if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP");
	
	$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "MFEP") if($tss == 2 and $fep > 0.0005);


	if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP"){
		if($nss==1){
			if($nrgb ne "." and $trgb ne "."){
				$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "RGBF") if(($nrgb < 10**-15 or $trgb <10**-15));
			}elsif($nrgb ne "."){
				$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "RGBF") if(($nrgb < 10**-15));
			}elsif($trgb ne "."){
				$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "RGBF") if(($trgb <10**-15));	
			}
		}else{
			if($trgb ne "."){
				$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "RGBF") if(($trgb <10**-15));	
			}
		}
	}


	if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} ne "SNP"){
		$indels{$chr}{$pos} = $vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"};
	}else{
		$snps{$chr}{$pos} = $vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"};
	}
	
	
}


#Run INDEL and SNP clustering filter
foreach my $l (sort {$a <=> $b} keys %vcf_order){
	my ($chr, $pos, $ref, $alt) = split(/\t/, $vcf_order{$l});
	my $vcf_line = $vcf_info{$chr}{$pos}{$ref}{$alt}{"LINE"};
	my @fields = split(/\t/, $vcf_line);
	my ($c, $p, $id, $rf, $at, $qual, $f, $info, $format, $sample1, $sample2) = split(/\t/, $vcf_line);
	
	
	if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP"){
		my ($indel_flag, $clust_flag) = snp_indel_cluster($chr, $pos);

		if($indel_flag == 1){
			$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "CNIF");
		}
		if($clust_flag == 1){
			$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "SNPCF");
		}
	}else{
		my ($indel_clust_flag) = indel_cluster($chr, $pos);
		$vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"} = add_filter($vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"}, "ICF") if($indel_clust_flag == 1);		
	}


}



#Add remaining FILTERS!!!
foreach my $l (sort {$a <=> $b} keys %vcf_order){
	my ($chr, $pos, $ref, $alt) = split(/\t/, $vcf_order{$l});
	my $vcf_line = $vcf_info{$chr}{$pos}{$ref}{$alt}{"LINE"};
	my @fields = split(/\t/, $vcf_line);
	my ($c, $p, $id, $rf, $at, $qual, $filt, $info, $format, $sample1, $sample2) = split(/\t/, $vcf_line);

	$filt = $vcf_info{$chr}{$pos}{$ref}{$alt}{"FILT"};

#### Hard Cutoff Filters ####
	my $vt = $vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"};
	my @infos = split(/;/, $info);
	my $mq;
	my $fep;
	foreach (@infos){
		if($_ =~ m/FEP=/){
			$fep = $_;
			$fep =~ s/FEP=//;
		}
		if($_ =~ m/MQ=/){
			$mq = $_;
			$mq =~ s/MQ=//;
		}
	}
	
	my $af;
	my $ss;
	my $ds;
	my $rgb;
	my $dp4;
	my $mq2;
	my @formats = split(/:/, $format);
	for(my $i = 0; $i <= $#formats; $i++){
		if($formats[$i] =~ m/DP4/){
			$dp4 = $i;
		}
		if($formats[$i] =~ m/AF/){
			$af = $i;
		}
		if($formats[$i] =~ m/SS/){
			$ss = $i;
		}
		if($formats[$i] =~ m/DS/){
			$ds = $i;
		}
		if($formats[$i] =~ m/RGB/){
			$rgb = $i;
		}		
		if($formats[$i] =~ m/MQ2/){
			$mq2 = $i;
		}
	}
	my @nother = split(/:/, $sample1);
	my @tother = split(/:/, $sample2);

	my ($nrfor, $nrrev, $nafor, $narev) = split(/,/, $nother[$dp4]);
	my ($trfor, $trrev, $tafor, $tarev) = split(/,/, $tother[$dp4]);
	
	my $nss = $nother[$ss];
	my $tss = $tother[$ss];

	my $trgb = $tother[$rgb] if($vt eq "SNP");
	my $nrgb = $nother[$rgb] if($vt eq "SNP");

	my $tds = $tother[$ds];
	
	if($vt eq "SNP"){
		if($tss == 2){
			if($trfor > 0 and $trrev > 0){
				$filt = add_filter($filt, "SBFILT") if($tafor == 0 or $tarev == 0);
			}
		}
		if($nss == 1){
			if($nrfor > 0 and $nrrev > 0){
				$filt = add_filter($filt, "SBFILT") if($nafor == 0 or $narev == 0);
			}
		}
	}
	
	
	$filt = add_filter($filt, "DSF") if($tss == 2 and $tds ne "." and $tds <= 0.5);
	$filt = add_filter($filt, "MAMQ") if($mq ne "." and $mq < 15);

	my ($trefmq, $taltmq) = split(/,/, $tother[$mq2]) if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP");	
	
	$filt = add_filter($filt, "MAAMQ") if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP" and $taltmq < 10);
	
	if($vcf_info{$chr}{$pos}{$ref}{$alt}{"VT"} eq "SNP" and $vcf_info{$chr}{$pos}{$ref}{$alt}{"ABQS"} < $min_aabq){
		$filt = add_filter($filt, "MAABQ");
	}

	
	$fields[6] = $filt;
	my $line = join("\t", @fields);

	print "$line\n";
}




sub snp_indel_cluster($$){
	my ($chr, $pos) = @_;
	my $min_ind_window = $pos - $window_indel_size;
	my $max_ind_window = $pos + $window_indel_size;
	
	my $min_snp_clust_window = $pos - $snp_cluster_window;
	my $max_snp_clust_window = $pos + $snp_cluster_window;
	
	my $min_window = $min_snp_clust_window;
	my $max_window = $max_snp_clust_window;
	
	if($min_ind_window < $min_snp_clust_window){
		$min_window = $min_ind_window;
	}
	if($max_ind_window > $max_snp_clust_window){
		$min_window = $max_ind_window;
	}
	
	my $indel_flag = 0;
	my $clust_flag = 0;
	
	my @positions;
	
	for(my $i = $min_window; $i <= $max_window; $i++){
		if(exists($indels{$chr}{$i}) and ($indels{$chr}{$i} eq "PASS" or ($indels{$chr}{$i} =~ m/HRUN/ and !($indels{$chr}{$i} =~ m/MBPF/)))){
			$indel_flag = 1 if($i >= $min_ind_window and $i <= $max_ind_window);
		}
		
		if(exists($snps{$chr}{$i}) and $snps{$chr}{$i} eq "PASS"){
			push(@positions, $i);
			while($i - $positions[0] > $snp_cluster_window){
				shift(@positions);
			}
			$clust_flag = 1 if($#positions >= ($snp_cluster_count-1) and ($positions[0] >= $min_snp_clust_window and $positions[$#positions] <= $max_snp_clust_window));
		}
		
	}
	return($indel_flag, $clust_flag);
}


sub indel_cluster($$){
	my ($chr, $pos) = @_;
	my $min_window = $pos - $indel_cluster_window;
	my $max_window = $pos + $indel_cluster_window;
	my $clust_flag = 0;
	my @positions;

	for(my $i = $min_window; $i <= $max_window; $i++){
		
		if(exists($indels{$chr}{$i})){
			push(@positions, $i);
			while($i - $positions[0] > $indel_cluster_window){
				shift(@positions);
			}
			$clust_flag = 1 if($#positions >= ($indel_cluster_count-1));
		}
	}
	return($clust_flag);	
}

sub add_filter($$){
	my ($filt, $name) = @_;
	if($filt eq "PASS"){
		return($name);
	}else{
		$filt .= ",$name";
		return($filt);
	}
}










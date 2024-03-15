#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);


sub print_header($$);
sub get_error_rates($\%);
sub get_sub_type($$);
sub get_black_list($\%);
sub get_dbsnp($\%\%);
sub get_bed($\%);

sub mean(\%\%\%$);
sub avg_error(\%\%);
sub get_pvalue($$$$);
sub rms(@);
sub get_sample_vcf($\%\%@);
sub set_chrom_order(\%);
sub find_indel_repeat($$);
sub calc_min_cov($$\%);
sub get_indel_information($$$$$$$);
sub add_filter($$$\%);
sub get_prob_b_given_a($$$);
sub convert_from_phred($);
sub log10($);

my $BQ = 0;
my $normal;
my $e = 0.005;
my $saaf = 0.5;
my $min_cov = 10;
my $min_perc_cov = 0.005;
my $gaaf = 15;
my $tumor;
my $MQ = 0;
my $het = 0.001;
my $hom_ref;
my $hom_alt;
my $indel_het = 0.000125;
my $indel_hom_ref;
my $indel_hom_alt;

my $error_rates;
my $blacklist;
my $dbsnp;
my $nxpileup;
my $txpileup;
my $pdir;
my $twoBit;
my $noBlacklist;
my $bed;

my @chroms = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY");

my @genotypes = ("AA", "AT", "AG", "AC", "TT", "TC", "TG", "CC", "CG", "GG");

my %normal_fields;
my %tumor_fields;
my %norm_indel_fields;
my %tum_indel_fields;

my %error_rate;
my %skip;
my %ref_error;
my %dbsnp;
my %chrom_order;

my %phred_to_pval;
my %locations;

GetOptions("error_rates=s" => \$error_rates,
						"blacklist=s" => \$blacklist,
						"dbsnp=s" => \$dbsnp,
						"nxpileup=s" => \$nxpileup,
						"txpileup=s" => \$txpileup,
						"pdir=s" => \$pdir,
						"twoBit=s" => \$twoBit,
						"BQ=i" => \$BQ,
						"MQ=i" => \$MQ,
						"normal=s" => \$normal,
						"tumor=s" => \$tumor,
						"saaf=f" => \$saaf,
						"min_cov=i" => \$min_cov,
						"e=f" => \$e,
						"het=f" => \$het,
						"indel_het=f" => \$indel_het,
						"gaaf=f" => \$gaaf,
						"min_perc_cov=f" => \$min_perc_cov,
						"noBlacklist" => \$noBlacklist,
						"bed=s" => \$bed);

$hom_ref = 1-3*$het;
$hom_alt = $het/2;

$indel_hom_ref = 1-3*$indel_het;
$indel_hom_alt = $indel_het/2;

for(my $i = 0; $i<=60; $i++){
	$phred_to_pval{$i} = convert_from_phred($i);
}

#################
##  Fills in the hash to sort by chrom order instead of alph.
#################
set_chrom_order(%chrom_order);

#################
##  Gets the calculated Error Rates
#################
get_error_rates($error_rates, %error_rate);

#################
##  Sets error rates for positions that have 0 error rate using same position and strand
#################
avg_error(%error_rate, %ref_error);

#################
##  Get BLACK listed mutations to skip
#################
get_black_list($blacklist, %skip) if(!$noBlacklist);

#################
##  Get bed info
#################
get_bed($bed, %locations);

#################
##  Get dbSNP info
#################
get_dbsnp($dbsnp, %dbsnp, %locations);


#################
##  Open tmp stats files
#################
open(BINOMT, ">$pdir/intermediate/$normal\_binom_tumor_tmp.txt");
open(BINOMN, ">$pdir/intermediate/$normal\_binom_normal_tmp.txt");
open(FISH, ">$pdir/intermediate/$normal\_fisher_tmp.txt");
open(BIAS, ">$pdir/intermediate/$normal\_rgbias_fisher_tmp.txt");
open(TBIAS, ">$pdir/intermediate/$tumor\_rgbias_fisher_tmp.txt");
open(INDFISH, ">$pdir/intermediate/$normal\_indel_fisher_test.txt");
open(INDBT, ">$pdir/intermediate/$normal\_indel_binom_tumor_tmp.txt");
open(INDBN, ">$pdir/intermediate/$normal\_indel_binom_normal_tmp.txt");


#################
##  Open the Tumor and Normal files
#################
if($nxpileup =~ m/\.gz$/){
	open(NFILE, "gunzip -c $nxpileup |") || die "Can't open pipe to $nxpileup";
}else{
	open(NFILE, $nxpileup) || die "Can't open $nxpileup";
}
if($txpileup =~ m/\.gz$/){
	open(TFILE, "gunzip -c $txpileup |") || die "Can't open pipe to $txpileup";
}else{
	open(TFILE, $txpileup) || die "Can't open $txpileup";
}

my $nline = <NFILE>;
my $tline = <TFILE>;
chomp $nline;
chomp $tline;
my @nfields = split(/\t/, $nline);
my @tfields = split(/\t/, $tline);


LINES: while(){
	
	if($nfields[0] eq $tfields[0] and $nfields[1] == $tfields[1]){
		my $chr = $nfields[0];
		my $pos = $nfields[1];
		
		get_sample_vcf("norm" ,%normal_fields, %norm_indel_fields, @nfields);
		get_sample_vcf("tum", %tumor_fields, %tum_indel_fields, @tfields);
		
		my ($nrf, $nrr, $naf, $nar) = split(/,/, $normal_fields{$chr}{$pos}{"FILTDP4"});
		my ($trf, $trr, $taf, $tar) = split(/,/, $tumor_fields{$chr}{$pos}{"FILTDP4"});
		
		my $talt = $taf + $tar;
		my $nalt = $naf + $nar;
		my $nref = $nrf + $nrr;
		my $tref = $trf + $trr;		
		
		if($normal_fields{$chr}{$pos}{"GT"} ne "0/0" or $tumor_fields{$chr}{$pos}{"GT"} ne "0/0"){
			$tumor_fields{$chr}{$pos}{"ER"} = $e if($tumor_fields{$chr}{$pos}{"ER"} == 0);
			$normal_fields{$chr}{$pos}{"ER"} = $e if($normal_fields{$chr}{$pos}{"ER"} == 0);
			
			print BINOMT "$chr\t$pos\t$talt\t".$tumor_fields{$chr}{$pos}{"FDP"}."\t".$tumor_fields{$chr}{$pos}{"ER"}."\n";
			print BINOMN "$chr\t$pos\t$nalt\t".$normal_fields{$chr}{$pos}{"FDP"}."\t".$normal_fields{$chr}{$pos}{"ER"}."\n";
			print FISH "$chr\t$pos\t$nalt\t$nref\t$talt\t$tref\n";
			print BIAS $normal_fields{$chr}{$pos}{"RGB"}."\n";
			print TBIAS $tumor_fields{$chr}{$pos}{"RGB"}."\n";
		}
		
		if(exists($norm_indel_fields{$chr}{$pos}{"REF"}) and exists($tum_indel_fields{$chr}{$pos}{"REF"})){
			if($norm_indel_fields{$chr}{$pos}{"REF"} ne "." or $tum_indel_fields{$chr}{$pos}{"REF"} ne "."){
				print INDBT "$chr\t$pos\t".$tum_indel_fields{$chr}{$pos}{"COV"}."\t".$tum_indel_fields{$chr}{$pos}{"DP"}."\t$e\n";
				print INDBN "$chr\t$pos\t".$norm_indel_fields{$chr}{$pos}{"COV"}."\t".$norm_indel_fields{$chr}{$pos}{"DP"}."\t$e\n";
				print INDFISH "$chr\t$pos\t".$norm_indel_fields{$chr}{$pos}{"COV"}."\t".$norm_indel_fields{$chr}{$pos}{"RC"}."\t".$tum_indel_fields{$chr}{$pos}{"COV"}."\t".$tum_indel_fields{$chr}{$pos}{"RC"}."\n";
			}
		}
		
		$nline = <NFILE>;
		$tline = <TFILE>;
		last LINES unless defined $nline;
		last LINES unless defined $tline;
		chomp $nline;
		chomp $tline;
		@nfields = split(/\t/, $nline);
		@tfields = split(/\t/, $tline);
		next LINES;
	}
	
	if($nfields[0] ne $tfields[0]){
		while($nfields[0] ne $tfields[0]){
			if($chrom_order{$nfields[0]} < $chrom_order{$tfields[0]}){
				$nline = <NFILE>;
				last LINES unless defined $nline;
				chomp $nline;
				@nfields = split(/\t/, $nline);					
			}else{
				$tline = <TFILE>;
				last LINES unless defined $tline;
				chomp $tline;
				@tfields = split(/\t/, $tline);
			}
		}
		next LINES;
	}
	
	if($nfields[1] > $tfields[1]){
		while($nfields[1] > $tfields[1]){
			$tline = <TFILE>;
			last LINES unless defined $tline;
			chomp $tline;
			@tfields = split(/\t/, $tline);
			next LINES if($nfields[0] ne $tfields[0]);
		}
		next LINES;
	}
	
	if($nfields[1] < $tfields[1]){
		while($nfields[1] < $tfields[1]){
			$nline = <NFILE>;
			last LINES unless defined $nline;
			chomp $nline;
			@nfields = split(/\t/, $nline);
			next LINES if($nfields[0] ne $tfields[0]);
		}
		next LINES;
	}
}


close(BINOMT);
close(BINOMN);
close(FISH);
close(INDFISH);
close(INDBT);
close(INDBN);
close(BIAS);
close(TBIAS);




####################
##  Get Min Coverage for Tumor and Normal
####################

my $min_norm_cov = calc_min_cov($min_cov, $min_perc_cov, %normal_fields);
my $min_tum_cov = calc_min_cov($min_cov, $min_perc_cov, %tumor_fields);

#################
##  Prints the VCF Header
#################
print_header($min_norm_cov, $min_tum_cov);

my $tmp_file = "$pdir/intermediate/$normal\_tmp_r_file.txt";
`R --slave --args $pdir/intermediate/$normal\_binom_normal_tmp.txt $tmp_file < $abs_path/calc_pval.r` if(!(-z "$pdir/intermediate/$normal\_binom_normal_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $nalt, $ncov, $error, $pval) = split(/\t/, $line);
	$normal_fields{$chr}{$pos}{"P"} = int($pval + 0.5);
	$normal_fields{$chr}{$pos}{"P"} = 100000 if($normal_fields{$chr}{$pos}{"P"} eq "inf");
}
close(TMPF);

`R --slave --args $pdir/intermediate/$normal\_binom_tumor_tmp.txt $tmp_file < $abs_path/calc_pval.r` if(!(-z "$pdir/intermediate/$normal\_binom_tumor_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $nalt, $ncov, $error, $pval) = split(/\t/, $line);
	$tumor_fields{$chr}{$pos}{"P"} = int($pval + 0.5);
	$tumor_fields{$chr}{$pos}{"P"} = 100000 if($tumor_fields{$chr}{$pos}{"P"} eq "inf");
}
close(TMPF);

`R --slave --args $pdir/intermediate/$normal\_fisher_tmp.txt $tmp_file < $abs_path/calc_fisher_pval.r` if(!(-z "$pdir/intermediate/$normal\_fisher_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $nalt, $nref, $talt, $tref, $pval_less, $pval_greater) = split(/\t/, $line);
	$normal_fields{$chr}{$pos}{"FEP"} = sprintf("%.3g",($pval_greater));
	$tumor_fields{$chr}{$pos}{"FEP"} = sprintf("%.3g",($pval_less));
}
close(TMPF);

`R --slave --args $pdir/intermediate/$normal\_indel_fisher_test.txt $tmp_file < $abs_path/calc_fisher_pval.r` if(!(-z "$pdir/intermediate/$normal\_indel_fisher_test.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $nalt, $nref, $talt, $tref, $pval_less, $pval_greater) = split(/\t/, $line);
	$norm_indel_fields{$chr}{$pos}{"FEP"} = sprintf("%.3g",($pval_greater));
	$tum_indel_fields{$chr}{$pos}{"FEP"} = sprintf("%.3g",($pval_less));
}
close(TMPF);

`R --slave --args $pdir/intermediate/$normal\_indel_binom_tumor_tmp.txt $tmp_file < $abs_path/calc_pval.r` if(!(-z "$pdir/intermediate/$normal\_indel_binom_tumor_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $nalt, $ncov, $error, $pval) = split(/\t/, $line);
	$tum_indel_fields{$chr}{$pos}{"P"} = int($pval + 0.5);
	$tum_indel_fields{$chr}{$pos}{"P"} = 100000 if($tum_indel_fields{$chr}{$pos}{"P"} eq "inf");
}
close(TMPF);

`R --slave --args $pdir/intermediate/$normal\_indel_binom_normal_tmp.txt $tmp_file < $abs_path/calc_pval.r` if(!(-z "$pdir/intermediate/$normal\_indel_binom_normal_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $nalt, $ncov, $error, $pval) = split(/\t/, $line);
	$norm_indel_fields{$chr}{$pos}{"P"} = int($pval + 0.5);
	$norm_indel_fields{$chr}{$pos}{"P"} = 100000 if($norm_indel_fields{$chr}{$pos}{"P"} eq "inf");
}
close(TMPF);

`R --slave --args $pdir/intermediate/$normal\_rgbias_fisher_tmp.txt $tmp_file < $abs_path/calc_rgbias_pval.r` if(!(-z "$pdir/intermediate/$normal\_rgbias_fisher_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $pval) = split(/\t/, $line);
	$pval = "." if($pval eq "NA");
	$normal_fields{$chr}{$pos}{"RGB"} = sprintf("%.3g",($pval)) if($pval ne ".");
	$normal_fields{$chr}{$pos}{"RGB"} = $pval if($pval eq ".");
}
close(TMPF);

`R --slave --args $pdir/intermediate/$tumor\_rgbias_fisher_tmp.txt $tmp_file < $abs_path/calc_rgbias_pval.r` if(!(-z "$pdir/intermediate/$tumor\_rgbias_fisher_tmp.txt"));
open(TMPF, $tmp_file);
while(my $line = <TMPF>){
	chomp $line;
	my ($chr, $pos, $pval) = split(/\t/, $line);
	$pval = "." if($pval eq "NA");
	$tumor_fields{$chr}{$pos}{"RGB"} = sprintf("%.3g",($pval)) if($pval ne ".");
	$tumor_fields{$chr}{$pos}{"RGB"} = $pval if($pval eq ".");
}
close(TMPF);


`rm $tmp_file`;
`rm $pdir/intermediate/$normal\_fisher_tmp.txt`;
`rm $pdir/intermediate/$normal\_binom_tumor_tmp.txt`;
`rm $pdir/intermediate/$normal\_binom_normal_tmp.txt`;
`rm $pdir/intermediate/$normal\_indel_fisher_test.txt`;
`rm $pdir/intermediate/$normal\_indel_binom_tumor_tmp.txt`;
`rm $pdir/intermediate/$normal\_indel_binom_normal_tmp.txt`;
`rm $pdir/intermediate/$normal\_rgbias_fisher_tmp.txt`;
`rm $pdir/intermediate/$tumor\_rgbias_fisher_tmp.txt`;



foreach my $chr (@chroms){
POS_LOOP:	foreach my $pos (sort {$a <=> $b} keys %{$normal_fields{$chr}}){

#INDEL CHECK
		if(exists($norm_indel_fields{$chr}{$pos}{"REF"}) and exists($tum_indel_fields{$chr}{$pos}{"REF"})){

			if(($norm_indel_fields{$chr}{$pos}{"REF"} ne "." or $tum_indel_fields{$chr}{$pos}{"REF"} ne ".") and ($norm_indel_fields{$chr}{$pos}{"GT"} ne "0/0" or $tum_indel_fields{$chr}{$pos}{"GT"} ne "0/0")){
				
				if($norm_indel_fields{$chr}{$pos}{"DP"} < $min_norm_cov or $tum_indel_fields{$chr}{$pos}{"DP"} < $min_tum_cov){
					add_filter($chr, $pos, "MCF", %norm_indel_fields);
					add_filter($chr, $pos, "MCF", %tum_indel_fields);
				}
				
				if(($norm_indel_fields{$chr}{$pos}{"REF"} ne "." and $tum_indel_fields{$chr}{$pos}{"REF"} ne ".") and ($norm_indel_fields{$chr}{$pos}{"REF"} ne $tum_indel_fields{$chr}{$pos}{"REF"} or $norm_indel_fields{$chr}{$pos}{"ALT"} ne $tum_indel_fields{$chr}{$pos}{"ALT"}) and ($norm_indel_fields{$chr}{$pos}{"GT"} ne "0/0" and $tum_indel_fields{$chr}{$pos}{"GT"} ne "0/0")){
					
					if(($norm_indel_fields{$chr}{$pos}{"REF"} ne $tum_indel_fields{$chr}{$pos}{"REF"}) and ($norm_indel_fields{$chr}{$pos}{"ALT"} ne $tum_indel_fields{$chr}{$pos}{"ALT"})){
						my $long_ref;
						my $short_alt;
						my $long_alt;
						
						if(length($norm_indel_fields{$chr}{$pos}{"REF"}) > length($tum_indel_fields{$chr}{$pos}{"REF"})){
							$long_ref = $norm_indel_fields{$chr}{$pos}{"REF"};
							$short_alt = $norm_indel_fields{$chr}{$pos}{"ALT"};
							$long_alt = $tum_indel_fields{$chr}{$pos}{"ALT"};
						}else{
							$long_ref = $tum_indel_fields{$chr}{$pos}{"REF"};
							$short_alt = $tum_indel_fields{$chr}{$pos}{"ALT"};
							$long_alt = $norm_indel_fields{$chr}{$pos}{"ALT"};							
						}
						
						my $tmp_alt = $long_ref;
						$long_alt =~ s/^$short_alt//;
						$tmp_alt =~ s/^$short_alt//;
						my $new_alt = $short_alt.$long_alt.$tmp_alt;
						
						$norm_indel_fields{$chr}{$pos}{"REF"} = $long_ref;
						$norm_indel_fields{$chr}{$pos}{"ALT"} = "$short_alt,$new_alt";
						$tum_indel_fields{$chr}{$pos}{"REF"} = $long_ref;
						$tum_indel_fields{$chr}{$pos}{"ALT"} = "$short_alt,$new_alt";
						
						
					}elsif($norm_indel_fields{$chr}{$pos}{"REF"} ne $tum_indel_fields{$chr}{$pos}{"REF"}){
						my $com_alt = $norm_indel_fields{$chr}{$pos}{"ALT"};
						my @n_ref = split(//, $norm_indel_fields{$chr}{$pos}{"REF"});
						my @t_ref = split(//, $tum_indel_fields{$chr}{$pos}{"REF"});
						my $new_alt;
						my $new_ref;
						
						shift(@n_ref);
						shift(@t_ref);
						
						my $new_n_ref = join("", @n_ref);
						my $new_t_ref = join("", @t_ref);
						
						if(length($new_n_ref) > length($new_t_ref)){
							$new_n_ref =~ s/^$new_t_ref//;
							$new_alt = $com_alt.$new_n_ref;
							$new_ref = $norm_indel_fields{$chr}{$pos}{"REF"};
						}else{
							$new_t_ref =~ s/^$new_n_ref//;
							$new_alt = $com_alt.$new_t_ref;
							$new_ref = $tum_indel_fields{$chr}{$pos}{"REF"};
						}
						
						$norm_indel_fields{$chr}{$pos}{"REF"} = $new_ref;
						$norm_indel_fields{$chr}{$pos}{"ALT"} = "$new_alt,$com_alt";
						$tum_indel_fields{$chr}{$pos}{"REF"} = $new_ref;
						$tum_indel_fields{$chr}{$pos}{"ALT"} = "$new_alt,$com_alt";
						
					}else{
						$norm_indel_fields{$chr}{$pos}{"ALT"} = $norm_indel_fields{$chr}{$pos}{"ALT"}.",".$tum_indel_fields{$chr}{$pos}{"ALT"};
						$tum_indel_fields{$chr}{$pos}{"ALT"} = $norm_indel_fields{$chr}{$pos}{"ALT"};
					}
					
					
					if($norm_indel_fields{$chr}{$pos}{"VT"} ne $tum_indel_fields{$chr}{$pos}{"ALT"}){
						$norm_indel_fields{$chr}{$pos}{"VT"}=$norm_indel_fields{$chr}{$pos}{"VT"}.",".$tum_indel_fields{$chr}{$pos}{"VT"};
						$tum_indel_fields{$chr}{$pos}{"VT"}=$norm_indel_fields{$chr}{$pos}{"VT"};
					}
					if($norm_indel_fields{$chr}{$pos}{"RSID"} ne $tum_indel_fields{$chr}{$pos}{"RSID"}){
						$norm_indel_fields{$chr}{$pos}{"RSID"}=$norm_indel_fields{$chr}{$pos}{"RSID"}.",".$tum_indel_fields{$chr}{$pos}{"RSID"};
						$tum_indel_fields{$chr}{$pos}{"RSID"}=$norm_indel_fields{$chr}{$pos}{"RSID"};
					}

					if($tum_indel_fields{$chr}{$pos}{"GT"} eq "1/1"){
						$tum_indel_fields{$chr}{$pos}{"GT"} = "2/2";
					}else{
						$tum_indel_fields{$chr}{$pos}{"GT"} = "0/2";
					}
					add_filter($chr, $pos, "MIF", %norm_indel_fields);
					add_filter($chr, $pos, "MIF", %tum_indel_fields);
				}
				
				if($norm_indel_fields{$chr}{$pos}{"GT"} ne "0/0" or $tum_indel_fields{$chr}{$pos}{"GT"} ne "0/0"){
					if($norm_indel_fields{$chr}{$pos}{"GT"} ne "0/0"){
	#GERMLINE or LOH
						print "$chr\t$pos\t".$norm_indel_fields{$chr}{$pos}{"RSID"}."\t".$norm_indel_fields{$chr}{$pos}{"REF"}."\t".$norm_indel_fields{$chr}{$pos}{"ALT"}."\t".$norm_indel_fields{$chr}{$pos}{"P"}."\t".$norm_indel_fields{$chr}{$pos}{"FILTER"}."\tDP=".($norm_indel_fields{$chr}{$pos}{"DP"}+$tum_indel_fields{$chr}{$pos}{"DP"}).";BQ=.;MQ=.;ER=$e;VT=".$norm_indel_fields{$chr}{$pos}{"VT"}.";FEP=".$norm_indel_fields{$chr}{$pos}{"FEP"}."\tGT:AF:DP4:AD:BQ2:NR:DS:RGB:SS\t".$norm_indel_fields{$chr}{$pos}{"GT"}.":".$norm_indel_fields{$chr}{$pos}{"AAF"}.":.:".$norm_indel_fields{$chr}{$pos}{"AD"}.":.:".$normal_fields{$chr}{$pos}{"NR"}.":".$norm_indel_fields{$chr}{$pos}{"DS"}.":.:1\t".$tum_indel_fields{$chr}{$pos}{"GT"}.":".$tum_indel_fields{$chr}{$pos}{"AAF"}.":.:".$tum_indel_fields{$chr}{$pos}{"AD"}.":.:".$tumor_fields{$chr}{$pos}{"NR"}.":".$tum_indel_fields{$chr}{$pos}{"DS"}.":.:1\n"

					}else{
	#SOMATIC
						print "$chr\t$pos\t".$tum_indel_fields{$chr}{$pos}{"RSID"}."\t".$tum_indel_fields{$chr}{$pos}{"REF"}."\t".$tum_indel_fields{$chr}{$pos}{"ALT"}."\t".$tum_indel_fields{$chr}{$pos}{"P"}."\t".$tum_indel_fields{$chr}{$pos}{"FILTER"}."\tDP=".($norm_indel_fields{$chr}{$pos}{"DP"}+$tum_indel_fields{$chr}{$pos}{"DP"}).";BQ=.;MQ=.;ER=$e;VT=".$tum_indel_fields{$chr}{$pos}{"VT"}.";FEP=".$tum_indel_fields{$chr}{$pos}{"FEP"}."\tGT:AF:DP4:AD:BQ2:NR:DS:RGB:SS\t".$norm_indel_fields{$chr}{$pos}{"GT"}.":".$norm_indel_fields{$chr}{$pos}{"AAF"}.":.:".$norm_indel_fields{$chr}{$pos}{"AD"}.":.:".$normal_fields{$chr}{$pos}{"NR"}.":".$norm_indel_fields{$chr}{$pos}{"DS"}.":.:0\t".$tum_indel_fields{$chr}{$pos}{"GT"}.":".$tum_indel_fields{$chr}{$pos}{"AAF"}.":.:".$tum_indel_fields{$chr}{$pos}{"AD"}.":.:".$tumor_fields{$chr}{$pos}{"NR"}.":".$tum_indel_fields{$chr}{$pos}{"DS"}.":.:2\n"
						
					}
				}
			}
		}
		
		
		
####
# SNPS
####

## Add coverage filter info
		if($normal_fields{$chr}{$pos}{"DP"} < $min_norm_cov or $tumor_fields{$chr}{$pos}{"DP"} < $min_tum_cov){
			add_filter($chr, $pos, "MCF", %normal_fields);
			add_filter($chr, $pos, "MCF", %tumor_fields);
		}

		if($normal_fields{$chr}{$pos}{"GT"} eq "0/0" and $tumor_fields{$chr}{$pos}{"GT"} eq "0/0"){
			next POS_LOOP;
#potential REF
		}
		if($normal_fields{$chr}{$pos}{"GT"} ne "0/0"){
#potential GERMLINE
			if($tumor_fields{$chr}{$pos}{"GT"} ne $normal_fields{$chr}{$pos}{"GT"}){
				print "$chr\t$pos\t".$normal_fields{$chr}{$pos}{"RSID"}."\t".$normal_fields{$chr}{$pos}{"REF"}."\t".$normal_fields{$chr}{$pos}{"ALT"}."\t".$normal_fields{$chr}{$pos}{"P"}."\t".$normal_fields{$chr}{$pos}{"FILTER"}."\tDP=".($normal_fields{$chr}{$pos}{"DP"}+$tumor_fields{$chr}{$pos}{"DP"}).";BQ=".(($normal_fields{$chr}{$pos}{"BQ"}+$tumor_fields{$chr}{$pos}{"BQ"})/2).";MQ=".(($normal_fields{$chr}{$pos}{"MQ"}+$tumor_fields{$chr}{$pos}{"MQ"})/2).";ER=".(($normal_fields{$chr}{$pos}{"ER"}+$tumor_fields{$chr}{$pos}{"ER"})/2).";VT=SNP;FEP=".$normal_fields{$chr}{$pos}{"FEP"}.";GQ=".$normal_fields{$chr}{$pos}{"GQ"}.";PL=".$normal_fields{$chr}{$pos}{"PL"}."\tGT:AF:DP4:AD:BQ2:MQ2:NR:ST:DS:RGB:SS\t".$normal_fields{$chr}{$pos}{"GT"}.":".$normal_fields{$chr}{$pos}{"UAA"}.":".$normal_fields{$chr}{$pos}{"DP4"}.":".$normal_fields{$chr}{$pos}{"AD"}.":".$normal_fields{$chr}{$pos}{"BQ2"}.":".$normal_fields{$chr}{$pos}{"MQ2"}.":".$normal_fields{$chr}{$pos}{"NR"}.":".$normal_fields{$chr}{$pos}{"ST"}.":".$normal_fields{$chr}{$pos}{"DS"}.":".$normal_fields{$chr}{$pos}{"RGB"}.":1\t".$tumor_fields{$chr}{$pos}{"GT"}.":".$tumor_fields{$chr}{$pos}{"UAA"}.":".$tumor_fields{$chr}{$pos}{"DP4"}.":".$tumor_fields{$chr}{$pos}{"AD"}.":".$tumor_fields{$chr}{$pos}{"BQ2"}.":".$tumor_fields{$chr}{$pos}{"MQ2"}.":".$tumor_fields{$chr}{$pos}{"NR"}.":".$tumor_fields{$chr}{$pos}{"ST"}.":".$tumor_fields{$chr}{$pos}{"DS"}.":".$tumor_fields{$chr}{$pos}{"RGB"}.":3\n";			
			}else{
				print "$chr\t$pos\t".$normal_fields{$chr}{$pos}{"RSID"}."\t".$normal_fields{$chr}{$pos}{"REF"}."\t".$normal_fields{$chr}{$pos}{"ALT"}."\t".$normal_fields{$chr}{$pos}{"P"}."\t".$normal_fields{$chr}{$pos}{"FILTER"}."\tDP=".($normal_fields{$chr}{$pos}{"DP"}+$tumor_fields{$chr}{$pos}{"DP"}).";BQ=".(($normal_fields{$chr}{$pos}{"BQ"}+$tumor_fields{$chr}{$pos}{"BQ"})/2).";MQ=".(($normal_fields{$chr}{$pos}{"MQ"}+$tumor_fields{$chr}{$pos}{"MQ"})/2).";ER=".(($normal_fields{$chr}{$pos}{"ER"}+$tumor_fields{$chr}{$pos}{"ER"})/2).";VT=SNP;FEP=".$tumor_fields{$chr}{$pos}{"FEP"}.";GQ=".$normal_fields{$chr}{$pos}{"GQ"}.";PL=".$normal_fields{$chr}{$pos}{"PL"}."\tGT:AF:DP4:AD:BQ2:MQ2:NR:ST:DS:RGB:SS\t".$normal_fields{$chr}{$pos}{"GT"}.":".$normal_fields{$chr}{$pos}{"UAA"}.":".$normal_fields{$chr}{$pos}{"DP4"}.":".$normal_fields{$chr}{$pos}{"AD"}.":".$normal_fields{$chr}{$pos}{"BQ2"}.":".$normal_fields{$chr}{$pos}{"MQ2"}.":".$normal_fields{$chr}{$pos}{"NR"}.":".$normal_fields{$chr}{$pos}{"ST"}.":".$normal_fields{$chr}{$pos}{"DS"}.":".$normal_fields{$chr}{$pos}{"RGB"}.":1\t".$tumor_fields{$chr}{$pos}{"GT"}.":".$tumor_fields{$chr}{$pos}{"UAA"}.":".$tumor_fields{$chr}{$pos}{"DP4"}.":".$tumor_fields{$chr}{$pos}{"AD"}.":".$tumor_fields{$chr}{$pos}{"BQ2"}.":".$tumor_fields{$chr}{$pos}{"MQ2"}.":".$tumor_fields{$chr}{$pos}{"NR"}.":".$tumor_fields{$chr}{$pos}{"ST"}.":".$tumor_fields{$chr}{$pos}{"DS"}.":".$tumor_fields{$chr}{$pos}{"RGB"}.":1\n";
			}
			next POS_LOOP;
		}

#potential Somatic		
		print "$chr\t$pos\t".$tumor_fields{$chr}{$pos}{"RSID"}."\t".$tumor_fields{$chr}{$pos}{"REF"}."\t".$tumor_fields{$chr}{$pos}{"ALT"}."\t".$tumor_fields{$chr}{$pos}{"P"}."\t".$tumor_fields{$chr}{$pos}{"FILTER"}."\tDP=".($normal_fields{$chr}{$pos}{"DP"}+$tumor_fields{$chr}{$pos}{"DP"}).";BQ=".(($normal_fields{$chr}{$pos}{"BQ"}+$tumor_fields{$chr}{$pos}{"BQ"})/2).";MQ=".(($normal_fields{$chr}{$pos}{"MQ"}+$tumor_fields{$chr}{$pos}{"MQ"})/2).";ER=".(($normal_fields{$chr}{$pos}{"ER"}+$tumor_fields{$chr}{$pos}{"ER"})/2).";VT=SNP;FEP=".$tumor_fields{$chr}{$pos}{"FEP"}."\tGT:AF:DP4:AD:BQ2:MQ2:NR:ST:DS:RGB:SS\t"."0/0:".$normal_fields{$chr}{$pos}{"UAA"}.":".$normal_fields{$chr}{$pos}{"DP4"}.":".$normal_fields{$chr}{$pos}{"AD"}.":".$normal_fields{$chr}{$pos}{"BQ2"}.":".$normal_fields{$chr}{$pos}{"MQ2"}.":".$normal_fields{$chr}{$pos}{"NR"}.":".$normal_fields{$chr}{$pos}{"ST"}.":".$normal_fields{$chr}{$pos}{"DS"}.":".$normal_fields{$chr}{$pos}{"RGB"}.":0\t".$tumor_fields{$chr}{$pos}{"GT"}.":".$tumor_fields{$chr}{$pos}{"UAA"}.":".$tumor_fields{$chr}{$pos}{"DP4"}.":".$tumor_fields{$chr}{$pos}{"AD"}.":".$tumor_fields{$chr}{$pos}{"BQ2"}.":".$tumor_fields{$chr}{$pos}{"MQ2"}.":".$tumor_fields{$chr}{$pos}{"NR"}.":".$tumor_fields{$chr}{$pos}{"ST"}.":".$tumor_fields{$chr}{$pos}{"DS"}.":".$tumor_fields{$chr}{$pos}{"RGB"}.":2\n";
		
		
	}
}










sub get_sample_vcf($\%\%@){
	my ($samp, $hash, $ind_hash, $chr, $pos, $ref, $cov, $indels, @rgs) = @_;
	$ref =~ tr/[a-z]/[A-Z]/;
	$$hash{$chr}{$pos}{"REF"} = $ref;
	$$hash{$chr}{$pos}{"ALT"} = ".";
	my $ref_cov=0;
	my $alt_cov=0;
	my $strand;
	my $read;
	my %errors;
	my %err_cov;
	my $alt = ".";
	my @tot_bq;
	my @tot_mq;
	my %alt_bq;
	my @ref_bq;
	my %alt_mq;
	my @ref_mq;
	my %strands;
	my %dp4;
	my $indel_cov = $cov;
	my %unfiltSNPs = ('A',0,'T',0,'C',0,'G',0);
	my $tot_cov = 0;
	my %filtdp4;
	my %rg_bias_info;
	my %loglikelihood;

	my %SNPs = ('A',0,'T',0,'C',0,'G',0);
	my $total_rgs = $#rgs+1;
	foreach my $rg_info (@rgs){
		my ($rg, $bases_info, $bq_info, $mq_info, $bp_info) = split(/\s/, $rg_info);
		
		my @bases = split(//, $bases_info);
		my @bq = split(//, $bq_info);
		my @bp = split(/,/, $bp_info);
		my @mq = split(//, $mq_info);
		$read = "R1" if($rg =~ m/R1/);
		$read = "R2" if($rg =~ m/R2/);
		
		for(my $i = 0; $i <= $#bases; $i++){
			my $base_qual = ord($bq[$i])-33;
			my $map_qual = ord($mq[$i])-33;
			
			if(exists($skip{$chr}{$pos}{uc($bases[$i])})){
				$cov--;
			}elsif($bases[$i] eq "*"){
				$cov--;
			}else{
				$tot_cov++;
				
				push(@tot_bq, $base_qual);
				push(@tot_mq, $map_qual);
				
				if($bases[$i] eq "." or $bases[$i] =~ m/[ATCG]/){
					$strand = "FOR";
				}else{
					$strand = "REV";
				}
				$strands{$strand} = 1;
				
				if(uc($bases[$i]) =~ m/[ATCG]/){
					push(@{$alt_bq{uc($bases[$i])}}, $base_qual);
					push(@{$alt_mq{uc($bases[$i])}}, $map_qual);
					$unfiltSNPs{uc($bases[$i])}++;
					$dp4{uc($bases[$i])}{$strand}++;
					$rg_bias_info{$rg}{uc($bases[$i])}++;
					
					if($base_qual >= $BQ and $map_qual >= $MQ){
						$filtdp4{uc($bases[$i])}{$strand}++;
						$SNPs{uc($bases[$i])}++;
						$err_cov{$read}{$bp[$i]}++;
						if(exists($error_rate{$read}{$bp[$i]}{get_sub_type($ref, uc($bases[$i]))})){
							$errors{uc($bases[$i])}{$read}{$bp[$i]} = $error_rate{$read}{$bp[$i]}{get_sub_type($ref, uc($bases[$i]))};
						}else{
							$errors{uc($bases[$i])}{$read}{$bp[$i]} = $e;
						}
					}else{
						$cov--;
					}
				}else{
					push(@ref_bq, $base_qual);
					push(@ref_mq, $map_qual);
					$ref_cov++;
					$dp4{"ref"}{$strand}++;
					$rg_bias_info{$rg}{"ref"}++;

					if($base_qual >= $BQ and $map_qual >= $MQ){
						$filtdp4{"ref"}{$strand}++;
					}else{
						$cov--;
					}
				}
			}
		}
	}
	
	
	my @sts = keys %strands;
	$$hash{$chr}{$pos}{"ST"} = $#sts+1;
	$$hash{$chr}{$pos}{"NR"} = $total_rgs; 
	$$hash{$chr}{$pos}{"RSID"} = ".";
	$$hash{$chr}{$pos}{"P"} = 0.000;
	$$hash{$chr}{$pos}{"FILTER"} = "PASS";
	$$hash{$chr}{$pos}{"DP"} = $tot_cov;
	$$hash{$chr}{$pos}{"FDP"} = $cov;
	$$hash{$chr}{$pos}{"BQ"} = 0;
	$$hash{$chr}{$pos}{"MQ"} = 0;
	$$hash{$chr}{$pos}{"ER"} = 0;
	$$hash{$chr}{$pos}{"DS"} = 0;
	$$hash{$chr}{$pos}{"AD"} = "$ref_cov,0";
	$$hash{$chr}{$pos}{"DP4"} = "0,0,0,0";
	$$hash{$chr}{$pos}{"BQ2"} = "0,0";
	$$hash{$chr}{$pos}{"MQ2"} = "0,0";
	$$hash{$chr}{$pos}{"GT"} = "0/0";
	$$hash{$chr}{$pos}{"UAA"} = 0;
	$$hash{$chr}{$pos}{"FILTDP4"} = "0,0,0,0";
	$$hash{$chr}{$pos}{"RGB"} = "$chr\t$pos\t0_0";
	$$hash{$chr}{$pos}{"GQ"} = 0;
	$$hash{$chr}{$pos}{"PL"} = "0,0,0";
	
	return () if($cov <= 0);	
	
#Determine first and second greatest covered bases
	my @sorted = sort {$b <=> $a} values %unfiltSNPs;
	my $sec_cov = $sorted[1];
	
	my @filtered_sorted = sort {$b <=> $a} values %SNPs;
	my $filtered_sec_cov = $filtered_sorted[1];
	my $filtered_alt_cov = $filtered_sorted[0];
	my $filtered_tot_nonRef_cov = $filtered_sorted[0] + $filtered_sorted[1] + $filtered_sorted[2];
	
	my $tot_nonRef_cov = 0;
	foreach my $b (keys %unfiltSNPs){
		$tot_nonRef_cov += $unfiltSNPs{$b};
		if($unfiltSNPs{$b}>$alt_cov){
			$alt_cov = $unfiltSNPs{$b} if($SNPs{$b} != 0);
			$alt = $b if($SNPs{$b} != 0);
		}
	}

#Get genotype for normal samples	
	if($alt ne "." and ($alt_cov/$tot_cov)*100 >= $gaaf){
		my $hom_ref_lls = 0;
		my $het_lls = 0;
		my $hom_alt_lls = 0;
		
		foreach my $rg_info (@rgs){
			my ($rg, $bases_info, $bq_info, $mq_info, $bp_info) = split(/\s/, $rg_info);
			$bases_info =~ tr/[a-z]/[A-Z]/;
			my @bases = split(//, $bases_info);
			my @bq = split(//, $bq_info);
			my @mq = split(//, $mq_info);
		
			for(my $i = 0; $i <= $#bases; $i++){
				my $base_qual = ord($bq[$i])-33;
				my $map_qual = ord($mq[$i])-33;
				
				if(!(exists($skip{$chr}{$pos}{uc($bases[$i])})) and uc($bases[$i]) ne "*"){
					if($base_qual >= $BQ and $map_qual >= $MQ){
						my $b = uc($bases[$i]);
						if($b eq "." or $b eq ","){
							$b = $ref;
						}
						$hom_ref_lls += log10(get_prob_b_given_g($ref, $ref, $b, $phred_to_pval{$base_qual}));
						$het_lls += log10(get_prob_b_given_g($ref, $alt, $b, $phred_to_pval{$base_qual}));
						$hom_alt_lls += log10(get_prob_b_given_g($alt, $alt, $b, $phred_to_pval{$base_qual}));
					}
				}
			}
		}
		
		$hom_ref_lls+=log10($hom_ref);
		$het_lls+=log10($het);
		$hom_alt_lls+=log10($hom_alt);

		$hom_ref_lls = -1*($hom_ref_lls);
		$het_lls = -1*($het_lls);
		$hom_alt_lls = -1*($hom_alt_lls);
					
		my $max_gt = "0/0";
		my $gt_qual = 0;
		my $pl = "";
		my $sec_lls = 0;

		if($hom_ref_lls >= $het_lls or $hom_ref_lls >= $hom_alt_lls){
			if($het_lls <= $hom_alt_lls){
				$sec_lls = $hom_ref_lls if($hom_ref_lls <= $hom_alt_lls);
				$sec_lls = $hom_alt_lls if($hom_ref_lls > $hom_alt_lls);

				$max_gt = "0/1";
				$gt_qual = int($sec_lls - $het_lls + 0.5);
				$pl = int($hom_ref_lls+0.5).",".int($het_lls+0.5).",".int($hom_alt_lls+0.5);
			}else{
				$sec_lls = $hom_ref_lls if($hom_ref_lls <= $het_lls);
				$sec_lls = $het_lls if($hom_ref_lls > $het_lls);
				
				$max_gt = "1/1";
				$gt_qual = int($sec_lls - $hom_alt_lls + 0.5);
				$pl = int($hom_ref_lls+0.5).",".int($het_lls+0.5).",".int($hom_alt_lls+0.5);				
			}
			$$hash{$chr}{$pos}{"GQ"} = $gt_qual;
			$$hash{$chr}{$pos}{"GT"} = $max_gt;
			$$hash{$chr}{$pos}{"PL"} = $pl;
		}
	}

	
	$$hash{$chr}{$pos}{"RGB"} = "$chr\t$pos\t";
	foreach my $t_rg (sort keys %rg_bias_info){
		if(exists($rg_bias_info{$t_rg}{"ref"})){
			$$hash{$chr}{$pos}{"RGB"} .= ($rg_bias_info{$t_rg}{"ref"}."_");
		}else{
			$$hash{$chr}{$pos}{"RGB"} .= "0_";
		}
		if(exists($rg_bias_info{$t_rg}{$alt})){
			$$hash{$chr}{$pos}{"RGB"} .= "$rg_bias_info{$t_rg}{$alt}_"; 
		}else{
			$$hash{$chr}{$pos}{"RGB"} .= "0_";
		}
	}

	$$hash{$chr}{$pos}{"UAA"} = sprintf("%.2f", (($unfiltSNPs{$alt}/$tot_cov) * 100)) if($alt ne ".");

	if($$hash{$chr}{$pos}{"UAA"} >= $saaf and $samp eq "tum"){
		$$hash{$chr}{$pos}{"GT"} = "0/1" if($$hash{$chr}{$pos}{"GT"} eq "0/0");
	}

	$$hash{$chr}{$pos}{"DS"} = sprintf("%.2f", (($filtered_alt_cov - $filtered_sec_cov)/($filtered_tot_nonRef_cov))) if($filtered_tot_nonRef_cov != 0 and $alt ne ".");
	$$hash{$chr}{$pos}{"DS"} = 1 if($$hash{$chr}{$pos}{"DS"} == 1);
	$$hash{$chr}{$pos}{"ALT"} = $alt;


	$$hash{$chr}{$pos}{"BQ"} = rms(@tot_bq);
	$$hash{$chr}{$pos}{"MQ"} = rms(@tot_mq);
	$$hash{$chr}{$pos}{"ER"} = sprintf("%.4e", mean(%errors, %err_cov, %ref_error, $alt));
	$$hash{$chr}{$pos}{"RSID"} = $dbsnp{$chr}{$pos}{$ref}{$alt} if(exists($dbsnp{$chr}{$pos}{$ref}{$alt}));
	my $alt_rms_val = rms(@{$alt_bq{$alt}});
	my $ref_rms_val = rms(@ref_bq);
	my $alt_mq_rms_val = rms(@{$alt_mq{$alt}});
	my $ref_mq_rms_val = rms(@ref_mq);
	$$hash{$chr}{$pos}{"BQ2"} = "$ref_rms_val,$alt_rms_val";
	$$hash{$chr}{$pos}{"MQ2"} = "$ref_mq_rms_val,$alt_mq_rms_val";
	
	my $ref_for = 0;
	my $ref_rev = 0;
	my $alt_for = 0;
	my $alt_rev = 0;
	$ref_for = $dp4{"ref"}{"FOR"} if(exists($dp4{"ref"}{"FOR"}));
	$ref_rev = $dp4{"ref"}{"REV"} if(exists($dp4{"ref"}{"REV"}));
	$alt_for = $dp4{$alt}{"FOR"} if(exists($dp4{$alt}{"FOR"}));
	$alt_rev = $dp4{$alt}{"REV"} if(exists($dp4{$alt}{"REV"}));
	$$hash{$chr}{$pos}{"DP4"} = "$ref_for,$ref_rev,$alt_for,$alt_rev";
	$$hash{$chr}{$pos}{"AD"} = ($ref_for+$ref_rev).",".($alt_for+$alt_rev);

	my $unrfor = 0;
	my $unrrev = 0;
	my $unafor = 0;
	my $unarev = 0;
	$unrfor = $filtdp4{"ref"}{"FOR"} if(exists($filtdp4{"ref"}{"FOR"}));
	$unrrev = $filtdp4{"ref"}{"REV"} if(exists($filtdp4{"ref"}{"REV"}));
	$unafor = $filtdp4{$alt}{"FOR"} if(exists($filtdp4{$alt}{"FOR"}));
	$unarev = $filtdp4{$alt}{"REV"} if(exists($filtdp4{$alt}{"REV"}));
	$$hash{$chr}{$pos}{"FILTDP4"} = "$unrfor,$unrrev,$unafor,$unarev";
	
	
	



	$$ind_hash{$chr}{$pos}{"FILTER"} = "NA";
	$$ind_hash{$chr}{$pos}{"VT"} = "NA";
	$$ind_hash{$chr}{$pos}{"IND"} = "NA";
	$$ind_hash{$chr}{$pos}{"COV"} = 0;
	$$ind_hash{$chr}{$pos}{"DS"} = ".";
	$$ind_hash{$chr}{$pos}{"ST"} = ".";
	$$ind_hash{$chr}{$pos}{"DP"} = $indel_cov;
	$$ind_hash{$chr}{$pos}{"AAF"} = 0;
	$$ind_hash{$chr}{$pos}{"RC"} = $indel_cov;
	$$ind_hash{$chr}{$pos}{"GT"} = "0/0";
	$$ind_hash{$chr}{$pos}{"REF"} = ".";
	$$ind_hash{$chr}{$pos}{"ALT"} = ".";
	$$ind_hash{$chr}{$pos}{"RSID"} = ".";
	$$ind_hash{$chr}{$pos}{"AD"} = "$indel_cov,0";
	if($indels =~ m/\w/){

		my (@indels_info) = split(/\s/, $indels);
		my %tmp_indel;
		foreach (@indels_info){
			my ($s1, $strnd, $loc, @inds) = split(/_/, $_);
			foreach my $ind_in (@inds){
				my ($ind, $count) = split(/:/, $ind_in);
				$tmp_indel{$ind}{"count"} += $count;
				$tmp_indel{$ind}{"strand"}{$strnd} += $count;
				$tmp_indel{$ind}{"rgs"}{"$strnd\_$loc"} = 1;
			}
		}
		
		my $max_ind;
		my $max_count = 0;
		my @counts;
		my $tot_ind_count = 0;
		my $ind_num_count = 0;
		foreach (sort keys %tmp_indel){
			$ind_num_count++ if($tmp_indel{$_}{"count"} >= 5);
			push(@counts, $tmp_indel{$_}{"count"});
			if($tmp_indel{$_}{"count"} > $max_count){
				$max_count = $tmp_indel{$_}{"count"};
				$max_ind = $_;
			}
			$tot_ind_count += $tmp_indel{$_}{"count"};
		}

		@counts = sort {$b <=> $a} @counts;
		my $sec_ind_c = 0;
		$sec_ind_c = $counts[1] if($#counts >= 1);
		
		my $ind_rg_count = 0;
		foreach (keys %{$tmp_indel{$max_ind}{"rgs"}}){ $ind_rg_count++; }

		$$ind_hash{$chr}{$pos}{"FILTER"} = "PASS";
		$$ind_hash{$chr}{$pos}{"VT"} = "DEL" if($max_ind =~ m/-/);
		$$ind_hash{$chr}{$pos}{"VT"} = "INS" if($max_ind =~ m/\+/);

		$max_ind =~ s/\+|-//;
		$max_ind =~ s/\d+//;
		$max_ind =~ tr/[a-z]/[A-Z]/;
		$$ind_hash{$chr}{$pos}{"IND"} = $max_ind;
		$$ind_hash{$chr}{$pos}{"COV"} = $max_count;
		$$ind_hash{$chr}{$pos}{"DS"} = sprintf("%.2f", (($max_count - $sec_ind_c)/($tot_ind_count)));
		$$ind_hash{$chr}{$pos}{"DP"} = $indel_cov;
		$$ind_hash{$chr}{$pos}{"AAF"} = sprintf("%.2f",($max_count/$indel_cov * 100));
		$$ind_hash{$chr}{$pos}{"RC"} = $indel_cov-$max_count;
		if($samp eq "norm"){
			$$ind_hash{$chr}{$pos}{"GT"} = "0/1" if($$ind_hash{$chr}{$pos}{"AAF"}>=$gaaf);
			$$ind_hash{$chr}{$pos}{"GT"} = "1/1" if($$ind_hash{$chr}{$pos}{"AAF"}>=75);		
		}else{
			$$ind_hash{$chr}{$pos}{"GT"} = "0/1" if($$ind_hash{$chr}{$pos}{"AAF"}>=$saaf);
			$$ind_hash{$chr}{$pos}{"GT"} = "1/1" if($$ind_hash{$chr}{$pos}{"AAF"}>=75);			
		}

		my ($ind_ref, $ind_alt, $ind_filt) = get_indel_information($chr, $pos, $ref, $max_ind, $$ind_hash{$chr}{$pos}{"VT"}, $$ind_hash{$chr}{$pos}{"FILTER"}, "$pdir/intermediate/$normal\_tmp_indel_fasta_file.fa");
		
		$$ind_hash{$chr}{$pos}{"REF"} = $ind_ref;
		$$ind_hash{$chr}{$pos}{"ALT"} = $ind_alt;
		$$ind_hash{$chr}{$pos}{"FILTER"} = $ind_filt;
		$$ind_hash{$chr}{$pos}{"RSID"} = $dbsnp{$chr}{$pos}{$ind_ref}{$ind_alt} if(exists($dbsnp{$chr}{$pos}{$ind_ref}{$ind_alt}));
		$$ind_hash{$chr}{$pos}{"AD"} = $$ind_hash{$chr}{$pos}{"RC"}.",".$$ind_hash{$chr}{$pos}{"COV"};
		
	}
	
}



sub get_sub_type($$){
	my ($ref, $alt) = @_;
	my $type;
	
	if($ref eq "A"){
		if($alt eq "T"){
			$type = "ATTA";
		}elsif($alt eq "G"){
			$type = "ATGC";
		}else{
			$type = "ATCG";
		}
	}elsif($ref eq "T"){
		if($alt eq "A"){
			$type = "ATTA";
		}elsif($alt eq "G"){
			$type = "ATCG";
		}else{
			$type = "ATGC";
		}		
	}elsif($ref eq "C"){
		if($alt eq "A"){
			$type = "CGAT";
		}elsif($alt eq "G"){
			$type = "CGGC";
		}else{
			$type = "CGTA";
		}
	}else{
		if($alt eq "A"){
			$type = "CGTA";
		}elsif($alt eq "C"){
			$type = "CGGC";
		}else{
			$type = "CGAT";
		}		
	}
	
	return $type;
}



sub mean(\%\%\%$){
	my ($err_hash, $cov_hash, $ref_err_hash, $alt) = @_;
	
	return($e) if($alt eq ".");
	
	my $total = 0;
	my $errors = 0;
	foreach my $st (keys %{$cov_hash}){
		foreach my $pos (keys %{$$cov_hash{$st}}){
			if(exists($$err_hash{$alt}{$st}{$pos})){
				$errors += $$err_hash{$alt}{$st}{$pos}*$$cov_hash{$st}{$pos};
			}else{
				if(exists($$ref_err_hash{$st}{$pos})){
					$errors += $$ref_err_hash{$st}{$pos}*$$cov_hash{$st}{$pos};
				}else{
					$errors += $e*$$cov_hash{$st}{$pos};
				}
			}
			$total += $$cov_hash{$st}{$pos};
		}
	}
	
	my $mean = $errors/$total;
	return($mean);	
}


sub get_error_rates($\%){
	my ($file, $hash) = @_;
	open(ERROR, $file);
	while(my $line = <ERROR>){
		chomp $line;
		my ($read, $pos, $type, $error, @other) = split(/\t/, $line);
		$$hash{$read}{$pos}{$type} = $error;
	}
	close(ERROR);
}


sub avg_error(\%\%){
	my ($err, $refe) = @_;
	foreach my $st (keys %$err){
		foreach my $pos (keys %{$$err{$st}}){
			my $tot = 0;
			my $errors = 0;
			foreach my $t (keys %{$$err{$st}{$pos}}){
				$tot++;
				$errors += $$err{$st}{$pos}{$t};
			}
			my $mean = $errors/$tot;
			$$refe{$st}{$pos} = $mean;
		}
	}
}

sub get_black_list($\%){
	my ($file, $hash) = @_;
	open(FILE, $file);
	while(my $line = <FILE>){
		chomp $line;
		my ($chr, $pos, $ref, $alt) = split(/\t/, $line);
		next if($ref eq "-");
		next if($alt eq "-");
		$$hash{$chr}{$pos}{$alt} = 1;
	}
	close(FILE);
}

sub get_dbsnp($\%\%){
	my ($file, $hash, $intervals) = @_;
	open(DBSNP, $file);
	while(my $line = <DBSNP>){
		chomp $line;
		next if($line =~ m/^#/);
		my ($chr, $pos, $rsid, $ref, $alt, @other) = split(/\t/, $line);
		
		INTLOOP: foreach my $start (sort {$a <=> $b} keys %{$$intervals{$chr}}){
			my $stop = $$intervals{$chr}{$start};
			if($pos >= $start and $pos <= $stop){
				if($alt =~ m/,/){
					my @alts = split(/,/, $alt);
					foreach my $a (@alts){
						$$hash{$chr}{$pos}{$ref}{$a} = $rsid;
					}
				}else{
					$$hash{$chr}{$pos}{$ref}{$alt} = $rsid;
				}
			}
			last INTLOOP if($pos < $start);
		}
		
	}
	close(DBSNP);
}


sub rms(@){
	my (@vals) = @_;
	
	if(not @vals){
		return(0);
	}
	my $tot=0;
	my $sq=0;
	foreach my $v (@vals){
		$tot++;
		$sq += ($v**2);
	}
	my $rms = ($sq/$tot)**0.5;
	return (int($rms+0.5));
}


sub set_chrom_order(\%){
	my ($chrom_order) = @_;
	$$chrom_order{"chr1"}=1;
	$$chrom_order{"chr2"}=2;
	$$chrom_order{"chr3"}=3;
	$$chrom_order{"chr4"}=4;
	$$chrom_order{"chr5"}=5;
	$$chrom_order{"chr6"}=6;
	$$chrom_order{"chr7"}=7;
	$$chrom_order{"chr8"}=8;
	$$chrom_order{"chr9"}=9;
	$$chrom_order{"chr10"}=10;
	$$chrom_order{"chr11"}=11;
	$$chrom_order{"chr12"}=12;
	$$chrom_order{"chr13"}=13;
	$$chrom_order{"chr14"}=14;
	$$chrom_order{"chr15"}=15;
	$$chrom_order{"chr16"}=16;
	$$chrom_order{"chr17"}=17;
	$$chrom_order{"chr18"}=18;
	$$chrom_order{"chr19"}=19;
	$$chrom_order{"chr20"}=20;
	$$chrom_order{"chr21"}=21;
	$$chrom_order{"chr22"}=22;
	$$chrom_order{"chrX"}=23;
	$$chrom_order{"chrY"}=24;
}



sub find_indel_repeat($$){
	my ($up_fa, $down_fa) = @_;
	my @up = split(//, $up_fa);
	my @down = split(//, $down_fa);
	my $repeat = 0;
	
	my $max_count = 0;
	for(my $i = 0; $i<=1; $i++){
		my $base = $down[$i];
		my $count = 1;
		my $j = $i+1;
		while($base eq $down[$j]){
			$count++;
			$j++;
			last if($j > $#down);
		}
		if($count > $max_count){
			$max_count = $count;
		}
	}
	if($max_count >= 5){
		return(1);
	}

	$max_count = 0;
	for(my $i = 9; $i>=8; $i--){
		my $base = $up[$i];
		my $count = 1;
		my $j = $i-1;
		while($base eq $up[$j]){
			$count++;
			$j--;
			last if($j < 0);
		}
		if($count > $max_count){
			$max_count = $count;
		}
	}
	if($max_count >= 5){
		return(1);
	}	
	
	return(0);
}


sub calc_min_cov($$\%){
	my ($min_c, $min_p, $hash) = @_;
	my @covs;
	
	foreach my $c (keys %$hash){
		foreach my $p (keys %{$$hash{$c}}){
			my $cov = $$hash{$c}{$p}{"DP"};
			if($cov >= 0){
				push(@covs, $cov);
			}else{
				push(@covs, 0);
			}
		}
	}

	my $quart = int($#covs*$min_p);
	my @sorted_c = sort {$a <=> $b} @covs;
	my $min = $sorted_c[$quart];
	if($min > $min_c){
		return $min;
	}else{
		return $min_c;
	}
}


sub get_indel_information($$$$$$$){
	my ($chr, $pos, $ref, $indel, $type, $filt, $tmpfile) = @_;
	my $alt_ind;
	my $location;
	my $ref_ind;
	my $ind_fasta = "";
	
	if($type eq "DEL"){
		$ref_ind = $ref.$indel;
		$alt_ind = $ref;
		$location = $chr.":".($pos-11)."-".($pos+10+length($ref_ind)-1);
	}else{
		$alt_ind = $ref.$indel;
		$ref_ind = $ref;
		$location = $chr.":".($pos-11)."-".($pos+10);
	}
						
	`$abs_path/twoBitToFa -seq=$location -noMask $twoBit $tmpfile`;
	open(TMPIN, $tmpfile);
	while(my $line = <TMPIN>){
		chomp $line;
		next if($line =~ m/^>/);
		$ind_fasta .= $line;
	}
	close(TMPIN);
	`rm $tmpfile`;
	my $first_fasta = substr($ind_fasta, 0, 10);
	my $last_fasta = substr($ind_fasta, 11, 10);

	my $ind_flag = find_indel_repeat($first_fasta, $last_fasta);
	$filt = "HRUN" if($ind_flag == 1);
	
	return($ref_ind, $alt_ind, $filt);
}	

sub add_filter($$$\%){
	my ($chr, $pos, $filt, $hash) = @_;
	if($$hash{$chr}{$pos}{"FILTER"} eq "PASS"){
		$$hash{$chr}{$pos}{"FILTER"} = $filt;
	}else{
		$$hash{$chr}{$pos}{"FILTER"} .= ",$filt";
	}
}




sub get_prob_b_given_g($$$$){
	my ($g1, $g2, $b, $error) = @_;
	my $p = (0.5*(get_prob_b_given_a($g1,$b,$error)) + 0.5*(get_prob_b_given_a($g2,$b,$error)));
	return ($p);
}


sub get_prob_b_given_a($$$){
	my ($a, $b, $error) = @_;
	if($a eq $b){
		return (1-$error);
	}else{
		return ($error/3);
	}
}

sub convert_from_phred($){
	my ($qual) = @_;
	return (10**(-1*$qual/10));
}

sub log10($){
	my ($n) = @_;
	return (log($n)/log(10));
}

sub get_bed($\%){
	my ($bed, $hash) = @_;
	open(BED, $bed);
	while(my $line = <BED>){
		chomp $line;
		next if($line =~ m/^track/);
		next if($line =~ m/^#/);
		next if($line =~ m/^>/);	
		my ($chr, $start, $stop, @other) = split(/\t/, $line);
		$$hash{$chr}{$start} = $stop;
	}
	close(BED);
}


sub print_header($$){
	my ($min_norm, $min_tum) = @_;
	print "##fileformat=VCFv4.0\n";
	my @time = localtime(time);
	my $year = $time[5] + 1900;
	my $month = $time[4]+1;
	print "##fileDate=$year$month$time[3]\n";
	print "##source=Mutascope.pl callSomatic\n";
	print "##reference=hg19\n";
	print "##phasing=partial\n";
	print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total coverage depth\">\n";
	print "##INFO=<ID=BQ,Number=1,Type=Integer,Description=\"Average RMS base quality of all samples\">\n";
	print "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average RMS mapping quality of all samples\">\n";
	print "##INFO=<ID=ER,Number=1,Type=Integer,Description=\"Error rate used in the callSomatic module\">\n";
	print "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type, can be SNP, INS or DEL\">\n";
	print "##INFO=<ID=FEP,Number=1,Type=Float,Description=\"Somatic P-value calculated in the callSomatic module\">\n";
	print "##INFO=<ID=GQ,Number=1,Type=Integer,Description=\"Bayesian likelihood Genotype Quality score = Second Best Score - Best Score\">\n";
	print "##INFO=<ID=PL,Number=3,Type=Integer,Description=\"Bayesian likelihood Genotype Quality Scores for REF,HET,ALT genotypes\">\n";
	print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	print "##FORMAT=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n";
	print "##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of ref-forward, ref-reverse, alt-forward and alt-reverse bases\">\n";
	print "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Ref,alt coverage\">\n";
	print "##FORMAT=<ID=BQ2,Number=2,Type=Integer,Description=\"REF RMS base quality and ALT RMS base quality\">\n";
	print "##FORMAT=<ID=MQ2,Number=2,Type=Integer,Description=\"REF RMS mapping quality and ALT RMS mapping quality\">\n";
	print "##FORMAT=<ID=NR,Number=1,Type=Integer,Description=\"The total Number of Read groups overlapping the position\">\n";
	print "##FORMAT=<ID=DS,Number=1,Type=Integer,Description=\"The distance to the second allele (1st - 2nd)/(SUM(Non-Ref))\">\n";
	print "##FORMAT=<ID=ST,Number=1,Type=Integer,Description=\"The number of strands with the alternate allele (1 or 2)\">\n";
	print "##FORMAT=<ID=RGB,Number=1,Type=Float,Description=\"Read Group Bias score for the alternate allele using a chi-squared test\">\n";
	print "##FORMAT=<ID=SS,Number=1,Type=Integer,Description=\"Final validation status relative to non-adjacent Normal,0=wildtype,1=germline,2=somatic,3=LOH,4=post transcriptional modification,5=unknown\">\n";
	print "##FILTER=<ID=MIF,Description=\"Multiple INDELs at one position\">\n";
	print "##FILTER=<ID=HRUN,Description=\"5bp or more homopolymer repeat within 1bp of the start or end of the INDEL\">\n";
	print "##FILTER=<ID=MCF,Description=\"The coverage was below the minimum Coverage (Normal = $min_norm and Tumor = $min_tum)\">\n";
	print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$normal\t$tumor\n";
}



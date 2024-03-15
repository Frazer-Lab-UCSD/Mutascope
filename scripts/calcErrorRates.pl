#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);

sub check_outDir($);
sub check_bed($);
sub get_sub_type($$);
sub mean_and_stdev(@);

my $usage = "\nUse:\n\tMutascope calcErrorRates [-bed BED] [-dbsnp VCF] [-pdir PATH] [-sample NAME] [-BQ INT] [-MQ INT] [-noBlacklist] [-h]

Description:
\tcalcErrorRates generates the sample.errorRates file used to call variants. This module requires the sample.xpileup.gz file.

Example:
\tMutascope calcErrorRates -bed BED_NAME.bed -dbsnp dbsnp131.vcf -pdir /home/projects/ -sample patient1Normal 

Options:
\t-bed\t\tBED\tBED format file containing target locations with primer locations (See BED format in the Manual.pdf for more detail) [REQUIRED]
\t-dbsnp\t\tVCF\tdbSNP file in VCF format [REQUIRED]
\t-pdir\t\tPATH\tPath to the PROJECT directory [REQUIRED]
\t-sample\t\tSTRING\tName of the sample (only use letters or numbers! No symbols!) [REQUIRED]
\t-BQ\t\tINT\tMinimum Base Quality score (default 10) [Optional]
\t-MQ\t\tINT\tMinimum Mapping Quality score (default 1) [Optional]
\t-noBlacklist\tFLAG to not use/generate the blacklist file [Optional]
\t-h\t\tFLAG to output help message [Optional]
\n";

my $BQ = 10;
my $MQ = 1;
my $bed;
my $dbsnp;
my $sample;
my $pdir;
my $h;
my $noBlacklist;

GetOptions("bed=s" => \$bed,
						"dbsnp=s" => \$dbsnp,
						"BQ=i" => \$BQ,
						"MQ=i" => \$MQ,
						"sample=s" => \$sample,
						"pdir=s" => \$pdir,
						"h" => \$h,
						"noBlacklist" => \$noBlacklist);

die "$usage" if($h);
die "\n$usage\nERROR: Must provide BED file!\n\n" if(!$bed);
die "\n$usage\nERROR: Bed file \"$bed\" doesn't exist!\n\n" if(!(-e $bed));
die "\n$usage\nERROR: Must provide the full path to the Output Directory!\n\n" if(!$pdir);
die "\n$usage\nERROR: The Directory \"$pdir\" doesn't exists!\n\n" if(!(-d $pdir));
die "\n$usage\nERROR: Must provide dbSNP file!\n\n" if(!$dbsnp);
die "\n$usage\nERROR: dbSNP file \"$dbsnp\" doesn't exist!\n\n" if(!(-e $dbsnp));
die "\n$usage\nERROR: Must provide a dbSNP file in VCF format!\n\n" if(!($dbsnp =~ m/\.vcf/));
die "\n$usage\nERROR: Please Specify the name of the sample. This should be consistent throughout all Mutascope modules..\n\n" if(!$sample);
die "\n$usage\nERROR: The sample name can NOT have a '_' in the name.  Please use a different sample name other than \"$sample\".\n\n"if($sample =~ m/_/);
die "\n$usage\nERROR: The file $pdir\/intermediate/$sample.xpileup.gz doesn't exists!\nCheck to make sure you put in the correct sample name ($sample) along with the correct directory ($pdir).\n\nThis module assumes you ran Mutascope xpileup.\n\n" if(!(-e "$pdir\/intermediate/$sample.xpileup.gz"));

check_outDir($pdir);
check_bed($bed);

my @tmps = split(/\//, $bed);
my $bed_name = "$pdir\/intermediate/".$tmps[$#tmps];
$bed_name =~ s/.bed//;



die "\n$usage\nERROR: The black list file \"$bed_name.blacklist\" does NOT exist. Please make sure that this file is in the same directory as your BED file. If you do not have this file use Mutascope makeBlackList to generate it OR specify \"-noBlacklist\" in the command options.\n\n" if(!(-e "$bed_name.blacklist") and !$noBlacklist);


my %SNPs = ('A',0,'T',0,'C',0,'G',0);
my %dbsnp;
my %intervals;
my %skip;

open(BED, $bed);
while(my $line = <BED>){
	chomp $line;
	next if($line =~ m/^track/);
	next if($line =~ m/^#/);
	next if($line =~ m/^>/);	
	my ($chr, $start, $stop, @other) = split(/\t/, $line);
	$intervals{$chr}{$start} = $stop;
}

open(DBSNP, $dbsnp);
while(my $line = <DBSNP>){
	chomp $line;
	next if($line =~ m/^#/);
	my ($chr, $pos, $rsid, $ref, $alt, @other) = split(/\t/, $line);
	if($ref ne "A" and $ref ne "T" and $ref ne "C" and $ref ne "G"){
		#indel
	}elsif(($alt ne "A" and $alt ne "T" and $alt ne "C" and $alt ne "G") and !($alt =~ m/,/)){
		#indel
	}else{
INTLOOP: foreach my $start (sort {$a <=> $b} keys %{$intervals{$chr}}){
			my $stop = $intervals{$chr}{$start};
			if($pos >= $start and $pos <= $stop){
				$dbsnp{$chr}{$pos} = $alt;
			}
			last INTLOOP if($pos < $start);
		}
	}
}
close(DBSNP);


if(!$noBlacklist){
	open(FILE, "$bed_name.blacklist");
	while(my $line = <FILE>){
		chomp $line;
		my ($chr, $pos, $ref, $alt) = split(/\t/, $line);
		next if($ref eq "-");
		next if($alt eq "-");
		$skip{$chr}{$pos}{$alt} = 1;
	}
	close(FILE);
}

my %sub_types;
my %mixed_error;
my %mixed_covs;

my $xpileup_file = "$pdir\/intermediate/$sample.xpileup.gz";

if($xpileup_file =~ m/\.gz$/){
	open(PVCF, "gunzip -c $xpileup_file |") || die "Can't open pipe to $xpileup_file";
}else{
	open(PVCF, $xpileup_file) || die "Can't open $xpileup_file";
}

LOOP: while(my $line = <PVCF>){
	chomp $line;
	my ($chr, $pos, $ref, $cov, $indels, @rgs) = split(/\t/, $line);

	if(!(exists($dbsnp{$chr}{$pos}))){
		$ref =~ tr/[a-z]/[A-Z]/;
		foreach my $rg_info (@rgs){
			my ($rg, $bases_info, $bq_info, $mq_info, $bp_info) = split(/\s/, $rg_info);
			$bases_info =~ tr/[a-z]/[A-Z]/;
			my @bases = split(//, $bases_info);
			my @bq = split(//, $bq_info);
			my @mq = split(//, $mq_info);
			my @bp = split(/,/, $bp_info);
			my ($samp, $read, $loc) = split(/_/, $rg);
			
			for(my $i = 0; $i <= $#bases; $i++){
				my $base_qual = ord($bq[$i])-33;
				my $mq_qual = ord($mq[$i])-33;

				if($base_qual >= $BQ and $mq_qual >= $MQ){
					if(!(exists($skip{$chr}{$pos}{$bases[$i]}))){
						if($bases[$i] =~ m/[ATCG]/){
							$mixed_error{$read}{$bp[$i]}{get_sub_type($ref, $bases[$i])}++;
							$sub_types{get_sub_type($ref, $bases[$i])}++;
						}else{
							$mixed_covs{$read}{$bp[$i]}{$ref}++;
						}
					}
				}
			}			
		}
	}
}
close(PVCF);


open(MIX, ">$pdir\/intermediate/$sample.errorRates");
	foreach my $read (sort keys %mixed_error){
		foreach my $bp (sort {$a <=> $b} keys %{$mixed_error{$read}}){
			foreach my $tp (sort keys %sub_types){
				my ($ref1, $ref2, @other) = split(//, $tp);
				my $error = 0;
				my $ref1_cov = 0;
				my $ref2_cov = 0;
				$error = $mixed_error{$read}{$bp}{$tp} if(exists($mixed_error{$read}{$bp}{$tp}));
				$ref1_cov = $mixed_covs{$read}{$bp}{$ref1}/3 if(exists($mixed_covs{$read}{$bp}{$ref1}));
				$ref2_cov = $mixed_covs{$read}{$bp}{$ref2}/3 if(exists($mixed_covs{$read}{$bp}{$ref2}));
				if($ref1_cov + $ref2_cov != 0){
					print MIX "$read\t$bp\t$tp\t".$error/($ref1_cov+$ref2_cov+$error)."\t$error\t".($ref1_cov+$ref2_cov+$error)."\n";
				}
			}
		}
	}
close(MIX);






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



sub mean_and_stdev(@){
	my @vals = @_;
	
	if(not @vals){
		return(0,0);
	}
	
	my $total = 0;
	foreach (@vals){ $total += $_; }
	my $mean = $total/($#vals + 1);
	my $sqtot = 0;
	foreach (@vals){ $sqtot += ($mean - $_) ** 2; }
	my $std = ($sqtot/($#vals)) ** 0.5;
	return($mean, $std);	
}



sub check_bed($){
	my ($bed) = @_;
	open(BED, $bed);
	while(my $line = <BED>){
		chomp $line;
		next if($line =~ m/^track/);
		next if($line =~ m/^#/);
		next if($line =~ m/^>/);
		
		my ($chr, $start, $stop, $name, $score, $strand, $pstart, $pstop, @other) = split(/\t/, $line);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$pstart);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$pstop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$start);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if(!$stop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($start > $pstart);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($stop < $pstop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($pstart > $pstop);
		die "The BED line \"$line\" is not the proper format!\nSee the Manual.pdf for BED format details.\n\n" if($start > $stop);
		
	}
	close(BED);
}

sub check_outDir($){
	my ($od) = @_;
	if(!(-d "$od/results")){
		`mkdir $od/results`;
		print STDERR "Making \"$od/results\" directory!\n";
	}
	if(!(-d "$od/intermediate")){
		`mkdir $od/intermediate`;
		print STDERR "Making \"$od/intermediate\" directory!\n";
	}
	if(!(-d "$od/quality")){
		`mkdir $od/quality`;
		print STDERR "Making \"$od/quality\" directory!\n";
	}
}




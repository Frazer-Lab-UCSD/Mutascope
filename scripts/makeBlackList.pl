#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);


sub check_bed($);
sub check_for_bwa();
sub check_bwa_index($);
sub check_for_fasta_twobit($);
sub check_outDir($);


my $usage = "\nUse:\n\tMutascope makeBlackList [-bed BED] [-fasta FA] [-pdir PATH] [-bwa_path FILE] [-length INT] [-t INT] [-h]


Description:
\tmakeBlackList generates a list of mutations to exclude from analysis (blacklist file). This module identifies alternate alignment possibilities for reads created from the list of PCR amplicons (BED). Specific base-pair mismatches resulting from these alternate alignments are added to the file for further exclusion from the analysis. The file generated is BED.blacklist

Example:
\tMutascope makeBlackList -bed BED_NAME.bed -fasta hg19.fa -pdir /home/projects/ 

Options:
\t-bed\t\tBED\tBED format file containing target locations with primer locations (See BED format in the Manual.pdf for more detail) [REQUIRED]
\t-fasta\t\tFA\tFASTA file of reference genome [REQUIRED]
\t-pdir\t\tPATH\tPath to the PROJECTS directory [REQUIRED]
\t-bwa_path\tFILE\tBWA executable (full path) [Optional]
\t-length\t\tINT\tLength of sequencing reads for the preparation of the BlackList file (default 151) [Optional]
\t-t\t\tINT\tNumber of threads for BWA (default 1) [Optional]
\t-h\t\tFLAG to output help message [Optional]
\n";


my $bed;
my $fasta;
my $bwa_path;
my $length = 151;
my $h;
my $n1=15;
my $n2=100;
my $t=1;
my $pdir;

GetOptions("bed=s" => \$bed,
						"bwa_path=s" => \$bwa_path,
						"pdir=s" => \$pdir,
						"fasta=s" => \$fasta,
						"length=i" => \$length,
						"h" => \$h,
						"t=i" => $t);

die "$usage" if($h);
die "\n$usage\nERROR: Must provide BED file!\n\n" if(!$bed);
die "\n$usage\nERROR: Must provide FASTA file!\n\n" if(!$fasta);
die "\n$usage\nERROR: Bed file \"$bed\" doesn't exist!\n\n" if(!(-e $bed));
die "\n$usage\nERROR: Fasta file \"$fasta\" doesn't exist!\n\n" if(!(-e $fasta));
die "\n$usage\nERROR: Must provide the full path to the Output Directory!\n\n" if(!$pdir);
die "\n$usage\nERROR: The Directory \"$pdir\" doesn't exists!\n\n" if(!(-d $pdir));

check_outDir($pdir);
check_bed($bed);
check_bwa_index($fasta);
my $fasta_2bit = check_for_fasta_twobit($fasta);
my @tmps = split(/\//, $bed);
my $bed_name = "$pdir\/intermediate/".$tmps[$#tmps];
$bed_name =~ s/.bed//;

if(!$bwa_path){
	$bwa_path = check_for_bwa();
}else{
	die "\n$usage\nERROR: BWA: \"$bwa_path\" doesn't exists!\n\n" if(!(-e $bwa_path));
}
print STDERR "Using bwa: $bwa_path\n\n";


`perl $abs_path/bed_to_fastq.pl -length $length -bed $bed -twoBit $fasta_2bit -sub_name $bed_name`;

`$bwa_path aln -N -n $n1 -t $t $fasta $bed_name\_Forward_reads.fastq > $bed_name\_forward.sai`;
`$bwa_path samse -n $n2 $fasta $bed_name\_forward.sai $bed_name\_Forward_reads.fastq > $bed_name\_Forward_reads_aligned.sam`;
`$bwa_path aln -N -n $n1 -t $t $fasta $bed_name\_Reverse_reads.fastq > $bed_name\_reverse.sai`;
`$bwa_path samse -n $n2 $fasta $bed_name\_reverse.sai $bed_name\_Reverse_reads.fastq > $bed_name\_Reverse_reads_aligned.sam`;

`perl $abs_path/create_altAlignment_matrix_file.pl -twoBit $fasta_2bit $bed_name\_Forward_reads_aligned.sam > $bed_name\_Forward_reads_altAlignPos.txt`;	
`perl $abs_path/create_altAlignment_matrix_file.pl -twoBit $fasta_2bit $bed_name\_Reverse_reads_aligned.sam > $bed_name\_Reverse_reads_altAlignPos.txt`;

my %hash;
open(IN, "$bed_name\_Forward_reads_altAlignPos.txt");
while(my $line = <IN>){
	chomp $line;
	my ($chr, $pos, $ref, $alt) = split(/\t/, $line);
	$hash{$chr}{$pos}{$ref}{$alt} = 1;
}
close(IN);

open(IN, "$bed_name\_Reverse_reads_altAlignPos.txt");
while(my $line = <IN>){
	chomp $line;
	my ($chr, $pos, $ref, $alt) = split(/\t/, $line);
	$hash{$chr}{$pos}{$ref}{$alt} = 1;
}
close(IN);


open(OUT, ">$bed_name.blacklist");

foreach my $c (sort keys %hash){
	foreach my $p (sort {$a <=> $b} keys %{$hash{$c}}){
		foreach my $ref (sort keys %{$hash{$c}{$p}}){
			foreach my $alt (sort keys %{$hash{$c}{$p}{$ref}}){
				print OUT "$c\t$p\t$ref\t$alt\n";
			}
		}
	}
}
close(OUT);


`rm $bed_name.fa`;
`rm $bed_name\_Forward_reads.fastq`;
`rm $bed_name\_forward.sai`;
`rm $bed_name\_Forward_reads_aligned.sam`;
`rm $bed_name\_Forward_reads_altAlignPos.txt`;
`rm $bed_name\_Reverse_reads.fastq`;
`rm $bed_name\_reverse.sai`;
`rm $bed_name\_Reverse_reads_aligned.sam`;
`rm $bed_name\_Reverse_reads_altAlignPos.txt`;



















sub check_for_fasta_twobit($){
	my ($f) = @_;
	$f =~ s/.fa//;
	die "\nERROR: The FASTA 2bit file '$f.2bit' must be in the same folder as the $f.fa file!!!\n\nTo generate the 2bit file run the command: scripts/faToTwoBit -noMask $f.fa $f.2bit\n\n" if(!(-e "$f.2bit"));
	return("$f.2bit");
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


sub check_for_bwa(){
	my $bwa_result = `which bwa`;
	chomp $bwa_result;
	if($bwa_result eq ""){
		die "\nERROR: Can not find BWA. Please install BWA (http://bio-bwa.sourceforge.net/).\n   If installed either specify the path (-bwa_path) or create a link to your local path (Example: /usr/bin/).\n"; 
	}else{
		return $bwa_result;
	}
}


sub check_bwa_index($){
	my ($fasta) = @_;
	my $index_flag = 0;
	my $missing = "";
	if(!(-e "$fasta.amb")){
		$index_flag = 1;
		$missing = "$fasta.amb";
	}elsif(!(-e "$fasta.ann")){
		$index_flag = 1;
		$missing = "$fasta.ann";
	}elsif(!(-e "$fasta.bwt")){
		$index_flag = 1;
		$missing = "$fasta.bwt";
	}elsif(!(-e "$fasta.pac")){
		$index_flag = 1;
		$missing = "$fasta.pac";
	}elsif(!(-e "$fasta.sa")){
		$index_flag = 1;
		$missing = "$fasta.sa";
	}
	
	if($index_flag == 1){
		die "\nERROR:\tThe fasta file \"$fasta\" is not indexed properly by BWA. You are missing the \"$missing\" file\n\nPlease index the fasta file with bwa! Run: bwa index -a bwtsw $fasta\n\n";
		#`bwa index -a bwtsw $fasta`;
		
	}
}




sub check_outDir($){
	my ($od) = @_;
	if(!(-d "$od/results")){
		`mkdir $od/results`;
		print STDERR "Making \"$od/results\" directory!\n\n";
	}
	if(!(-d "$od/intermediate")){
		`mkdir $od/intermediate`;
		print STDERR "Making \"$od/intermediate\" directory!\n\n";
	}
	if(!(-d "$od/quality")){
		`mkdir $od/quality`;
		print STDERR "Making \"$od/quality\" directory!\n\n";
	}
}





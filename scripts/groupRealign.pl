#! /usr/bin/perl -w

use warnings;
use strict;
use Getopt::Long;
use Cwd 'abs_path';
my $abs_path = abs_path(__FILE__);
my @tmp = split(/\//, $abs_path);
pop(@tmp);
$abs_path = join("/", @tmp);


sub check_for_samtools();
sub check_outDir($);
sub get_intervalList($$);

my $usage = "\nUse:\n\tMutascope groupRealign [-bed BED] [-fasta FA] [-pdir PATH] [-normal NAME] [-tumor NAME] [-samtools_path FILE] [-h]


Description:
\tgroupRealign realigns the normal and tumor samples together using GATK.

Example:
\tMutascope groupRealign -bed BED_NAME.bed -fasta hg19.fa -pdir /home/projects/ -normal patient1Normal -tumor patient1Tumor  

Options:
\t-bed\t\tBED\tBED format file containing target locations with primer locations (See BED format in the Manual.pdf for more detail) [REQUIRED]
\t-fasta\t\tFA\tFASTA file of reference genome [REQUIRED]
\t-pdir\t\tPATH\tPath to the PROJECT directory [REQUIRED]
\t-normal\t\tSTRING\tName of the NORMAL sample (only use letters or numbers! No symbols!) [REQUIRED]
\t-tumor\t\tSTRING\tName of the TUMOR sample (only use letters or numbers! No symbols!) [REQUIRED]
\t-samtools_path\tFILE\tSAMTools executable (full path) [Optional]
\t-h\t\tFLAG to output help message [Optional]
\n";

my $fasta;
my $bed;
my $samtools_path;
my $pdir;
my $normal;
my $tumor;
my $h;


GetOptions("fasta=s" => \$fasta,
						"bed=s" => \$bed,
						"samtools_path=s" => \$samtools_path,
						"pdir=s" => \$pdir,
						"normal=s" => \$normal,
						"tumor=s" => \$tumor,
						"h" => \$h);

die "$usage" if($h);
die "\n$usage\nERROR: Must provide BED file!\n\n" if(!$bed);
die "\n$usage\nERROR: Bed file \"$bed\" doesn't exist!\n\n" if(!(-e $bed));
die "\n$usage\nERROR: Must provide FASTA file!\n\n" if(!$fasta);
die "\n$usage\nERROR: Fasta file \"$fasta\" doesn't exist!\n\n" if(!(-e $fasta));
die "\n$usage\nERROR: Must provide the full path to the Output Directory!\n\n" if(!$pdir);
die "\n$usage\nERROR: The Directory \"$pdir\" doesn't exists!\n\n" if(!(-d $pdir));
die "\n$usage\nERROR: Please Specify the name of the normal sample. This should be consistent throughout all Mutascope modules.\n\n" if(!$normal);
die "\n$usage\nERROR: Please Specify the name of the tumor sample. This should be consistent throughout all Mutascope modules.\n\n" if(!$tumor);
die "\n$usage\nERROR: The normal sample name can NOT have a '_' in the name.  Please use a different sample name other than \"$normal\".\n\n" if($normal =~ m/_/);
die "\n$usage\nERROR: The tumor sample name can NOT have a '_' in the name.  Please use a different sample name other than \"$tumor\".\n\n" if($tumor =~ m/_/);
die "\n$usage\nERROR: The file $pdir\/intermediate/$normal\_merged_unClipped_primersClipped.bam doesn't exists!\nCheck to make sure you put in the correct sample name ($normal).\n\nThis module assumes you ran Mutascope runBWA. If you did not you need to first run 'runBWA' before running this module\n\n" if(!(-e "$pdir\/intermediate/$normal\_merged_unClipped_primersClipped.bam"));
die "\n$usage\nERROR: The file $pdir\/intermediate/$tumor\_merged_unClipped_primersClipped.bam doesn't exists!\nCheck to make sure you put in the correct sample name ($tumor).\n\nThis module assumes you ran Mutascope runBWA. If you did not you need to first run 'runBWA' before running this module\n\n" if(!(-e "$pdir\/intermediate/$tumor\_merged_unClipped_primersClipped.bam"));


check_outDir($pdir);

if(!$samtools_path){
	$samtools_path = check_for_samtools();
}else{
	die "\n$usage\nERROR: SAMTools: \"$samtools_path\" doesn't exists!\n\n" if(!(-e $samtools_path));
}
print STDERR "Using SAMTools: $samtools_path\n";

my $interval_list = get_intervalList($bed, $pdir);


`$samtools_path view -H $pdir\/intermediate/$normal\_merged_unClipped_primersClipped.bam | perl $abs_path/get_readGroups.pl - > $pdir\/intermediate/$normal\_readGroups.txt`;
`$samtools_path view -H $pdir\/intermediate/$tumor\_merged_unClipped_primersClipped.bam | perl $abs_path/get_readGroups.pl - > $pdir\/intermediate/$tumor\_readGroups.txt`;

`mkdir $pdir/intermediate/$$/`;
`java -Xmx2048M -Djava.io.tmpdir=$pdir/intermediate/$$/ -jar $abs_path/MergeSamFiles.jar INPUT=$pdir\/intermediate/$normal\_merged_unClipped_primersClipped.bam INPUT=$pdir\/intermediate/$tumor\_merged_unClipped_primersClipped.bam OUTPUT=$pdir\/intermediate/$normal\_$tumor.bam ASSUME_SORTED=TRUE USE_THREADING=TRUE `;
`$samtools_path index $pdir\/intermediate/$normal\_$tumor.bam`;

`java -Xmx8G -Djava.io.tmpdir=$pdir/intermediate/$$/ -jar $abs_path/GenomeAnalysisTK.jar -T IndelRealigner --targetIntervals $interval_list -I $pdir\/intermediate/$normal\_$tumor.bam -R $fasta --out $pdir\/intermediate/$normal\_$tumor\_realigned.bam -maxInMemory 300000 --maxReadsForRealignment 500000 --LODThresholdForCleaning 2.5 --maxConsensuses 60 --maxReadsForConsensuses 240`;
`rm -rf $pdir/intermediate/$$`;

`$samtools_path view -bhR $pdir\/intermediate/$normal\_readGroups.txt $pdir\/intermediate/$normal\_$tumor\_realigned.bam > $pdir\/intermediate/$normal\_merged_unClipped_primersClipped_realigned.bam`;
`$samtools_path index $pdir\/intermediate/$normal\_merged_unClipped_primersClipped_realigned.bam`;

`$samtools_path view -bhR $pdir\/intermediate/$tumor\_readGroups.txt $pdir\/intermediate/$normal\_$tumor\_realigned.bam > $pdir\/intermediate/$tumor\_merged_unClipped_primersClipped_realigned.bam`;
`$samtools_path index $pdir\/intermediate/$tumor\_merged_unClipped_primersClipped_realigned.bam`;


`rm $pdir\/intermediate/$normal\_readGroups.txt`;
`rm $pdir\/intermediate/$tumor\_readGroups.txt`;
`rm $pdir\/intermediate/$normal\_$tumor.bam`;
`rm $pdir\/intermediate/$normal\_$tumor.bam.bai`;
`rm $pdir\/intermediate/$normal\_$tumor\_realigned.bam`;
`rm $pdir\/intermediate/$normal\_$tumor\_realigned.bai`;






sub get_intervalList($$){
	my ($bed, $pdir) = @_;
	
	my @tmps = split(/\//, $bed);
	my $interval_list = "$pdir\/intermediate/".$tmps[$#tmps];
	$interval_list =~ s/.bed//;
	$interval_list .= ".interval_list";
	
	if(-e $interval_list){
		return ($interval_list);
	}else{
		`perl $abs_path/make_intervalList_from_BED.pl $bed > $interval_list`;
		return ($interval_list);
	}
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

sub check_for_samtools(){
	my $samtools_results = `which samtools`;
	chomp $samtools_results;
	if($samtools_results eq ""){
		die "ERROR: Can not find SAMTools. Please install SAMTools (http://samtools.sourceforge.net/).\n   If installed either specify the path (-samtools_path) or create a link to your local path (Example: /usr/bin/).\n"; 
	}else{
		return $samtools_results;
	}
}









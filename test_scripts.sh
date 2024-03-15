#!/bin/sh

#########################################################
# AUTHOR: Shawn Yost <yostshawn@gmail.com>
#
# To run type "./test_scripts.sh"
# The script must be run in the directory containing the 
#     perl programs
#
# See the manual.pdf for more details on file formats and
#   help using the programs.
#
#########################################################

echo "#######  Try runBWA ########"
echo "# Running: ./Mutascope.pl runBWA -fasta test_files/hg19_chr22.fa -pdir test_output -sample NORM -r1_fastq test_files/normal_reads1.fastq -r2_fastq test_files/normal_reads2.fastq "

./Mutascope.pl runBWA -fasta test_files/hg19_chr22.fa -pdir test_output -sample NORM -r1_fastq test_files/normal_reads1.fastq -r2_fastq test_files/normal_reads2.fastq

echo "#"
echo "# indexing the test genome"

bwa index -a bwtsw test_files/hg19_chr22.fa

if [ "${?}" -ne "0" ] ; then
    echo "BWA is not found!!
    Please make sure BWA is installed and in your path!
    "
    exit ${?}
fi

echo "#
# Running 'runBWA' again!"

./Mutascope.pl runBWA -fasta test_files/hg19_chr22.fa -pdir test_output -sample NORM -r1_fastq test_files/normal_reads1.fastq -r2_fastq test_files/normal_reads2.fastq

if [ "${?}" -ne "0" ] ; then
    echo "SAMTools is not found!!
    Please make sure SAMTools is installed and in your path!
    "
    exit ${?}
fi

echo "#
# Running Tumor sample through 'runBWA'"

./Mutascope.pl runBWA -fasta test_files/hg19_chr22.fa -pdir test_output -sample TUM -r1_fastq test_files/tumor_reads1.fastq -r2_fastq test_files/tumor_reads2.fastq

echo "#
#
# Mutascope.pl runBWA works great!"
echo "#
# Running Mutascope.pl makeBlackList -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output"

./Mutascope.pl makeBlackList -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output

echo "#
# Need to make 2bit file of the genome!"

faToTwoBit -noMask test_files/hg19_chr22.fa test_files/hg19_chr22.2bit

echo "# Running Mutascope.pl makeBlackList Again!"

./Mutascope.pl makeBlackList -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output


echo "#
# Running Mutascope.pl refinement -bed test_files/target_locations_withPrimers_chr22.bed -pdir test_output -sample NORM"

./Mutascope.pl refinement -bed test_files/target_locations_withPrimers_chr22.bed -pdir test_output -sample NORM

echo "#
# Running Mutascope.pl refinement -bed test_files/target_locations_withPrimers_chr22.bed -pdir test_output -sample TUM"

./Mutascope.pl refinement -bed test_files/target_locations_withPrimers_chr22.bed -pdir test_output -sample TUM

echo "#
# Running: ./Mutascope.pl groupRealign -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output -normal NORM -tumor TUM"

./Mutascope.pl groupRealign -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output -normal NORM -tumor TUM

echo "#
# Running: ./Mutascope.pl xpileup -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output -sample NORM"

./Mutascope.pl xpileup -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output -sample NORM

echo "#
# Running: ./Mutascope.pl xpileup -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output -sample TUM"

./Mutascope.pl xpileup -bed test_files/target_locations_withPrimers_chr22.bed -fasta test_files/hg19_chr22.fa -pdir test_output -sample TUM


echo "#
# Running: ./Mutascope.pl calcErrorRates -bed test_files/target_locations_withPrimers_chr22.bed -dbsnp test_files/dbsnp_135_chr22.vcf -pdir test_output -sample NORM"

./Mutascope.pl calcErrorRates -bed test_files/target_locations_withPrimers_chr22.bed -dbsnp test_files/dbsnp_135_chr22.vcf -pdir test_output -sample NORM



echo "#
# Running: ./Mutascope.pl callSomatic -bed test_files/target_locations_withPrimers_chr22.bed -dbsnp test_files/dbsnp_135_chr22.vcf -pdir test_output -fasta test_files/hg19_chr22.fa -normal NORM -tumor TUM "

./Mutascope.pl callSomatic -bed test_files/target_locations_withPrimers_chr22.bed -dbsnp test_files/dbsnp_135_chr22.vcf -pdir test_output -fasta test_files/hg19_chr22.fa -normal NORM -tumor TUM 


echo "#
#
#######################################
#                                     #
#    MUTASCOPE is ready to run!!!     #
#                                     #
#######################################"


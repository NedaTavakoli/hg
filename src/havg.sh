#!/bin/bash
# $1 - chromosome

# Change this part according to your settings
project_dir=/storage/coda1/p-saluru8/0/ntavakoli6/hg
software_dir=/storage/coda1/p-saluru8/0/ntavakoli6/software
DATA=/storage/coda1/p-saluru8/0/ntavakoli6/hg/data

cd ${project_dir}
bcftools=/storage/coda1/p-saluru8/0/ntavakoli6/software/bcftools/bcftools
vcftools=/storage/coda1/p-saluru8/0/ntavakoli6/software/vcftools-0.1.16/bin/vcftools
tabix=/storage/coda1/p-saluru8/0/ntavakoli6/software/tabix/tabix

cd ${DATA}

# https://www.biostars.org/p/298361/

id=${1}

############################################################################################

# keep only snp positions
$vcftools --vcf  chr${id}.vcf --remove-indels --recode --recode-INFO-all --out chr${id}_snps_only

#output  out.frq.count
$vcftools --vcf chr${id}.vcf --counts --remove-indels 
mv out.frq.count chr${id}.out.frq.count


paste <(cat chr${id}_snps_only.recode.vcf |awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
        !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
    \
    <(cat chr${id}.out.frq.count |\
        awk -F"\t"  '{print $3"\t"$5"\t"$6"\t"$7"\t"$8}') \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat.txt

#*******************************************************************************************************
#The obtained file has these colomns: 
#  $1   $2  $3  $4   $5     $6          $7                $8                $9            $10
# CHR  POS  ID  REF  ALT  N_ALLELE  {REF:count}  {Allele1:count}  {Allele2:count}   {Allele2:count}
#********************************************************************************************************

paste <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nALT1_sample1"} {print gsub(/1\|0|1\/0|1\|1|1\/1|1\|2|1\/2|1\|3|1\/3/, "")}') \
        \
    <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nALT1_sample2"} {print gsub(/0\|1|0\/1|1\|1|1\/1|2\|1|2\/1|3\|1|3\/1/, "")}') \
        \
    <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nALT2_sample1"} {print gsub(/2\|0|2\/0|2\|1|2\/1|2\|2|2\/2|2\|3|2\/3/, "")}') \
        \
    <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nALT2_sample2"} {print gsub(/0\|2|0\/2|1\|2|1\/2|2\|2|2\/2|3\|2|3\/2/, "")}') \
        \
     <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nALT3_sample1"} {print gsub(/3\|0|3\/0|3\|1|3\/1|3\|2|3\/2|3\|3|3\/3/, "")}') \
        \
     <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nALT3_sample2"} {print gsub(/0\|3|0\/3|1\|3|1\/3|2\|3|2\/3|3\|3|3\/3/, "")}') \
        \
     <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nREF_sample1"} {print gsub(/0\|0|0\/0|0\|1|0\/1|0\|2|0\/2|0\|3|0\/3/, "")}') \
        \
     <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.recode.vcf |\
        awk 'BEGIN {print "nREF_sample2"} {print gsub(/0\|0|0\/0|1\|0|1\/0|2\|0|2\/0|3\|0|3\/0/, "")}') \
        \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat_2.txt

#*********************************************************************************************************************************************************************
#The obtained file has these colomns: 
#     $11               $12               $13            $14               $15           $16                $17          $18
# nALT1_sample1    nALT1_sample2    nALT2_sample1   nALT2_sample2    nALT3_sample1   nALT3_sample2    nREF_sample1     nREF_sample2
#   1|0                0|1               2|0             0|2              3|0            0|3              0|0              0|0
#   1|1                1|1               2|1             1|2              3|1            1|3              0|1              1|0
#   1|2                2|1               2|2             2|2              3|2            2|3              0|2              2|0
#   1|3                3|1               2|3             3|2              3|3            3|3              0|3              3|0
#*********************************************************************************************************************************************************************    

paste <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_1_ALT1"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|0|1\/0|1\|1|1\/1|1\|2|1\/2|1\|3|1\/3/, "", $(i))==1) {printf header[i]"_1,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_2_ALT1"} \
    !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|0\/1|1\|1|1\/1|2\|1|2\/1|3\|1|3\/1/, "", $(i))==1) {printf header[i]"_2,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_1_ALT2"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/2\|0|2\/0|2\|1|2\/1|2\|2|2\/2|2\|3|2\/3/, "", $(i))==1) {printf header[i]"_1,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_2_ALT2"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|2|0\/2|1\|2|1\/2|2\|2|2\/2|3\|2|3\/2/, "", $(i))==1) {printf header[i]"_2,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_1_ALT3"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/3\|0|3\/0|3\|1|3\/1|3\|2|3\/2|3\|3|3\/3/, "", $(i))==1) {printf header[i]"_1,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_2_ALT3"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|3|0\/3|1\|3|1\/3|2\|3|2\/3|3\|3|3\/3/, "", $(i))==1) {printf header[i]"_2,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_1_REF"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0|0\|1|0\/1|0\|2|0\/2|0\|3|0\/3/, "", $(i))==1) {printf header[i]"_1,"}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "Samples_2_REF"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0|1\|0|1\/0|2\|0|2\/0|3\|0|3\/0/, "", $(i))==1) {printf header[i]"_2,"}; if (i==NF) {printf "\n"}}}') \
        \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat3.txt

#*********************************************************************************************************************************************************************
#The obtained file has these colomns: 
#   $19                    $20                 $21              $22              $23                $24                 $25           $26
# Samples_1_ALT1      Samples_2_ALT1     Samples_1_ALT2     Samples_2_ALT2    Samples_1_ALT3    Samples_2_ALT3    Samples_1_REF   Samples_2_REF
#*********************************************************************************************************************************************************************    

paste chr${id}_concat.txt chr${id}_concat_2.txt > chr${id}_concat_concat2.txt
paste chr${id}_concat.txt chr${id}_concat_2.txt chr${id}_concat3.txt > chr${id}_concat_concat2_concat3.txt
# remove the first line 
sed '1d' chr${id}_concat_concat2_concat3.txt > tmpfile; mv tmpfile  chr${id}_concat_concat2_concat3_no_header.txt
cat chr${id}_concat_concat2_concat3_no_header.txt | awk -F"\t" 'BEGIN {print "POS\tN_ALLELE\tREF\tALT\t{SamplesRef}\t{SamplesALT1}\t{SamplesALT2}\t{SamplesALT3}"}\
{print $2"\t"$6"\t"$4"\t"$5"\t"$25,$26"\t"$19,$20"\t"$21,$22"\t"$23,$24}' > chr${id}_haplotypes_frq_all_strings.txt


#*********************************************************************************************************************************************************************
#The obtained file has these colomns: 
#
# POS   N_ALLELLE   REF   ALT {SamplesRef} {SamplesALT1}  {SamplesALT2} {SamplesALT3}
#*********************************************************************************************************************************************************************


# keeping one copys of this file to make sure not loosing it
cp chr${id}_haplotypes_frq_all_strings.txt chr${id}_haplotypes_frq_all_strings_bc.txt

##### remove the first line 
sed '1d'  chr${id}_haplotypes_frq_all_strings.txt  > tmpfile; mv tmpfile chr${id}_haplotypes_frq_all_strings.txt

# Add line numbers to file
# nl chr${id}_concat.txt > chr${id}_concat_numbered.txt
# nl chr${id}_concat_2.txt > chr${id}_concat_2_numbered.txt


#### count the number of coloms, to find N

# N: number of all the edges
awk '{s+=$2}END{print s}' chr${id}_haplotypes_frq_all_strings.txt
#2123135  for chr22 ---> N-n = 1063618
#  12449093  for chr1 ----> N-n =12449093 -6215039 =6234054


# n : number of snp positions
cat chr22_concat_concat2_concat3_no_header.txt| wc -l
#1059517 for chr22
# 6215039 for chr1

####################################

# List of samples  #printing a list of samples from a VCF:
$bcftools query -l chr${id}.vcf > chr${id}_sample.txt

# Add line numbers to file
nl chr${id}_sample.txt > chr${id}_sample_numbered.txt

# Copy each line of sample and rename to sample_1 and smaple_2
while read line; do for i in {1..1}; do echo -e "$line"_1"\n"$line"_2";done; done < chr${id}_sample.txt  > chr${id}_sample_dup.txt

# Add line numbers to file
nl chr${id}_sample_dup.txt > chr${id}_sample_dup_numbered.txt

#####################################

# get list os SNP positions
# list of snps posiions
cat  chr${id}.out.frq.count | cut -f2 >chr${id}_snp_positions.txt

# remove the first line (POS)
sed '1d'  chr${id}_snp_positions.txt  > tmpfile; mv tmpfile chr${id}_snp_positions.txt


cat chr${id}_snp_positions.txt |wc -l
#n = 1059517

# Add line numbers to file
nl chr${id}_snp_positions.txt > chr${id}_snp_positions_numbered.txt

#####################################
 # get number of allels per position
cat  chr${id}.out.frq.count | cut -f3 >chr${id}_num_allels_per_positions.txt
# remove the first line
sed '1d' chr${id}_num_allels_per_positions.txt  > tmpfile; mv tmpfile chr${id}_num_allels_per_positions.txt

# Add line numbers to file
nl chr${id}_num_allels_per_positions.txt > chr${id}_num_allels_per_positions_numbered.txt


###******************************************

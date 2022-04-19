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

# Take care of multi alleles snps
$bcftools norm -m-any chr${id}.vcf -Ov > chr${id}.Norm.vcf


# keep only snp positions
$vcftools --vcf  chr${id}.Norm.vcf --remove-indels --recode --recode-INFO-all --out chr${id}_snps_only.Norm


#output  out.frq.count
$vcftools --vcf chr${id}.Norm.vcf --counts --remove-indels 

mv out.frq.count chr${id}.out.frq.count

paste <(cat chr${id}_snps_only.Norm.recode.vcf |awk -F"\t" 'BEGIN {print "CHR\tPOS\tID\tREF\tALT"} \
        !/^#/ {print $1"\t"$2"\t"$3"\t"$4"\t"$5}') \
    \
    <(cat out.frq.count |\
        awk -F"\t"  '{print $3"\t"$5"\t"$6}') \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat.txt
        
        
 
paste <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.Norm.recode.vcf |\
    awk 'BEGIN {print "nHet"} {print gsub(/0\|1|1\|0|0\/1|1\/0/, "")}') \
    \
    <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.Norm.recode.vcf |\
        awk 'BEGIN {print "nHomAlt"} {print gsub(/1\|1|1\/1/, "")}') \
        \
    <($bcftools query -f '[\t%SAMPLE=%GT]\n' chr${id}_snps_only.Norm.recode.vcf |\
        awk 'BEGIN {print "nHomRef"} {print gsub(/0\|0|0\/0/, "")}') \
    \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat_2.txt
        

paste <($bcftools view chr${id}_snps_only.Norm.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HetSamples"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|1|1\|0|0\/1|1\/0/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.Norm.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesAlt"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/1\|1|1\/1/, "", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
        \
    <($bcftools view chr${id}_snps_only.Norm.recode.vcf | awk -F"\t" '/^#CHROM/ {split($0, header, "\t"); print "HomSamplesRef"} \
        !/^#CHROM/ {for (i=10; i<=NF; i++) {if (gsub(/0\|0|0\/0/,"", $(i))==1) {printf header[i]","}; if (i==NF) {printf "\n"}}}') \
        \
    | sed 's/,\t/\t/g' | sed 's/,$//g'   > chr${id}_concat3.txt



paste chr${id}_concat.txt chr${id}_concat_2.txt > chr${id}_concat_concat2.txt
paste chr${id}_concat.txt chr${id}_concat_2.txt chr${id}_concat3.txt > chr${id}_concat_concat2_concat3.txt
cat chr${id}_concat_concat2_concat3.txt | awk -F"\t" '{print $2"\t"$6"\t"$4":{"$14"}""\t"$5":{"$12","$13"}"}' > chr${id}_haplotypes_frq_all_strings.txt


#######
#### count the number of coloms, to find N

# N
awk '{s+=$2}END{print s}' chr${id}_haplotypes_frq_all_strings.txt
#2128158

# n 
cat chr${id}_concat_concat2.txt | wc -l
#1064080

# List of samples  #printing a list of samples from a VCF:
$bcftools query -l chr${id}.vcf > chr${id}_sample.txt


# Add line numbers to file
nl chr${id}_sample.txt > chr${id}_sample_numbered.txt

# I am keeping two copies of this file to make sure of not loosing it
cp chr${id}_haplotypes_frq_all_strings.txt chr${id}_haplotypes_frq_all_strings_bc.txt
cp chr${id}_haplotypes_frq_all_strings.txt chr${id}_haplotypes_frq_all_strings_bc2.txt



#### replace sample names with number

readarray index < <(cat chr${id}_sample_numbered.txt | cut -f1)
#echo ${index[5]}
readarray  haplotypes < <(cat chr${id}_sample_numbered.txt | cut -f2)


for k in $(seq 0 2503);
do 
    num=$(($k+1))
    r=${haplotypes[$k]}
    echo $num
    rnospace=$(echo $r | sed -e 's/^[[:space:]]*//')   # To remove space at the end
    sed -i "s/${rnospace}/$num/g" chr${id}_haplotypes_frq_all_strings.txt
done

mv chr${id}_haplotypes_frq_all_strings.txt chr${id}_haplotypes_frq_all_numbers.txt

# get list os SNP positions
# list of snps posiions
cat  chr${id}.out.frq.count | cut -f2 >chr${id}_snp_positions.txt

# remove the first line (POS)
sed '1d'  chr${id}_snp_positions.txt  > tmpfile; mv tmpfile chr${id}_snp_positions.txt


cat chr${id}_snp_positions.txt |wc -l
#n = 1064080


# bulid matrix A

# 1. extract required informations from chr${id}_haplotypes_frq_all_numbers file the 4th colomn
#keep a copy of the file that has 4 colomns
cp chr${id}_haplotypes_frq_all_numbers.txt  chr${id}_haplotypes_frq_all_numbers_bc.txt 


# Extract colomn 4 form haplotypes as it represents the list of haplotypes for alt
cat chr${id}_haplotypes_frq_all_numbers.txt | cut -f4 > col4_chr${id}_haplotypes_frq_all_numbers.txt
cp col4_chr${id}_haplotypes_frq_all_numbers.txt col4_chr${id}_haplotypes_frq_all_numbers_bc.txt


# https://www.theunixschool.com/2014/08/sed-examples-remove-delete-chars-from-line-file.html
# To remove 1st n characters of every line:
sed -r 's/.{3}//' col4_chr${id}_haplotypes_frq_all_numbers.txt > col4_chr${id}_haplotypes_frq_all_numbers_1.txt

# To remove 1st n characters of every line:
sed -r 's/.{2}$//' col4_chr${id}_haplotypes_frq_all_numbers_1.txt > col4_chr${id}_haplotypes_frq_all_numbers_2.txt

# Each line of this file show a list of haplotypes per position, and the associated values 


# remove the first line from the file
sed '1d'  col4_chr${id}_haplotypes_frq_all_numbers_2.txt  > tmpfile; mv tmpfile col4_chr${id}_haplotypes_frq_all_numbers_2.txt

# https://linuxize.com/post/how-to-read-a-file-line-by-line-in-bash/
# Reading a File Line By Line Syntax #
while IFS= read -r line; do printf '%s\n' "$line"; done < col4_chr${id}_haplotypes_frq_all_numbers_2.txt

# Remove fist comma on each line
sed -i.bak 's/,//' col4_chr${id}_haplotypes_frq_all_numbers_2.txt






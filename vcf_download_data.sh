#!/bin/bash
# ssh ntavakoli6@login-phoenix-3.pace.gatech.edu
# qsub -I -q inferno -A GT-saluru8-CODA20 -l nodes=1:ppn=1,mem=300gb,walltime=92:00:00

# Change this part according to your settings
cd /storage/coda1/p-saluru8/0/ntavakoli6/hg
VG=/storage/coda1/p-saluru8/0/ntavakoli6/software/vg
DATA=/storage/coda1/p-saluru8/0/ntavakoli6/hg/data
REFERENCE=${DATA}/hs37d5.fa
HAPLOTYPEPATH=/storage/coda1/p-saluru8/0/ntavakoli6/hg/haplotypes
VGgraph=/storage/coda1/p-saluru8/0/ntavakoli6/hg/vg

cd data

# Get the reference
REFERENCE=hs37d5.fa
rm -f ${REFERENCE} ${REFERENCE}.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/${REFERENCE}.gz
gunzip ${REFERENCE}

# Get the variants
VARIANTS=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
SHORTV=all.vcf.gz
rm -f ${SHORTV} ${SHORTV}.tbi
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${VARIANTS}
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${VARIANTS}.tbi
mv ${VARIANTS} ${SHORTV}
mv ${VARIANTS}.tbi ${SHORTV}.tbi

# Get the phasings for chromosomes 1-22
PREFIX=ALL.chr
SPREFIX=chr
SUFFIX=.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
SSUFFIX=.vcf.gz
for i in $(seq 1 22; echo X; echo Y)
do
  NAME=${PREFIX}${i}${SUFFIX}
  SHORT=${SPREFIX}${i}${SSUFFIX}
  rm -f ${SHORT} ${SHORT}.tbi
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${NAME}
  wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/${NAME}.tbi
  mv ${NAME} ${SHORT}
  mv ${NAME}.tbi ${SHORT}.tbi
done

# Get the phasings for chromosomes X and Y
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi
mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz chrX.vcf.gz
mv ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz.tbi chrX.vcf.gz.tbi
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi
mv ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz chrY.vcf.gz
mv ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz.tbi chrY.vcf.gz.tbi


#------------------------------Data Done-------------------------------------------

# Unzip all the chromosomes
# create .vcf for all chromosomes
for id in $(seq 1 22; echo X; echo Y)
do
  cp chr${id}.vcf.gz chr${id}_bc.vcf.gz
  gzip -d chr${id}.vcf.gz
  mv chr${id}_bc.vcf.gz chr${id}.vcf.gz
done
#!/bin/bash
#Purpose: download and install softwares, datasets, install and run VF, FORGe, and HISAT2

mainwd=$(pwd)  #project top-level directory
mkdir -p software && cd software
softwarewd=$(pwd)

#get vcftools 
echo "downloading vcftools"
wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz
tar xzf vcftools-0.1.16.tar.gz && cd vcftools-0.1.16
./autogen.sh
./configure --prefix=$(pwd)
make -j -s && make install
vcftools=$(pwd)
rm -f "vcftools-0.1.16.tar.gz"
echo "vcftools download and compilation finished"

#get Gurobi
echo "downloading gurobi"
cd $softwarewd
wget  https://packages.gurobi.com/9.1/gurobi9.1.0_linux64.tar.gz
tar xzf gurobi9.1.0_linux64.tar.gz
make -j -C gurobi910/linux64/src/build #re-compile gurobi cpp files using user's c++ compiler
cp gurobi910/linux64/src/build/libgurobi_c++.a gurobi910/linux64/lib
gurobi=$(pwd)
rm -f "gurobi9.1.0_linux64.tar.gz"
echo "gurobi download and compilation finished"

#check if appropriate files exist
if [ ! -f gurobi910/linux64/include/gurobi_c++.h ]; then echo "gurobi download failed"; fi
if [ ! -f gurobi910/linux64/src/build/libgurobi_c++.a ]; then echo "gurobi compilation failed"; fi
if [ ! -f vcftools-0.1.16/bin/vcftools ]; then echo "vcftools compilation failed"; fi

echo "Looks like it went okay, now run <make>"
#Next, run make
#get VF

#get htslib
#Note: tabix, htslib and bgzip2 will be installed in the bin directory and rhe main directory
echo "downloading htslib"
cd $softwarewd
wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2
tar xvjf htslib-1.12.tar.bz2
cd htslib-1.12
autoreconf -i  # Build the configure script and install files it uses
./configure  --prefix=$softwarewd/htslib-1.12  # Optional but recommended, for choosing extra functionality
make
make install
htslib=$(pwd)
bgzip2=$(pwd)
tabix=$(pwd)
rm -f "htslib-1.12.tar.bz2"
echo "htslib download and compilation finished"

#get bcftools
#Note that bcftools will be installed in the bin directory of bcftools folder
cd $softwarewd
git clone https://github.com/samtools/bcftools.git
cd bcftools
autoreconf -i  # Build the configure script and install files it uses
./configure  --prefix=$softwarewd/bcftools  # Optional but recommended, for choosing extra functionality
make
make install 
bcftools=$(pwd)
echo "bcftools download and compilation finished"

#get SAMtools
echo "downloading SAMtools"
cd $softwarewd
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar -xvf samtools-1.12.tar.bz2
cd samtools-1.12
autoheader            # Build config.h.in (this may generate a warning about # AC_CONFIG_SUBDIRS - please ignore it).
autoconf -Wno-syntax  # Generate the configure script
./configure   --prefix=$softwarewd/samtools-1.12        # Needed for choosing optional functionality
make
make install
samtools=$(pwd)
rm -f "samtools-1.12.tar.bz2"
echo "samtools download and compilation finished"

#get Vg executable
echo "downloading Vg"
cd $softwarewd
wget https://github.com/vgteam/vg/releases/download/v1.31.0/vg
chmod +x vg
vg=$(pwd)
echo "vg download and compilation finished"

#get mason executable
cd $softwarewd
wget http://packages.seqan.de/mason/mason-0.1.2-Linux-x86_64.tar.bz2
tar -xvf mason-0.1.2-Linux-x86_64.tar.bz2 && cd mason-0.1.2-Linux-x86_64/bin
mason=$(pwd)
rm -f "mason-0.1.2-Linux-x86_64.tar.bz2"
echo "mason download and compilation finished"

#get Jellyfish
#note: python is needed
cd $softwarewd
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz
tar -xvf jellyfish-2.2.6.tar.gz && cd Jellyfish
#module load anaocnda2
./configure --prefix=$HOME --enable-python-binding
make -j 4
make install
jellyfish=$(pwd)
rm -f "jellyfish-2.2.6.tar.gz"
echo "jellyfish download and compilation finished"

#get k8
mainwd=$softwarewd  #project top-level directory
mkdir -p k8 && cd k8
# download the k8 executable
wget -O- http://sourceforge.net/projects/biobin/files/devtools/k8-0.2.1.tar.bz2/download \
	| bzip2 -dc | tar xf - 
k8=$(pwd)
echo "k8 download and compilation finished"

#get FORGe
cd $softwarewd
git clone --recursive https://github.com/langmead-lab/FORGe.git
cd FORGe
FORGe=$(pwd)

#get HISAT2
echo "downloading HISAT2"
cd $softwarewd
git clone https://github.com/infphilo/hisat2 hisat-genotype-top
cd hisat-genotype-top
hisat2=$(pwd)
git checkout hisat2_v2.2.0_beta
make hisat2-align-s hisat2-build-s hisat2-inspect-s
export PATH=hisat-genotype-top:=$hisat/hisatgenotype_scripts:$PATH
export PYTHONPATH=$hisat/hisatgenotype_modules:$PYTHONPATH
source ~/.bashrc
echo "HISAT2 download and compilation finished"

#get FORGe experiments
#FORGe experiments
echo "downloading FORGe experiments"
cd $softwarewd
git clone --recursive https://github.com/langmead-lab/FORGe-experiments.git
HSA=$softwarewd/hisat-genotype-top/hisat2
EXP_HOME =$softwarewd/FORGe-experiments
MASON =$softwarewd/mason-0.1.2-Linux-x86_64/bin
VIS_HOME =$softwarewd/FORGe/src
SCRIPT_HOME =$softwarewd/FORGe-experiments/scripts
echo "FORGe experiment download and compilation finished"
#LD_LIBRARY_PATH=/home/langmead/.linuxbrew/lib:$LD_LIBRARY_PATH
#PYTHONPATH=/usr/lib/oracle/11.2/lib64/python2.7/site-packages:$PYTHONPATH
#ln /usr/bin/sbatch sbatch

#get VF
echo "downloading VF"
cd $softwarewd
git clone git@github.com:AT-CG/VF.git
cd VF
./dependencies.sh
make
VF=$(pwd)
echo "VF download and compilation finished"


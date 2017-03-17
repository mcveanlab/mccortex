sudo apt-get update
sudo apt-get install -y g++ libncurses5-dev python-dev python3-dev emacs cmake autoconf

# Stampy
cd
curl -O http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz
tar xfz stampy-latest.tgz
cd stampy
make

# VCFTools
cd
wget https://downloads.sourceforge.net/project/vcftools/vcftools_0.1.13.tar.gz
tar xfz vcftools_0.1.13.tar.gz
cd vcftools_0.1.13
make

# Cortex
cd
git clone --recursive https://github.com/iqbal-lab/cortex.git
cd cortex
bash install.sh
for k in 31 63 95 127; do
  for ncol in 1 2 3 9 10 11; do
    make cortex_var MAXK=$k NCOLS=$ncol
  done
done
echo 'export PERL5LIB="${HOME}/cortex/scripts/analyse_variants/bioinf-perl/lib/:${HOME}/cortex/scripts/calling/:${PERL5LIB}"' >> .profile
echo 'export PATH="${HOME}/cortex/scripts/analyse_variants/needleman_wunsch/:${PATH}"' >> .profile

# McCortex
cd
git clone --recursive -b develop https://github.com/mcveanlab/mccortex.git
cd mccortex
cd libs && make all && cd ..
for k in 31 63 95 127; do
  make all test MAXK=31
done

# Freebayes
cd
git clone --recursive https://github.com/ekg/freebayes.git
cd freebayes
make

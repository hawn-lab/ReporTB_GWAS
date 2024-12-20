#### QC and ANALYSIS ####
## m4.16xlarge
## AWS AMI Linux 2
## 16GB base, 1TB EBS
ssh -i ~/Documents/AWS/keys/reportTB.pem ec2-user@

# Define access keys
AWS_ACCESS_KEY=""
AWS_SECRET_ACCESS_KEY=""
AWS_REGION=""

#### Basic AWS update ####
sudo yum upgrade -y
sudo yum update -y

## Install AWS command line client
export AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY
export AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY
export AWS_DEFAULT_REGION=$AWS_REGION

### Setup fuse
sudo amazon-linux-extras install -y epel
sudo yum install -y s3fs-fuse
## Setup fuse keys
echo $AWS_ACCESS_KEY:$AWS_SECRET_ACCESS_KEY > ~/.passwd-s3fs
chmod 600 ~/.passwd-s3fs

# Setup EBS
ebs_name="xvdb"
## Format volume
sudo mkfs -t ext4 /dev/$ebs_name
## Attach SEAsnake directory to volume
sudo mkdir -p ~/project
sudo mount /dev/$ebs_name ~/project
## Change permissions to read-write
sudo chmod 777 -R ~/project/

# Data
sudo mkdir -m 777 -p ~/project/data
s3fs hawn-report-tb-vantage ~/project/data -o passwd_file=~/.passwd-s3fs \
    -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007 -o compat_dir
  
## If you need to switch machines
# sudo mkdir -m 777 -p ~/project/result
# s3fs hawn-reporttb-results ~/project/result -o passwd_file=~/.passwd-s3fs_2 \
#     -o default_acl=public-read -o uid=1000 -o gid=1000 -o umask=0007 -o compat_dir
# fusermount -u ~/project/result

#############

# Install GNU parallel
# Tange, O. (2023, March 22). GNU Parallel 20230322 ('Arrest Warrant'). Zenodo. https://doi.org/10.5281/zenodo.7761866
mkdir ~/apps
cd ~/apps
wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2
sudo tar xjf parallel-latest.tar.bz2
cd parallel-20241022
sudo ./configure && make
sudo make install
## Test
parallel --citation
cd

## Make temp directory 
mkdir ~/project/tmp

# Install PLINK1
cd ~/apps
sudo wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
sudo unzip -o plink_linux_x86_64_20230116.zip
sudo rm plink_linux_x86_64_20230116.zip
## Check install version
~/apps/plink --version #PLINK v1.90b7 64-bit (16 Jan 2023)
cd

# Install PLINK2
cd ~/apps
sudo wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_avx2_20240818.zip
sudo unzip -o plink2_linux_avx2_20240818.zip 
sudo rm plink2_linux_avx2_20240818.zip 
## Check install version
~/apps/plink2 --version #PLINK v2.00a5LM AVX2 Intel (7 Jun 2023)
cd

# Install CMake
sudo yum install -y gcc gcc-c++ gcc-gfortran openssl-devel 
cd ~/apps
wget https://cmake.org/files/v3.18/cmake-3.18.0.tar.gz
tar -xvzf cmake-3.18.0.tar.gz
cd cmake-3.18.0
./bootstrap
gmake
sudo gmake install
cd

# Install R
#### R dependencies ####
sudo yum install -y readline-devel \
  zlib-devel bzip2 bzip2-devel xz xz-devel \
  libcurl libcurl-devel libcurl.i686 libcurl-devel.x86_64 \
  findutils libffi-devel \
  libxml2-devel pcre pcre2 java nlopt nlopt-devel \
  fontconfig-devel.x86_64 fribidi-devel.x86_64 harfbuzz-devel.x86_64 \
  freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel \
  harfbuzz harfbuzz-devel fribidi fribidi-devel fftw fftw-devel
# sudo yum update -y

# Older R's to autosetup some background programs
# R v3
sudo yum remove -y libicu60 libicu60-devel
sudo yum install -y R
# R R v4.0.2
sudo amazon-linux-extras install -y R4


# Install latest R
mkdir -p ~/apps/
sudo chmod 777 -R ~/apps
cd ~/apps/
wget https://cran.r-project.org/src/base/R-4/R-4.4.0.tar.gz
tar xf R-4.4.0.tar.gz
cd R-4.4.0/

./configure --prefix=$HOME/R-4.4.0/ --with-x=no --with-pcre1
make
# Add the following to ~/.bash_profile. 
echo export PATH=~/apps/R-4.4.0/bin:$PATH >> ~/.bash_profile
# Restart 

#Update libicu
sudo yum install -y libicu60 libicu60-devel

# Install R packages
install.packages("tidyverse", Ncpus=10)
install.packages(c("data.table","foreach","doParallel",
                   "BiocManager", "devtools","GGally"),
BiocManager::install(c("SNPRelate","SNPassoc","SeqArray","SeqVarTools",
                       "GENESIS"))
install.packages(c("GMMAT"), Ncpus=10)

# devtools::install_github("BIGslu/kimma")

#### GDAL ####
# sudo yum install -y gdal-devel geos-devel proj-dev
# 
# sudo yum install -y gcc-c++.x86_64 cpp.x86_64 \
#   sqlite-devel.x86_64 libtiff.x86_64 cmake3.x86_64 \
#   proj.x86_64 proj-devel.x86_64 proj-epsg.x86_64 proj-nad.x86_64
# sudo yum install -y proj.x86_64 
# 
# cd ~/apps/
# wget https://download.osgeo.org/proj/proj-6.1.1.tar.gz
# tar -xvf proj-6.1.1.tar.gz
# cd proj-6.1.1
# ./configure
# sudo make
# sudo make install
# 
# cd ~/apps/
# wget https://github.com/OSGeo/gdal/releases/download/v3.2.1/gdal-3.2.1.tar.gz
# tar -xvf gdal-3.2.1.tar.gz
# cd gdal-3.2.1
# ./configure --with-proj=/usr/local --with-python
# sudo make
# sudo make install
# sudo yum update -y
# 
# which gdalinfo; gdalinfo --version /usr/bin/gdalinfo

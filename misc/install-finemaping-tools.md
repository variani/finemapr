## FINEMAP

finemap="finemap_v1.1_x86_64"
appdir="$HOME/.local/apps/finemap/"

mkdir install_finemap
cd install_finemap

# make a binary file 
wget http://www.christianbenner.com/${finemap}.tgz
tar -xvzf ${finemap}.tgz 
cd ${finemap}/

# install in a local directory
mkdir -p ${appdir}
cp ${finemap} ${appdir}/finemap
cp -r example ${appdir}/

## CAVIAR

- https://github.com/fhormoz/caviar
- http://genetics.cs.ucla.edu/caviar/manual.html

```
mkdir install_caviar
cd install_caviar

# compile
git clone https://github.com/fhormoz/caviar.git
cd caviar/CAVIAR-C++/
make

# install in a local directory
mkdir -p ~/apps/caviar
cp *CAVIAR ~/apps/caviar
cp -r sample_data/ ~/apps/caviar/sample_data/

# run example 1
./CAVIAR -c 2 -z sample_data/50_LD.txt -l sample_data/50_Z.txt -o log
# Segmentation fault: 11

# run example 2
./CAVIAR -c 2 -z sample_data/DDB1.top100.sig.SNPs.ZScores -l sample_data/DDB1.top100.sig.SNPs.ld -o log
# output files: log_set, log_post, log.

# run example 3
./CAVIAR -c 2 -z ~/apps/finemap/example/region1.z -l ~/apps/finemap/example/region1.ld -o log
head log_set # prints two snps: rs15 & rs47
```

## PAINTOR

- https://github.com/gkichaev/PAINTOR_V3.0
- running without annotations: https://github.com/gkichaev/PAINTOR_V3.0/issues/11#issuecomment-303135031

```
mkdir install_paintor
cd install_paintor

# compile
git clone https://github.com/gkichaev/PAINTOR_V3.0.git
cd PAINTOR_V3.0
bash install.sh

# create example file with 1 region
mkdir SampleData1
cp SampleData/Locus1 SampleData/Locus1.* SampleData1/
head -n 1 SampleData/input.files > SampleData1/input.files

# install in a local directory
mkdir -p ~/apps/paintor
cp PAINTOR ~/apps/paintor
cp -r CANVIS ~/apps/paintor/CANVIS
cp -r PAINTOR_Utilities ~/apps/paintor/PAINTOR_Utilities
cp -r SampleData ~/apps/paintor/SampleData
cp -r SampleData1 ~/apps/paintor/SampleData1

# run example on Locus 1
./PAINTOR -input SampleData1/input.files -in SampleData1/ -out SampleData1/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS
```




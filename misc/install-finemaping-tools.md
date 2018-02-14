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
cp sample_data/ ~/apps/caviar/sample_data/

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

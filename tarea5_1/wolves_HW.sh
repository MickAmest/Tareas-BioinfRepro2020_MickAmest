#### Script to perform HW test on wolves data using Plink
#### Miguel Amaro, February 2020

## Download data
mkdir ../data
cd ../data
wget https://datadryad.org/stash/downloads/file_stream/6226
mv 6226 wolves.vcf
cd ../bin

## Convert vcf file to plink format using VCFtools
# Create shortcut to run VCFtools on a Docker container and remove volume afterwards
vcftools="docker run -u 1600 --rm -v /home/cirio/miguelamaro/BioinfinvRepro/Unidad5/Prac_Uni5/data:/data biocontainers/vcftools:0.1.15 vcftools"
#Convert vcf file to Plink format
$vcftools --vcf ../data/wolves.vcf --plink --out ../data/wolves

## Run Plink to test for HW equilibrium. Before runing this, make sure to download the Plink and place it in the bin directory
./plink --file ../data/wolves --hardy --allow-extra-chr --chr-set 38 --out ../data/wolves
head ../data/wolves.hwe

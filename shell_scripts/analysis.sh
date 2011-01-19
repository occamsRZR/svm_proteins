#!/bin/sh
### Job name 
#PBS -N msa_phy
### Output files
#PBS -e test.err
#PBS -o test.log
### The variables need to be imported  
#PBS -v num,dir
#PBS -V


cd /home/bfulk/public_html/gPanda_v0.1/temp_dir/$dir/$num
touch gotcha$num
### run the MSA using mafft's linsi
/export/SOFTWARE/mafft-6.240/scripts/mafft $num.fa > $num.aln
### run phylogenic analysis with 100 bootstraps, using our .phy as input
/export/SOFTWARE/FastTree/FastTree $num.aln > $num.tree

exit 0
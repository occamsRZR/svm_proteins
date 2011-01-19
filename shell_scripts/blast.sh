#!/bin/bash
### Job name 
#PBS -N BLSTjorb
### Output files
#PBS -e test.err
#PBS -o test.log
### We need some variables
#PBS -v num,dir
### Work in the current working directory



cd /home/bfulk/public_html/gPanda_v0.1/temp_dir/$dir/
/export/SOFTWARE/ncbi-blast-2.2.24+/bin/blastp -db /export/DATABASES/blast/nr -query temp_blast$num -p blastp -out temp_results$num

exit 0

#!/arch/bin/perl
#
# svm.pl
#
# This program utilizes support vector machines to mine for distantly related proteins.
# 
# Author        : Brandon Fulk
# Date created  : 10/15/2010
# Last modified : 12/14/2010
BEGIN{
    unshift @INC, "./modules";
    unshift @INC, "./modules/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi";
}
use strict;
use warnings;
use Bio::SearchIO;
use IO::File;
use Getopt::Long;
use SVMProteins;
#use SVMProteins::SVMTrainingSet;
#use SVMProteins::SVMTestingSet;

my($pos_infile, $neg_infile, $outfile, $descriptor, $optimize, $help, $gamma, $c, $all_genomes, $genome, $all_proteins);

my $result = GetOptions ("pos_infile=s" 	=> \$pos_infile,
						"neg_infile=s"	=> \$neg_infile,
						"outfile:s"		=> \$outfile,
						"descriptor:s"	=> \$descriptor,
						"optimize"		=> \$optimize,
						"gamma:s"		=> \$gamma,
						"tradeoff:s"	=> \$c,
						"genome:s"      => \$genome,  
						"all_genomes!"	=> \$all_genomes,
						"all_proteins!" => \$all_proteins,
						"help"			=> \$help);
						
						

# If the user sets the help flag or $result returns with an error
# print the usage of this program
if((!$result) || $help || ! (defined $pos_infile && defined $neg_infile) ){
	print "\n
Usage of this program:
	
    Parameters:
	    --pos_infile=   supply the location of your positive data
	    --neg_infile=   supply the location of your negative data
	    --outfile=      supply the location of where you want the positive and negative data to be stored
	    --descriptor=   supply the descriptor you want to use on this SVM
	            aa  :   amino acid composition 
	            dip :   dipeptide composition
	            trip:   tripeptide composition
	    --gamma=        explicitly set the gamma value for the SVM
	    --tradeoff=     explicitly set the tradeoff value (C) for the SVM
	    --genome=       supply the three letter code of which genome you want to mine
	            ath :   Arabidopsis thaliana
	            chl :   Chlamydomonas reinhardtii	    
	    
    Flags:
        --optimize      set the optimize flag if you want this program to optimize the gamma and tradeoff values foryou
	    
	    --help          set this flag and return to this usage screen
	    
	Example:
	    perl svm.pl --pos_infile=infile_pos --neg_infile=infile_neg --outfile=data_set --optimize

";
	exit;
}
###
# Here is the most used command line for this file
# perl svm.pl --pos_infile=infile_pos --neg_infile=infile_neg --outfile=data_set 
#
###

# create our args hash and set the pos and neg datasets
my %args = (pos_location =>  $pos_infile,
            neg_location =>  $neg_infile);
    
# if the outfile is defined, put it in the args hash (default is 'outfile')
if(defined $outfile){
    $args{outfile} = $outfile;
}
# if the descriptor is set, put it in the args hash (default is dipeptide)
if(defined $descriptor){
    $args{descriptor} = $descriptor;
}
# set the default here, or else the method lookup table won't work
else{
    $args{descriptor} = "dip";
}
# if the genome is set, put it in the args hash (default is ath)
if(defined $genome){
    $args{genome}   = $genome;
}
else{
    $args{genome}   = "ath";
}


# if optimize is set, put it in the args hash (default is not to optimize)
if(defined $optimize){
    $args{optimize} = $optimize;
}
# if gamma is set, put it in the args hash (default is 8)
if(defined $gamma){
    $args{gamma} = $gamma;
}
# if tradeoff is set, put it in the args hash (default is 8)
if(defined $c){
    $args{tradeoff} = $c;
}



# create a new SVMProteins object with the args
my $new_svm = SVMProteins->new(%args);
$new_svm->start;

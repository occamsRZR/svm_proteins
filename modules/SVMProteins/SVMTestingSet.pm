package SVMTestingSet;
use base 'SVMProteins';
use strict;
use warnings;
use Algorithm::SVM;
use Bio::SeqIO;
use Bio::SearchIO; 
use Data::Dumper;
use IO::File;
use DBI;
use constant {
	#BLAST_BLOCKS is the number of blocks you want your blast runs split into
	BLAST_BLOCKS	=> 10,	
        TOP_DESC        => 10	    
};



##################################################
# new
#
# DESCRIPTION:  this method is called whenever a new
#               SVMTestingSet object is created.
#               It will accept any %args passed.
#
# IN
#       %arguments  : parameters set in the .pl program
#
# OUT
#       $self       : will return the $self blessed with 
#                   the class
##################################################
sub new {
    my $class = shift;
    my %arguments = @_;
    
    my $self = {
        pos_location => "pos_infile",
        neg_location => "neg_infile",
        outfile      => "outfile",
        descriptor   => "dip",
        optimize     => 0,
        gamma        => 8,
        tradeoff     => 8,
        %arguments,
    };
    # bless this with the class
    bless($self, $class);
    return $self;
}

##################################################
# test
#
# DESCRIPTION:  this method will use all the arguments
#               provided in the constructor and will test
#               the genomes of interest.  a hash will be returned
#               with the sequences >0 SVM value
#
# IN
#       nothing is needed as the input
#
# OUT
#       $svm_hash   : hash of hashes that is returned
#           KEY         VALUE
#           id          id of protein
#           svm_score   score predicted from the svm
#           seq         sequence of protein of interest
#           length      length of sequence of interest
##################################################
sub test {
    my $self        = shift;
    my $trained_svm = shift;
    my $descriptor  = $self->{descriptor};
    my $genome      = $self->{genome};
    my %results     = ();
        
    # DBI variables
    my $platform  = "mysql";
    my $database  = "test";
    my $host      = "bioservweb";
    my $port      = "3306";
    my $tablename = "results";
    my $user      = "bfulk";
    my $password  = "TEc4E3Ud";
# DATA SOURCE NAME                                                                                                         
    my $dsn = "dbi:mysql:$database:$host:$port";

# create a database handle                                                                                                 
my $dbh = DBI->connect($dsn, $user, $password)
    or die ("cannot connect to database");

    # create a new BioSeqIO object with the genome location
    my $genome_obj = Bio::SeqIO->new(   -file   => "$genome",
                                        -format => "fasta");
                                        
    # this will loop through all the sequences in the genome                             
    while(my $in = $genome_obj->next_seq){
        my $description = $descriptor->($in->seq);
        my $refer       = SVMProteins::attribute_ref($description);
        my $test        = Algorithm::SVM::DataSet->new(Label    => 0,
                                                       Data     => $refer);
        my $res         = $trained_svm->predict_value($test);
        if($res > 0){
            my %temp_hash           = ();         
            $temp_hash{id}          = $in->id;
	    my @match = ($in->id =~ /\w{2}\|\d+\|\w+\|(\w{3}\d{5})/);
            $temp_hash{svm_score}   = $res;
            $temp_hash{seq}         = $in->seq;
            $temp_hash{length}      = $in->length;
	    my $length              = $in->length;
            $results{$in->id}       = \%temp_hash;
    # DEFINE a MySQL Query for this entry                                                                                 
                                                                                                                          
    my $query = "INSERT INTO $tablename ( score, length, analysis_id, job_id, gi) VALUES ( $res, $length, 43, 192, '@match')";
	    my $query_handle = $dbh->prepare($query);
	    $query_handle->execute();
        }
    }
    return \%results;
}

##################################################
# analysis
#
# DESCRIPTION:  This method will perform analysis on each of the results.
#               It will take each result, perform a multiple sequence alignment (MSA) 
#				 and subsequent phylogenetic analysis with the known positive proteins.
#
#               
#
# IN
#       $results       : this is a reference to a hash
#                   that contains all the results
#
# OUT
#       nothing is returned
#       
##################################################
sub analysis{
    my $self    	   = shift;
    my $results      = shift;
    my $pos     	   = $self->{pos_location};
    my $dir 		   = $self->{dir};
    
    # open the temp_blast file, generated in rBLAST
    my $seqio_obj   = Bio::SeqIO->new(  -file   => "./temp_dir/$dir/temp_blast.all",
                                        -format => "fasta");
                                        
    my $i = 0;
    # while there are still sequences in the temp_blast file
    while(my $in = $seqio_obj->next_seq){
        # we need to make a temp directory to house each hit
        `mkdir ./temp_dir/$dir/$i`;
        # we need to copy the positive dataset to a new temp file.
		`cp $pos ./temp_dir/$dir/$i/$i.fa`;
        # open the copied file you just created
        open(my $FH, '>>', "./temp_dir/$dir/$i/$i.fa") or die "The temp file in the analysis method could not be opened/appened to. \n $!";
        # add the sequence from the temp_blast file into our positive data
        print $FH ">QUERY$i ID" . $in->id . "\n" . $in->seq;
		$results->{$in->id}->{analysis_num} = $i;
		# submit this job to qsub queue
		`qsub -v num=$i,dir=$dir ./shell_scripts/analysis.sh`;
		
        $i++;
		close($FH) or die "The temp file in the analysis method could not be closed. \n $!";
    }
    return $results;
}

##################################################
# printResults
#
# DESCRIPTION:  This method will take the final results hash reference and print
#               out a results summary into an excel spreadsheet.  This spreadsheet
#               will have the ID, svm_score, length and analysis number (this will show where
#               the folder is that contains the MSA and PHY files) as well as a word
#               count of the reciprocal BLAST
#
# IN
#       $results    : this is a reference to a hash
#                  	  please see test method description for a 
#                        more thorough description
#
# OUT
#       $return     : nothing is returned
##################################################
sub printResults {
    my $self    = shift;
	my $results = shift;
	my $dir     = $self->{dir};
	
	open(my $SUM, '>', "./temp_dir/$dir/results_summary.xls") or die "The summary file could not be opened.";
	
	print $SUM "Summary of Results:\n";
	print $SUM "The results can be found in the temp_dir directory in this folder: $dir \nThe multiple sequence alignments and phylogenies can be found in the respective ANALYSIS NUM folder\n\n";
	print $SUM "ID:\tSVM_SCORE\tLENGTH:\tANALYSIS NUM:\tDESCRIPTION COUNT(Top TOP_DESC):\n";
	
    # for each of the word counts, we need to find the top number default is 10 set by TOP_DESC
    foreach my $id (keys %$results){
        print $SUM "$id \t";
        print $SUM $results->{$id}->{svm_score} . "\t";
        print $SUM $results->{$id}->{length} . "\t";
        print $SUM $results->{$id}->{analysis_num} . "\t";
        
        # we then need to get our desc_count and print out the word and the count
        # this hash isn't too big, so we can fully dereference it.
        my %desc_count  = %{$results->{$id}->{desc_count_table}};
        # for each of the word counts, we need to find the top number default is 10 set by TOP_DESC
		foreach my $word (sort {$desc_count{$b} <=> $desc_count{$a}} keys %desc_count){
			print $SUM "$word:$desc_count{$word}\t";
		}
		print $SUM "\n";
    }
    close($SUM) or die "The summary file could not be closed!";
	return;
}



#################################################
# rBLAST
#
# DESCRIPTION:  This is a reciprocal BLAST method.  It will split the sequences
#               from the hash generated in the test method.  The number of files
#               the sequences are split into is determined by BLAST_BLOCKS.  Each
#				 block is submited to the queue using qsub.  The goal here is to have the 	
#               each block running in parallel to reduce the amount of time the blast will take.           
#               
# IN
#       $results    : this is a reference to a hash
#                  	  please see test method description for a 
#                        more thorough description
#
# OUT
#				nothing is returned
#       		
##################################################
sub rBLAST {
    my $self    = shift;
	my $results = shift;
	my $dir = $self->{dir};
    
	# make a temporary directory to store all the blast and analysis files
    `mkdir ./temp_dir/$dir`;

    # open a file that includes all the sequences
	open(my $TEMP_ALL, '>', "./temp_dir/$dir/temp_blast.all") or die "The temp BLAST containing all the sequences could not be created $!";
    my @ids     = keys %$results;
    my $block = 0;
	# foreach of the IDs, write each id and sequence to separate temp files
    foreach my $id (@ids){

		# open a temp file that is just a block of all the results
		open(my $TEMP, '>>', "./temp_dir/$dir/temp_blast$block") or die "The temp BLAST files could not be created";
        # print the id and seq in this temp block in FASTA format
		print $TEMP ">" . $id;
        print $TEMP "\n" . $results->{$id}->{seq} . "\n";
		
		# print the id and seq in the TEMP_ALL file in FASTA format
		print $TEMP_ALL ">" . $id;
		print $TEMP_ALL "\n" . $results->{$id}->{seq} . "\n";
		
		
		$block++;
		# if this reaches the number of blocks we want, reset $ to 0
		if($block == BLAST_BLOCKS){
			$block = 0;
		}
		close($TEMP) or die "The temp file could not be closed";
    }
	close($TEMP_ALL) or die "The temp file including all the sequences could not be closed";
        
    sleep(30);
	# This for loop will submit each block using a shell script with qsub
    for($block = 0; $block < BLAST_BLOCKS; $block++){
		`qsub -v num=$block,dir=$dir ./shell_scripts/blast.sh`;
	#	# wait for 5 seconds to submit new job
		sleep(5);
	}

    print "Waiting for BLAST to finish...";
	# This block is needed to see if there are still submitted jobs
	my $stat = `qstat`;
	while($stat =~ /BLSTjorb/g){
		$stat = `qstat`;
		print ".";
		# check again in 15 seconds
		sleep (15);
	}

	# After all the jobs have completed...
	# We then need to "stitch" our reports of our blocks together
	for(my $block = 0; $block < BLAST_BLOCKS; $block++){
		`cat ./temp_dir/$dir/temp_results$block >> ./temp_dir/$dir/temp_results.all`;
	}
	
	return $results;
}

##################################################
# parseBLAST
#
# DESCRIPTION:  This is a method to parse each of the BLAST reports. The number
#               of BLAST reports is determined by BLAST_BLOCKS.  Each report is parsed
#				 and the $result hash table is changed to include a word count of the results. 	
#               
# IN
#       $results    : this is a reference to a hash
#                  	  please see test method description for a 
#                        more thorough description
#
# OUT
#		 $results	:	the reference to the hash is returned with the reflected changes
#       		
##################################################
sub parseBLAST{
	 my $self = shift;
	 my $results = shift;
	 my $dir = $self->{dir};
	 my $in = new Bio::SearchIO(-format  => 'blast', 
													-file   => "./temp_dir/$dir/temp_results.all");
    while( my $result = $in->next_result ) {
	## $result is a Bio::Search::Result::ResultI compliant object
		my $query = $result->query_name;
		# create a new hash for each result
		my %desc_count = (); 
        while( my $hit = $result->next_hit ) {
        ## $hit is a Bio::Search::Hit::HitI compliant object


            while( my $hsp = $hit->next_hsp ) {
            ## $hsp is a Bio::Search::HSP::HSPI compliant object
                if( $hsp->length('total') > 20 ) {
                    if ( $hsp->percent_identity >= 75 ) {
                        # create a new hash, a word count for the description
						my $desc = $hit->description();
						my @words = ($desc =~ /(\w+)/g);
						# foreach of the words captured, count them in our hash
						foreach my $word (@words){
							$desc_count{$word}++;
						}
                    }
                }
            }  
   	    }
		my $i = 1;
		my %temp = ();
		# for each of the word counts, we need to find the top number default is 10 set by TOP_DESC
		foreach my $word (sort {$desc_count{$b} <=> $desc_count{$a}} keys %desc_count){
			#print "$word: $desc_count{$word}\n";
			$temp{$word} = $desc_count{$word};
			if($i == TOP_DESC){
				last;
			}
			$i++;
		}
		
		# We then want to find each of the elements in our results hash so we can add the desc_count hash
		my @ids = keys %$results;
		foreach my $id(@ids){
			# if the query name matches the $id, add the temp hash to the results hash
			if($id eq $query){
				$results->{$id}->{desc_count_table} = \%temp ;
			}
		}
    }
	 return($results);
}


1;

package SVMTrainingSet;
use base "SVMProteins";
use strict;
use warnings;
use IO::File;
use Algorithm::SVM;
use Algorithm::SVM::DataSet;
use Data::Dumper;


##################################################
# new
#
# DESCRIPTION:  this method is called whenever a new
#               SVMTrainingSet object is created.
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
# train
#
# DESCRIPTION:  this method will use all the arguments
#               provided in the constructor and will train
#               the dataset.  an SVM object will be returned
#               so it can be used with a SVMTestingSet object
#
# IN
#       nothing is needed as the input
#
# OUT
#       $svm    : svm that is returned 
##################################################
sub train {
    my $self         = shift;
    my $optimize     = $self->{optimize};
    my @training_set = $self->prepare();
    my $tradeoff     = $self->{tradeoff};
    my $gamma        = $self->{gamma};
    
    
    # if the optimize flag was set, run svm_parameters
    if($optimize){
        # this will overwrite the parameters set above
        ( $tradeoff, $gamma ) = $self->svm_parameters();
    }
    # make a new Algorithm::SVM object, train and return
    my $svm = new Algorithm::SVM(Type 	=> 'C-SVC',
                                 Kernel => 'radial');
    
    # we can then set the parameters and train our training set
    $svm->C($tradeoff);
    $svm->gamma($gamma);
    $svm->train(@training_set);

    my $accuracy = $svm->validate(10);

    print "\nThis is the 10-fold validation accuracy:  $accuracy \n";


    my $test_number = 0;
    foreach my $test(@training_set){                                  
        $test_number++;
        my $res = $svm->predict_value($test);
        print "\n\n \t The prediction for this data $test_number is: \t" . $res . "\n";
    }
    
    return $svm;
}


##################################################
# prepare($pos_infile, $neg_infile, $outfile, $descriptor)
#
# DESCRIPTION: uses the data_set made in the prepare
#				to run a python script in order
#  				to find the optimal parameters
#	
# IN	
#       $pos_infile		: this is a string that points to the
#						location that contains the positive data
#		$neg_infile		: this is a string that points to the 
#						location that contains the negative data
#		$outfile		: this is where all the SVM descriptor
#						data will reside (pos and neg datasets)
#		$descriptor		: this is a string which tells which
#						descriptor we need to run
#
# OUT   
#       @training_set	: array containing the training set 
#						
##################################################
sub prepare{
    my $self  = shift;
	my $pos_infile = $self->{pos_location};
	my $neg_infile = $self->{neg_location};
	my $outfile    = $self->{outfile};
	my $descriptor = $self->{descriptor};
	my @training_set = ();
	
	open(my $FH, '>', "$outfile") or die "The outfile could not be written!\n";
	
	#While loop to get the pos data
	my $seqio_obj_pos = Bio::SeqIO->new(  -file => "$pos_infile",
									 	-format => "fasta");
	while(my $in = $seqio_obj_pos->next_seq){
		# This uses the %descriptor_table to get the sub reference
		# and the $descriptor variable to determine which subroutine to call
		my $description = $descriptor->($in->seq);
		print $FH  "+1 " . $description;
		my $refer = SVMProteins::attribute_ref($description);
		my $pos_ds = new Algorithm::SVM::DataSet(Label => +1,
									         Data  => $refer );

		push @training_set, $pos_ds;
	}


	#While loop to get the neg data
	my $seqio_obj_neg = Bio::SeqIO->new(  -file => "$neg_infile",
								 		-format => "fasta");
	while(my $in = $seqio_obj_neg->next_seq){
		# This uses the %descriptor_table to get the sub reference
		# and the $descriptor variable to determine which subroutine to call
		my $description = $descriptor->($in->seq);
		print $FH "-1 " . $description;
		my $refer = SVMProteins::attribute_ref($description);
		my $neg_ds = new Algorithm::SVM::DataSet(Label => -1,
									         Data  => $refer );
		push @training_set, $neg_ds;
	}
	$FH->close;
	
	return @training_set;
}

##################################################
# svm_parameters($outfile)
#
# DESCRIPTION: uses the data_set made in the prepare
#				to run a python script in order
#  				to find the optimal parameters
#	
# IN	
#       $data_set			: this is the outfile that was
#						generated by the prepare() sub
#
# OUT   
#       @para  			: array containing the C 
#						and gamma parameters as well
#						as their accuracy
#						
##################################################
sub svm_parameters{
    my $self    = shift;
	my $data_set = $self->{outfile};	
	#run grid.py and extract the best parameter data
	my $line = `python grid.py $data_set | tail -1`;
	my @para = split (' ', $line);
	print "grid.py determined that the best parameters for this model are: \n \tC: $para[0] \n\tGamma: $para[1] \n This produces an Accuracy of $para[2]\%\n";
	return @para;
}


1;

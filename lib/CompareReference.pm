##### Module to Compare Reference Sequences#####
##### Author: Shalabh Thakur ################
##### Date: 18-July-2013 #####################

#!/usr/bin/perl
package CompareReference;
use strict;
use Exporter;
use FindBin qw($Bin);
use File::Basename;
use Parallel::ForkManager;
use Configuration;
use AdjustSequence;
use SequenceHash;
use Hmmer;
use FilterPair;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(compareReference buildReferenceSequenceDB);


sub compareReference {

    my(%config_param)=%{(shift)};
    my(%out_dir)=%{(shift)};
    my(@reference_file)=@{(shift)};
    my(%seq_len_table)=%{(shift)};
    my(%seq_id_table)=%{(shift)};
    my($inseq_dir)=(shift);
    my($hmmer)=(shift);
    my $hmm_program="phmmer";
    my $mcl_output=undef;
    my $cluster_file=undef;
    my $cluster_num=1;
    my $similarity_file="SimilarityPair.txt";

    #### BUILD REFERENCE DATABASE ####
    print STDOUT "BUILD SEQUENCE DATABASE FOR REFERENCE GENOME\n";
    my ($db_file,$count_ref)=buildReferenceSequenceDB($inseq_dir,\@reference_file,\%out_dir,\%config_param);

    if($config_param{REFERENCE}->{NUMBER_REFERENCE}!=$count_ref){
         print STDOUT "Number of Reference files do not match with number specified in configuration file.\n";
         exit;
    }
    
    open(SIM_FILE,">$out_dir{mcl_dir}/$similarity_file");
    close SIM_FILE;

    ### START COMPARING EACH REFERENCE GENOME ######

    foreach(@reference_file){

           if($_=~/^\.+$/){next;}
           
           my $genome_name=undef;
           my $reference_file=undef;
           my $program=$hmmer."/".$hmm_program;
           
           
           $genome_name=$_;
           $genome_name=~s/\.(\w+)//g;   
           $reference_file=$inseq_dir."/".$_;            

          
           print STDOUT "Start Scanning $reference_file\n";
          
           Hmmer::Run($program,$reference_file,$out_dir{hmm_db_dir}."/".$config_param{DATABASE}->{SEQUENCE_DB},
                      $out_dir{hmm_fullout_dir},
                      $out_dir{hmm_tblout_dir},
                      $out_dir{hmm_domout_dir},
                      $config_param{HMMER});

          
           print STDOUT "Parsing $reference_file scan result\n";

           my($read_phmmer_output,$hit_feature)=Hmmer::Read_domain_table($hmm_program,$reference_file,$out_dir{hmm_domout_dir});

           my($sequence_alignment)=Hmmer::Read_aligned_sequence($hmm_program,$reference_file,$out_dir{hmm_fullout_dir},$hit_feature); 
     
           FilterPair::getHitForQuery($hmm_program,$genome_name,$config_param{REFERENCE},$config_param{PARSE_HMMER},$config_param{HMMER},$out_dir{similarity_dir},$out_dir{all_similarity_dir},$out_dir{chimeric_similarity_dir},$read_phmmer_output,$sequence_alignment,$seq_id_table{$genome_name},$seq_len_table{$genome_name});           
           
           system("cat $out_dir{similarity_dir}/besthit_$genome_name.txt >> $out_dir{mcl_dir}/$similarity_file");                 
    }
    return($similarity_file);
}


##### Subroutin to build sequence database for reference #####

sub buildReferenceSequenceDB {

    my($inseq_dir)=(shift);
    my(@reference_file)=@{(shift)};
    my(%out_dir)=%{(shift)};
    my(%config_param)=%{(shift)};   

       my $db_file=$out_dir{hmm_db_dir}."/$config_param{DATABASE}->{SEQUENCE_DB}";
       open(db_file,">".$db_file); close db_file;
       my $count_ref=0;

       foreach(@reference_file){

           if($_=~/^\.+$/){next;}

           my $file=$inseq_dir."/".$_;
           system("cat $file >> $db_file");
           $count_ref++; 
       }
    return($db_file,$count_ref);  
}


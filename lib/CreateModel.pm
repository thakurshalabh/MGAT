##### Module to create HMM Models from alignment #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package CreateModel;
use strict;
use Exporter;
use File::Basename;
use File::Path qw(remove_tree);
use Hash::Merge qw( merge );
use Bio::AlignIO;
use Bio::SeqIO;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(createHMMModel);

##### Subroutin to create fasta, alignment and hmm model #####
sub createHMMModel {

    my(%model_gene)=%{(shift)};
    my(%group_for_query)=%{(shift)};   
    my(%out_dir)=%{(shift)};
    my($hmmer)=(shift);
    my $process=(shift);    
    my $group_modeled=0;

    #### Initiate parallel processes ####    

    #my $fork=Parallel::ForkManager->new($process);   

     #### Print Profiled Alignment for the Group #####
    foreach my $group_id(keys %model_gene){

            if(!defined($group_for_query{$group_id})){
               next;
            }   
            $group_modeled++;           
            my $out=Bio::AlignIO->new(-file=>">$out_dir{aln_dir}/$group_id.afa",
                                      -format=> "fasta"
                                     );
            my $alnObj=$model_gene{$group_id};  
               $out->write_aln($alnObj);           
     }

      #### Build HMM Model for group alignment ##### 
     foreach my $group_id(keys %model_gene){

            if(!defined($group_for_query{$group_id})){
               next;
            }  
            my $total_gene_in_model=`grep '>' $out_dir{aln_dir}/$group_id.afa | wc -l` ;                       
               chomp($total_gene_in_model);
         
            if($total_gene_in_model==1){      
               system("cp $out_dir{aln_dir}/$group_id.afa $out_dir{singleton_group}/$group_id.fasta");
               system("perl -pi -e 's/^>/>$group_id\t/g' $out_dir{singleton_group}/$group_id.fasta");
               next;
            }else{            
               #$fork->start and next;               
               system("$hmmer/hmmbuild --informat afa $out_dir{hmm_file_dir}/$group_id.hmm $out_dir{aln_dir}/$group_id.afa > $out_dir{tmp_log}/hmmbuild.log");
               unlink("$out_dir{singleton_group}/$group_id.fasta");                          
               #$fork->finish(); 
            }                                
     }

    #$fork->wait_all_children; 
   
    return($group_modeled);      
}

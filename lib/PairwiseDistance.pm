##### Module to calculate Pairwise Distance between protein pairs#####
##### Author: Shalabh Thakur ################
##### Date: 18-AUG-2013 #####################

#!/usr/bin/perl
package PairwiseDistance;
use strict;
use Exporter;
use FindBin qw($Bin);
use File::Basename;
use File::Path qw(remove_tree);
use Hash::Merge qw( merge );
use List::Util qw(sum);
use Bio::Tools::Phylo::Phylip::ProtDist;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(findDistanceForUnModelPair);


##### Subroutin to select gene for creating hmm models #####

sub findDistanceForUnModelPair {

    my(%sequence)=%{(shift)};
    my(%homologue_group)=%{(shift)};
    my(%model_gene)=%{(shift)};
    my(%out_dir)=%{(shift)};
    my($muscle)=(shift);
    my($protdist)=(shift);
    my($process)=(shift);
    my($tmp)=(shift);
    my(%pair_distance)=();      
    my $bin=$Bin;
    
    if(-e "$tmp/forkmanager"){
          remove_tree("$tmp/forkmanager");
     }
     mkdir("$tmp/forkmanager");

     #### Initiate parallel processes ####
     my $fork=Parallel::ForkManager->new($process,"$tmp/forkmanager");  

         $fork->run_on_finish(
                  #### Run this when each child process finish #####
                  sub {                    
                      my($pid, $exit_code, $ident, $exit_signal, $core_dump, $send_parent) = @_;

                      my @send_parent=@{$send_parent};                  
                      my %distance_pair=%{$send_parent[0]};                                    
                      my $group_id=$send_parent[1];                                
                      %pair_distance=%{merge(\%pair_distance,\%distance_pair)};                                                                         
                  }              
         );

         foreach my $group_id(keys %homologue_group){ 
                                                                    
                    my $gene_list=$homologue_group{$group_id};
                    my @gene_list=split("\t",$gene_list);
 
                    if(scalar(@gene_list) eq 1){
                       next;
                    }
                    $fork->start and next;                                                                           
                    my @send_parent=();                
                    my $ref_pair_distance={}; 
                                   
		    ($ref_pair_distance)=GeneticDistanceForPair(\%sequence,\@gene_list,\%model_gene,$group_id,\%out_dir,$muscle,$protdist,$tmp);            		      
         
                    push(@send_parent,$ref_pair_distance);
                    push(@send_parent,$group_id);    	     

                    $fork->finish(0,\@send_parent);
               }
               $fork->wait_all_children;

   return(\%pair_distance);
}     

###### Subroutin to calculate Genetic Distance for whole Group ##############

sub GeneticDistanceForPair{
   
    my(%sequence)=%{(shift)};
    my(@gene_list)=@{(shift)};
    my(%model_gene)=%{(shift)};
    my($group_id)=(shift);   
    my(%out_dir)=%{(shift)};
    my($muscle)=(shift);
    my($protdist)=(shift);
    my($tmp)=(shift);
    my(%pair_distance)=();
    my(%paired_gene)=();
    my $bin=$Bin;    

     mkdir("$out_dir{tmp_log}/$group_id"); 

     ##### PROTDIST PARAMETER FILE ###################
     open(PARAMS,">$out_dir{tmp_log}/$group_id/param");
     print PARAMS "$bin/$out_dir{tmp_log}/$group_id/tmp_pair.phy\n";
     print PARAMS "Y";
     close PARAMS;         

    foreach my $gene_id(@gene_list){
               
         my $seq_counter="0001"; 
         my $total_dist_for_group=0;
         my $count_pair_for_group=0; 
         my $average_dist_for_group=0;         
         my %temp_id_map=();    
               
         $gene_id=~/(\w+)(\|)(\w+)/;          
         my $gene_seq=$sequence{$1}->{$gene_id};             
         $gene_seq=~s/\*$//g; 
         my $temp_gene_id="Seq_".$seq_counter;
         $temp_id_map{$temp_gene_id}=$gene_id;                     
         $seq_counter++;           
       
         if($paired_gene{$group_id}){
 
            my %paired_id=%{$paired_gene{$group_id}}; 
                         
            foreach my $paired_id(keys %paired_id){

               if(!defined($model_gene{$group_id}->{$paired_id}) and ($gene_id ne $paired_id)){                 
                  
                  open(TMP_PAIR,">$out_dir{tmp_log}/$group_id/tmp_pair.fasta");
                  print TMP_PAIR ">$temp_gene_id\n$gene_seq\n";                                  
                  $paired_id=~/(\w+)(\|)(\w+)/;          
                  my $paired_seq=$sequence{$1}->{$paired_id};             
                  $paired_seq=~s/\*$//g; 
                  my $temp_pair_id="Seq_".$seq_counter;
                  $temp_id_map{$temp_pair_id}=$paired_id;
                  print TMP_PAIR ">$temp_pair_id\n$paired_seq\n";                                  
                  close TMP_PAIR;

                  ##### MULTIPLE ALIGNMENT OF Query Gene with Model Proteins ############
                  system("$muscle/muscle -in $out_dir{tmp_log}/$group_id/tmp_pair.fasta -phyiout $out_dir{tmp_log}/$group_id/tmp_pair.phy -log $out_dir{tmp_log}/muscle.log -verbose -quiet");

                  #### PROTDIST PAIRWISE DISTANCE CALCULATION #######
                  chdir("$out_dir{tmp_log}/$group_id");
                  system("$bin/$protdist/protdist < $bin/$out_dir{tmp_log}/$group_id/param > $bin/$out_dir{tmp_log}/protdist.log");
                   
                  #### PARSE PROTDIST OUTPUT FILE ################### 
                  my($pair_distance)=ParseDistanceMatrix("outfile",$group_id,\%temp_id_map);         
                  unlink("outfile");
                  chdir($bin);                            
                  %pair_distance=%{merge(\%pair_distance,$pair_distance)};               
              }

            }            
                      
         }elsif(!defined($paired_gene{$group_id})){
              $paired_gene{$group_id}->{$gene_id}=$gene_id;
              next;
         }
         $paired_gene{$group_id}->{$gene_id}=$gene_id;             
    }  

    remove_tree("$out_dir{tmp_log}/$group_id"); 
  
    return(\%pair_distance);
}

####### Parse Distance Matrix #######################

sub ParseDistanceMatrix {

    my($outfile)=(shift);
    my($group_id)=(shift);
    my(%temp_id_map)=%{(shift)};

    my %group_distance=(); 
    my %pair_distance=();   

    my $dist = Bio::Tools::Phylo::Phylip::ProtDist->new(-file=> $outfile,                  
	                                                -format=> 'phylip',
                                                        -program=>"ProtDist");  
    my $matrix=$dist->next_matrix;

    foreach my $temp_idA(keys %temp_id_map){
          foreach my $temp_idB(keys %temp_id_map){
                  if($temp_id_map{$temp_idA} eq $temp_id_map{$temp_idB}) {next;}                                   
                  my $distance_value = $matrix->get_entry($temp_idA,$temp_idB);                
                  if($group_distance{$temp_idB."_".$temp_idA}){
                     next;
                  }                
                  $group_distance{$temp_idA."_".$temp_idB}=$distance_value;
                  $pair_distance{$group_id}->{$temp_id_map{$temp_idA}}->{$temp_id_map{$temp_idB}}=$distance_value;                              
          }
    }
    return(\%pair_distance);
}
 

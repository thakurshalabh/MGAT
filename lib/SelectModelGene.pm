##### Module to Select Gene for HMM Models#####
##### Author: Shalabh Thakur ################
##### Date: 13-AUG-2013 #####################

#!/usr/bin/perl
package SelectModelGene;
use strict;
use Exporter;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use File::Path qw(remove_tree);
use Hash::Merge qw( merge );
use List::Util qw(sum);
use Bio::Tools::Phylo::Phylip::ProtDist;
use Bio::AlignIO;
use Bio::SeqIO;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(findGeneForHMMModel);


##### Subroutin to select gene for creating hmm models #####

sub findGeneForHMMModel {

    my(%sequence)=%{(shift)};
    my(%homologue_group)=%{(shift)};
    my(%group_for_query)=%{(shift)};
    my(%group_gene)=%{(shift)}; 
    my(%group_distance)=%{(shift)};    
    my(%out_dir)=%{(shift)};
    my($muscle)=(shift);
    my($protdist)=(shift);
    my($process)=(shift);
    my($tmp)=(shift);
    my(%model_gene)=();
    my(%processed_genome)=();      
    my $bin=$Bin;
    
     if(-e "$tmp/forkmanager"){
          remove_tree("$tmp/forkmanager");
     }
     mkdir("$tmp/forkmanager");

     #### Initiate parallel processes ####

=for comment
         my $fork=Parallel::ForkManager->new($process,"$tmp/forkmanager");  
 
         $fork->run_on_finish(
                  #### Run this when each child process finish #####
                  sub {                    
                      my($pid, $exit_code, $ident, $exit_signal, $core_dump, $send_parent) = @_;

                      my @send_parent=@{$send_parent};
                      my %distance=%{$send_parent[0]}; 
                      my $model_aln_object=${$send_parent[1]}; 
                      my %gene_family=%{$send_parent[2]};
                      my %distance_pair=%{$send_parent[3]};                                    
                      my $group_id=$send_parent[4]; 

                      $group_distance{$group_id}=\%distance;
                      $model_gene{$group_id}=$model_aln_object; 
                      %group_gene=%{merge(\%group_gene,\%gene_family)};
                      %pair_distance=%{merge(\%pair_distance,\%distance_pair)};                                                                                                                                    
                  }              
         );
=cut
         foreach my $group_id(keys %group_for_query){
                                                                     
                    my $gene_list=$homologue_group{$group_id};
                    my @gene_list=split("\t",$gene_list);                     
                      
                    if(scalar(@gene_list) eq 1){
                         
                         mkdir("$out_dir{tmp_log}/$group_id");                         

                         my $gene_id=shift(@gene_list);                       
                         $gene_id=~/(\w+)(\|)(\w+)/;                            
                         my $gene_seq=$sequence{$1}->{$gene_id};             
                         $gene_seq=~s/\*$//g;
                        
                         open(ALN_FILE,">$out_dir{tmp_log}/$group_id/$group_id.afa");
                         print ALN_FILE ">$gene_id\n$gene_seq";
                         close ALN_FILE; 

                         my $in=Bio::AlignIO->new(-file=>"$out_dir{tmp_log}/$group_id/$group_id.afa",
                                                  -format=> "fasta"
                                                 );

                         my $alnObj=$in->next_aln();  
                         $model_gene{$group_id}=$alnObj;                                            
                         $group_distance{$group_id}->{total_group_distance}=0;
                         $group_distance{$group_id}->{average_group_distance}=0;
                         $group_distance{$group_id}->{number_pair}=0; 
                         remove_tree("$out_dir{tmp_log}/$group_id");                       
                         next;                                        
                    }         
                    #$fork->start and next; 
                    #my @send_parent=();             
                    #my $ref_pair_distance={};
                    my $ref_group_distance={};
                    my $ref_model_gene={};
                    my $ref_processed_genome={};  
                    my $alnObj=undef; 

                    if(-e "$out_dir{aln_dir}/$group_id.afa"){
                         my $in=Bio::AlignIO->new(-file=>"$out_dir{aln_dir}/$group_id.afa",
                                                  -format=> "fasta"
                                                  );
                         $alnObj=$in->next_aln();
                    }               
                                               
		    ($ref_group_distance,$ref_model_gene,$ref_processed_genome)=GeneticDistanceForGroup(\%group_gene,\%sequence,\@gene_list,\$alnObj,\$group_distance{$group_id},$group_id,\%out_dir,$muscle,$protdist,$tmp);            		                  			
                      
                    $model_gene{$group_id}=${$ref_model_gene};
                    $group_distance{$group_id}=$ref_group_distance;
                    %processed_genome=%{merge(\%processed_genome,$ref_processed_genome)};
                    
                    #push(@send_parent,$ref_group_distance);
                    #push(@send_parent,$ref_model_gene);
                    #push(@send_parent,$ref_group_gene);
                    #push(@send_parent,$ref_pair_distance);
                    #push(@send_parent,$group_id);   	                        
                    #$fork->finish(0,\@send_parent);
        }
        #$fork->wait_all_children;
        %group_gene=%{merge(\%group_gene,\%processed_genome)};

   return(\%model_gene,\%group_gene,\%group_distance);
}     

###### Subroutin to calculate Genetic Distance for whole Group ##############

sub GeneticDistanceForGroup{

    my(%group_gene)=%{(shift)};
    my(%sequence)=%{(shift)};
    my(@gene_list)=@{(shift)};
    my($model_gene_group)=${(shift)};
    my($group_distance_ref)=${(shift)};
    my($group_id)=(shift);   
    my(%out_dir)=%{(shift)};
    my($muscle)=(shift);
    my($protdist)=(shift);
    my($tmp)=(shift);
    my(%processedgenome)=();
    my $bin=$Bin;     

    my %group_distance=();
    
    if(!defined($group_distance_ref)){
       $group_distance{total_group_distance}=0;
       $group_distance{average_group_distance}=0;
       $group_distance{number_pair}=0;       
    }else{
      %group_distance=%{$group_distance_ref};
    }       

    mkdir("$out_dir{tmp_log}/$group_id");     

    foreach my $gene_id(@gene_list){

         $gene_id=~/(\w+)(\|)(\w+)/;

         if($group_gene{$1} or $gene_id eq ''){
              next;
         }         
         my $total_dist_for_group=0;
         my $count_pair_for_group=0; 
         my $average_dist_for_group=0;  
                      
         $gene_id=~/(\w+)(\|)(\w+)/;  
         $processedgenome{$1}=$1;        
         my $gene_seq=$sequence{$1}->{$gene_id};             
         $gene_seq=~s/\*$//g;                       
        
         if(defined($model_gene_group)){

            ###### MODEL ALIGNMENT #######
            my $out=Bio::AlignIO->new(-file=>">$out_dir{tmp_log}/$group_id/group_aln.afa",
                                      -format=> "fasta"
                                     );
            my $alnObj=$model_gene_group;  
               $out->write_aln($alnObj);             
            my $count_aln_seq=`grep '>' $out_dir{tmp_log}/$group_id/group_aln.afa | wc -l` ;
               chomp($count_aln_seq);     

            ##### MULTIPLE ALIGNMENT OF Query Protein with Model Protein Alignment ############  
            if($count_aln_seq>1){              
               open(GENE_SEQ,">$out_dir{tmp_log}/$group_id/gene_seq.afa");
               print GENE_SEQ ">$gene_id\n$gene_seq";
               close GENE_SEQ;       
               system("$muscle/muscle -profile -in1 $out_dir{tmp_log}/$group_id/group_aln.afa -in2 $out_dir{tmp_log}/$group_id/gene_seq.afa -out $out_dir{tmp_log}/$group_id/$group_id.afa -log $out_dir{tmp_log}/muscle.log -verbose -quiet");
            
            }elsif($count_aln_seq eq 1){             
               open(GROUP_ALN,">>$out_dir{tmp_log}/$group_id/group_aln.afa");
               print GROUP_ALN ">$gene_id\n$gene_seq";
               close GROUP_ALN;       
               system("$muscle/muscle -in $out_dir{tmp_log}/$group_id/group_aln.afa -out $out_dir{tmp_log}/$group_id/$group_id.afa -log $out_dir{tmp_log}/muscle.log -verbose -quiet");
            }
   
            ##### PARSE PHYLIP FORMAT #####                   
            my ($phylip_alignment_file,$id_map)=convertFastaToPhylip("$out_dir{tmp_log}/$group_id/$group_id.afa","$out_dir{tmp_log}/$group_id/tmp_aln.afa");

            ##### PROTDIST PARAMETER FILE ###################
            open(PARAMS,">$out_dir{tmp_log}/$group_id/param");
            print PARAMS "$bin/$phylip_alignment_file\n";
            print PARAMS "Y";
            close PARAMS;   

            #### PROTDIST PAIRWISE DISTANCE CALCULATION #######
            chdir("$out_dir{tmp_log}/$group_id");
            system("$bin/$protdist/protdist < $bin/$out_dir{tmp_log}/$group_id/param > $bin/$out_dir{tmp_log}/protdist.log");
                   
            #### PARSE PROTDIST OUTPUT FILE ################### 
            my($dist_array_ref)=ParseDistanceMatrix("outfile",$group_id,$gene_id,$id_map);         
            unlink("outfile");
            chdir($bin);
                                     
            $count_pair_for_group=(scalar(@{$dist_array_ref}));
            $total_dist_for_group=sprintf "%.5f",(sum @{$dist_array_ref});                                            
         }
         elsif(!defined($model_gene_group)){

             open(GROUP_ALN,">$out_dir{tmp_log}/$group_id/$group_id.afa");
             print GROUP_ALN ">$gene_id\n$gene_seq";
             close GROUP_ALN; 
             
             my $in=Bio::AlignIO->new(-file=>"$out_dir{tmp_log}/$group_id/$group_id.afa",
                                      -format=> "fasta"
                                     );

             my $alnObj=$in->next_aln();  
             $model_gene_group=$alnObj;                    
             $group_distance{total_group_distance}=0;
             $group_distance{average_group_distance}=0;
             $group_distance{number_pair}=0; 
             next;                                         
         } 

         #### CALCULATE TOTAL DISTANCE AND AVERAGE DISTANCE FOR PROTEIN_ID and MODEL PROTEINS #####             
         $total_dist_for_group=$group_distance{total_group_distance} + $total_dist_for_group;
         $count_pair_for_group=$group_distance{number_pair} + $count_pair_for_group;           
              
         $average_dist_for_group=($total_dist_for_group/$count_pair_for_group);
         $average_dist_for_group=sprintf "%.5f",$average_dist_for_group;                  

         if($average_dist_for_group > $group_distance{average_group_distance}){                     

               my $in=Bio::AlignIO->new(-file=>"$out_dir{tmp_log}/$group_id/$group_id.afa",
                                        -format=> "fasta"
                                        );
               my $alnObj=$in->next_aln();  
               $model_gene_group=$alnObj;                       
               $group_distance{total_group_distance}=$total_dist_for_group;
               $group_distance{average_group_distance}=$average_dist_for_group;
               $group_distance{number_pair}=$count_pair_for_group;            
         }          
    }  
    remove_tree("$out_dir{tmp_log}/$group_id");   
    return(\%group_distance,\$model_gene_group,\%processedgenome);
}

##### PARSE FASTA TO PHYLIP FORMAT ######

sub convertFastaToPhylip {

    my($fasta_aln_file)=(shift);
    my($tmp_aln_file)=(shift);
    my %id_map=();

    my $in_fasta_aln=Bio::SeqIO->new(-file =>"$fasta_aln_file",-format => "fasta");

    open(TMP_ALN,">$tmp_aln_file");
    my $seq_counter="0001";
  
    while(my $aln_seq=$in_fasta_aln->next_seq()){
        my $seq_id=$aln_seq->id;
           $seq_id=~s/\/(.*)//g;
        my $temp_id="SEQ_".$seq_counter;
        print TMP_ALN ">$temp_id\n",$aln_seq->seq(),"\n";
        $id_map{$temp_id}=$seq_id;
        $seq_counter++; 
    }
    close TMP_ALN;

    my $phylip_file=$fasta_aln_file;
       $phylip_file=~s/\.afa/\.phy/g;

    system("perl Fasta2Phylip.pl $tmp_aln_file $phylip_file");

    return($phylip_file,\%id_map);
}

####### Parse Distance Matrix #######################

sub ParseDistanceMatrix {

    my($outfile)=(shift);
    my($group_id)=(shift);
    my($gene_id)=(shift);
    my(%temp_id_map)=%{(shift)};

    my %group_distance=(); 
    my @dist_array=(); 

    my $dist = Bio::Tools::Phylo::Phylip::ProtDist->new(-file=> $outfile,                  
	                                                -format=> 'phylip',
                                                        -program=>"ProtDist");  
    my $matrix=$dist->next_matrix;

    foreach my $temp_idA(keys %temp_id_map){
                
            if($temp_id_map{$temp_idA} ne $gene_id) {
               next;
            }
          foreach my $temp_idB(keys %temp_id_map){
                  if($temp_id_map{$temp_idA} eq $temp_id_map{$temp_idB}) {next;}                                   
                  my $distance_value = $matrix->get_entry($temp_idA,$temp_idB);                
                  if($group_distance{$temp_idB."_".$temp_idA}){
                     next;
                  }                
                  $group_distance{$temp_idA."_".$temp_idB}=$distance_value;
                  push(@dist_array,$distance_value);               
          }
    }
    return(\@dist_array);
}
 

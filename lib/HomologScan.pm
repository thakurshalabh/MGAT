##### Module to Perform Homolog Prediction Analysis#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package HomologScan;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Hmmer;
use CompareReference;
use FilterPair;
use SelectModelGene;
use PairwiseDistance;
use CreateModel;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use File::Basename;
use File::Copy;
use File::Path qw(remove_tree);
use Parallel::ForkManager;
use Hash::Merge qw( merge );
use MclClustering;
use AddGroupID;
use Initialize;
use Cluster;
use PartialMap;
use Ortholog;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(run_homologue_scan);


sub run_homolog_scan {

##### read input variables ######
my %config_param=%{(shift)};
my %out_dir=%{(shift)};
my %sequence=%{(shift)};
my %seq_len_table=%{(shift)};
my %seq_id_table=%{(shift)};
my $inseq_dir=(shift);
my @reference_file=@{(shift)};
my @query_file=@{(shift)};
my %reference_genome=%{(shift)};
my $hmmer=(shift);
my $mcl=(shift);
my $protdist=(shift);
my $muscle=(shift);
my $process=(shift);
my $tmp=(shift);

##### DECLARE GLOBAL VARAIBLES ####
my %query_genome=();
my %group_gene=();
my %group_distance=();
my %model_gene=();
my $cluster_num=1;
my $cluster_file='';
my $hmm_program='';
my $query_status=0;
my $group_for_query={};
my $jump_to_step=0;
my $count_query=0;

open(LOG_TIME,">$out_dir{tmp_log}/execution_time.log");

##### PHASE I: BUILD SEED / REFERENCE PROTEIN FAMILY USING REFERENCE GENOMES ########


if($config_param{GROUP}->{INPUT_GROUP} eq "NO"){          
      
          my %homologue_group=();
          my $start_reference_scan = time;

          my $similarity_file=CompareReference::compareReference(\%config_param,\%out_dir,\@reference_file,\%seq_len_table,\%seq_id_table,$inseq_dir,$hmmer); 
          ##### MCL CLUSTERING ######          
          my $mcl_output=MclClustering::run_mcl($similarity_file,$mcl,$out_dir{mcl_dir});
          ##### ADD GROUP ID ########
          $hmm_program="phmmer";
          ($cluster_file,$group_for_query)=AddGroupID::add_id($mcl_output,$cluster_num,$hmm_program,\%homologue_group,$out_dir{cluster_dir},$out_dir{mcl_dir});   
          $config_param{GROUP}->{INPUT_GROUP}="YES";
          my $count_groups=keys %{$group_for_query};
          my $end_reference_scan = time;
          my $run_time_reference=$end_reference_scan-$start_reference_scan;
          print LOG_TIME "$cluster_file\t$count_groups\tReference Scan\t$run_time_reference\t",int($run_time_reference/ 3600),"h:",int(($run_time_reference % 3600) / 60),"m:",int($run_time_reference % 60),"s\n";            
} 

##### PHASE II: ##### DENOVO PROTEIN FAMILY PREDICTION USING SEED / REFERENCE PROTEIN FAMILY ######

if($config_param{GROUP}->{INPUT_GROUP} eq "YES"){
      
      my %homologue_group=(); 

      if($cluster_file eq ''){
           print "CHECKING INPUT PARAMTERS.....\n";
           my($group_file,$ref_group_gene,$ref_group_distance,$homologue_group,$ref_model_gene)=Initialize::check_input(\%config_param,\%out_dir,\%sequence); 
           %group_gene=%{$ref_group_gene};
           %group_distance=%{$ref_group_distance}; 
           %homologue_group=%{$homologue_group};
           $group_for_query=$homologue_group;
           %model_gene=%{$ref_model_gene};
           $cluster_file=$group_file;                                      
      }else{
           print "Read Cluster\n";
           my($homologue_group,$query_status)=Cluster::read_cluster($cluster_file,'xxx'); 
           %homologue_group=%{$homologue_group};
           $group_for_query=$homologue_group;     
      }


#       my $count_groups=keys(%{$group_for_query});
#       print "$count_groups\n";
#
#      if((!defined($config_param{GROUP}->{GROUP_ALIGNMENT_DIR})) or (!defined($config_param{GROUP}->{GROUP_DISTANCE_FILE})) or (!defined($config_param{GROUP}->{GROUP_MODEL_DIR})) or (!defined($config_param{GROUP}->{GROUP_SINGLETON_DIR}))){
#            
#            %group_gene=();  
#            ########## Call Group Model Function ######
#            my($ref_group_gene,$ref_group_distance)=GroupModel($cluster_file,$count_groups,\%sequence,\%homologue_group,$group_for_query,\%group_gene,\%group_distance,\%out_dir,$hmmer,$muscle,$protdist,$process,$tmp);
#            
#            %group_gene=%{$ref_group_gene};
#            %group_distance=%{$ref_group_distance};                   
#      }
#      $jump_to_step=4;
       
      foreach(@query_file){

             if($_=~/^\.+$/ or $_=~/\~/){             
                 $count_query++;
                 next;
             }elsif($reference_genome{$_}){               
                  $count_query++;
                  next;
             } 
             my $query_file=$inseq_dir."/".$_;
             my $genome_name=$_;
             $genome_name=~s/\.(\w+)//g;                         
=for comment             
             undef %homologue_group;
             print "READ CLUSTER $cluster_file\n";           
             my($homologue_group,$query_status)=Cluster::read_cluster($cluster_file,$genome_name);
             %homologue_group=%{$homologue_group};

             #### count number of clusters #####
              my $count_groups=keys(%{$group_for_query});   
             print "Number of groups: $count_groups\n";

             ### If any sequence from query file is already present in cluster, do not process that query further #####
           #  if($query_status==1){ 
           #       print "Warning: Cluster file contains sequences from $genome_name, skip this genome\n";
           #       $count_query++;
           #       next;
           #  }   
             if($jump_to_step==4){
                 #$jump_to_step=0;
                 goto(BUILD_DB);                
             }

             ##### CALL Group Model Function ####### 
             my($ref_group_gene,$ref_group_distance)=GroupModel($cluster_file,$count_groups,\%sequence,\%homologue_group,$group_for_query,\%group_gene,\%group_distance,\%out_dir,$hmmer,$muscle,$protdist,$process,$tmp);        
             
             #### UPDATE GLOBAL VARIABLES ##### 
             undef %group_gene;
             undef %group_distance;              

             %group_gene=%{$ref_group_gene};             
             %group_distance=%{$ref_group_distance};

             ###### JUMP HERE IF JUMP_STEP IS EQUAL TO 4 ########
             
             ######### MODEL DATABASE ###########################
             print "BUILD HMM MODEL DATABASE\n";
             if(-e "$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3i"){
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3i");
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3f");
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3m");
                   unlink("$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}.h3p");
             }
             chdir("$out_dir{hmm_file_dir}");
             system("ls | xargs cat > $Bin/$out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}");
             chdir($Bin);

             system("$hmmer/hmmpress $out_dir{hmm_db_dir}/$config_param{DATABASE}->{MODEL_DB}"); 

             ######## SINGLETON DATABASE #########################
             print "BUILD HMM SINGLETON DATABASE\n";
             open(SINGLETON_DB,">$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQUENCE_DB}");
             close SINGLETON_DB;
             
             chdir("$out_dir{singleton_group}");             
             system("ls | xargs cat > $Bin/$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQUENCE_DB}");             
             chdir($Bin);
             system("cat $query_file >> $Bin/$out_dir{hmm_db_dir}/$config_param{DATABASE}->{SEQUENCE_DB}");

             ######## HMMSCAN ANALYSIS ##########################
             print "RUN HMMSCAN for $genome_name\n";
             my $start_hmmscan = time;
             my $program=$hmmer."/"."hmmscan";

             Hmmer::Run($program,$query_file,$out_dir{hmm_db_dir}."/".$config_param{DATABASE}->{MODEL_DB},
                        $out_dir{tmp_log},
                        $out_dir{tmp_log},
                        $out_dir{tmp_log},
                        $config_param{HMMER});  
             
             if(-z "$out_dir{tmp_log}/dom_$genome_name.out" or -z "$out_dir{tmp_log}/full_$genome_name.out"){
               print STDERR "HMMSCAN process failed to generate output for $genome_name\n";
               exit;
             }

             print "Compiling HMM Search Output\n";
             system("cat $out_dir{tmp_log}/dom_$genome_name.out > $out_dir{hmm_domout_dir}/dom_$genome_name.out");
             system("cat $out_dir{tmp_log}/full_$genome_name.out > $out_dir{hmm_fullout_dir}/full_$genome_name.out");
             unlink("$out_dir{tmp_log}/dom_$genome_name.out");
             unlink("$out_dir{tmp_log}/full_$genome_name.out"); 

             my $end_hmmscan = time;
             my $run_hmmscan=$end_hmmscan-$start_hmmscan;
             print LOG_TIME "$cluster_file\thmmscan\t$run_hmmscan\t",int($run_hmmscan/ 3600),"h:",int(($run_hmmscan % 3600) / 60),"m:",int($run_hmmscan % 60),"s\n"; 

             ######## PHMMER ANALYSIS FOR SINGLETON FAMILY ##########
             print "Run Phmmer for Singlton groups\n";
             my $start_phmmer = time; 
             $program=$hmmer."/"."phmmer";     
                               
             Hmmer::Run($program,$query_file,$out_dir{hmm_db_dir}."/".$config_param{DATABASE}->{SEQUENCE_DB},
                        $out_dir{tmp_log},
                        $out_dir{tmp_log},
                        $out_dir{tmp_log},
                        $config_param{HMMER});

             print "Compiling PHMMER Search Output\n";
 
             if(-z "$out_dir{tmp_log}/dom_$genome_name.out" or -z "$out_dir{tmp_log}/full_$genome_name.out"){
               print STDERR "PHMMER process failed to generate output for $genome_name\n";
               exit;
             }
             ##### SHUFFLE COLUMNS IN PHMMER OUTPUT (15: querystart with 17:targetstart) and (16:queryend with 18:targetend) TO MATCH OUTPUT FORMAT FROM HMMSCAN #####
             open(TMP_DOM,"$out_dir{tmp_log}/dom_$genome_name.out");
             my @tmp_dom=<TMP_DOM>;
             close TMP_DOM;
             open(TMP_DOM,">$out_dir{tmp_log}/dom_$genome_name.out");
             
             foreach(@tmp_dom){               
                 if($_=~/\#/) {next;}
                 my @column=split(' ',$_);
                 my $temp=$column[15];
                 $column[15]=$column[17];
                 $column[17]=$temp;
                 $temp=$column[16];
                 $column[16]=$column[18];
                 $column[18]=$temp;

                 foreach my $value(@column){
                    chomp($value);
                    print TMP_DOM "$value\t";
                 }
                 print TMP_DOM "\n";              
              }
              close TMP_DOM;
              undef @tmp_dom;

              system("cat $out_dir{tmp_log}/dom_$genome_name.out >> $out_dir{hmm_domout_dir}/dom_$genome_name.out");
              system("cat $out_dir{tmp_log}/full_$genome_name.out >> $out_dir{hmm_fullout_dir}/full_$genome_name.out");
              unlink("$out_dir{tmp_log}/dom_$genome_name.out");
              unlink("$out_dir{tmp_log}/full_$genome_name.out");

              my $end_phmmer= time;
              my $run_phmmer=$end_phmmer-$start_phmmer;
              print LOG_TIME "$cluster_file\tphmmer\t$run_phmmer\t",int($run_phmmer/ 3600),"h:",int(($run_phmmer % 3600) / 60),"m:",int($run_phmmer % 60),"s\n";
              
              BUILD_DB:
              my $similarity_file="SimilarityPair.txt";
              $hmm_program="hmmscan"; 
              $group_for_query={};

              ####### PARSE OUTPUT FILE FROM HMMSCAN AND PHMMER SEARCH #########
              print STDOUT "PARSING OUTPUT FOR $query_file\n"; 
              Hmmer::Sort_table($query_file,$out_dir{hmm_domout_dir});
              my($read_hmmscan_output,$hit_feature)=Hmmer::Read_domain_table($hmm_program,$query_file,$out_dir{hmm_domout_dir});
              my($sequence_alignment)=Hmmer::Read_aligned_sequence($hmm_program,$query_file,$out_dir{hmm_fullout_dir},$hit_feature);    
              FilterPair::getHitForQuery($hmm_program,$genome_name,$config_param{REFERENCE},$config_param{PARSE_HMMER},$config_param{HMMER},$out_dir{similarity_dir},$out_dir{all_similarity_dir},$out_dir{chimeric_similarity_dir},$read_hmmscan_output,$sequence_alignment,$seq_id_table{$genome_name},$seq_len_table{$genome_name});       
              system("cat $out_dir{similarity_dir}/besthit_$genome_name.txt > $out_dir{mcl_dir}/$similarity_file");

              ####### MCL CLUSTERING ########
              $cluster_num++;
              my $mcl_output=MclClustering::run_mcl($similarity_file,$mcl,$out_dir{mcl_dir});
              ($cluster_file,$group_for_query)=AddGroupID::add_id($mcl_output,$cluster_num,$hmm_program,\%homologue_group,$out_dir{cluster_dir},$out_dir{mcl_dir});
              $count_query++;
=cut
              $query_genome{$genome_name}=$genome_name;            
      }   
      ###### MODEL BUILDING AFTER PROCESSING LAST QUERY IN QUEUE #####
#      print "Computing Cluster and Distance for analysis\n";   
#  
#      my($homologue_group,$query_status)=Cluster::read_cluster($cluster_file,'xxx');
#      %homologue_group=%{$homologue_group};
#      $count_groups=keys(%{$group_for_query}); 
#
#      ##### CALL GROUP MODEL FUNCTION ########   
#      my($ref_group_gene,$ref_model_gene,$ref_group_distance)=GroupModel($cluster_file,$count_groups,\%sequence,\%homologue_group,$group_for_query,\%group_gene,\%group_distance,\%out_dir,$hmmer,$muscle,$protdist,$process,$tmp);
#      
#      system("cp $out_dir{cluster_dir}/Cluster_$cluster_num.txt $out_dir{seed_dir}/model_cluster.txt");
#   
#      $cluster_file="$out_dir{seed_dir}/model_cluster.txt";
#     
#      undef %group_gene;
#      undef %group_distance;      
}

##### MAP PARTIAL SEQUENCE GROUPS ON HOMOLOG GROUPS ####  
print "MAP TRUNCATED GROUPS\n";    

#my($homolog_cluster_file)=PartialMap::mapPartialSequence($cluster_file,\%out_dir,\%seq_len_table,$mcl);

##### PREDICT ORTHOLOG GROUPS ##########################
print "Predict Ortholog Pairs\n";

my $homolog_cluster_file="$out_dir{homolog_dir}/homolog_cluster.txt";

Ortholog::predict_ortholog($homolog_cluster_file,$cluster_file,$muscle,$protdist,\%out_dir,\%sequence,\%query_genome,\%reference_genome);
      
}


###### GROUP MODEL FUNCTION ######
sub GroupModel {
    
    my($cluster_file)=(shift);
    my($count_groups)=(shift);
    my($sequence)=(shift);
    my($homologue_group)=(shift);
    my($group_for_query)=(shift);
    my($group_gene)=(shift);    
    my($group_distance)=(shift);
    my(%out_dir)=%{(shift)};
    my($hmmer)=(shift);
    my($muscle)=(shift);
    my($protdist)=(shift);
    my($process)=(shift);
    my($tmp)=(shift);

     print "Select Proteins for HMM MODEL\n";
                       
     my $start_protein_selection = time;
     my($ref_model_gene,$ref_group_gene,$ref_group_distance)=SelectModelGene::findGeneForHMMModel($sequence,$homologue_group,$group_for_query,$group_gene,$group_distance,\%out_dir,$muscle,$protdist,$process,$tmp);        
     my $end_protein_selection = time;
     my $run_time_divergence=$end_protein_selection-$start_protein_selection;            
     print LOG_TIME "$cluster_file\t$count_groups\tFilter_Diverged_Sequence\t$run_time_divergence\t",int($run_time_divergence/ 3600),"h:",int(($run_time_divergence % 3600) / 60),"m:",int($run_time_divergence % 60),"s\n";
     
     print "Build Protein Family HMM Models\n";

     my $start_model_build = time;    
     my($group_modeled)=CreateModel::createHMMModel($ref_model_gene,$group_for_query,\%out_dir,$hmmer,$process); 
     my $end_model_build = time;
     my $run_time2=$end_model_build-$start_model_build;
     print LOG_TIME "$cluster_file\t$group_modeled\tBuild_Model\t$run_time2\t",int($run_time2/ 3600),"h:",int(($run_time2 % 3600) / 60),"m:",int($run_time2 % 60),"s\n";

     print "GROUP DISTANCE AND PAIRWISE DISTANCE \n";

     printGroupDistance($ref_group_distance,$out_dir{distance_file}); 

    # printPairDistance($ref_pair_distance,$out_dir{distance_file}); 

     return($ref_group_gene,$ref_group_distance);
}

#### PRINT GROUP DISTANCE #####
sub printGroupDistance {

    my(%group_distance)=%{(shift)};
    my($distance_dir)=(shift); 

    open(DISTANCE,">$distance_dir/GroupDistance.txt");
    foreach my $group_id(keys %group_distance){
         print DISTANCE "$group_id\t$group_distance{$group_id}->{total_group_distance}\t$group_distance{$group_id}->{average_group_distance}\t$group_distance{$group_id}->{number_pair}\n";
    }
    close DISTANCE;
}

#### PRINT PAIR DISTANCE ####
=for comment
sub printPairDistance{

    my(%pair_distance)=%{(shift)};
    my($output_dir)=(shift);

    open(PAIR_DISTANCE,">>$output_dir/pairwise_distance.txt");

    foreach my $group_id(sort {$a cmp $b} keys %pair_distance){

            my $ref_geneA=$pair_distance{$group_id};
            my %geneA_hash=%{$ref_geneA};

        foreach my $geneA_id(keys %geneA_hash){
            my $ref_geneB=$geneA_hash{$geneA_id};
            my %geneB_hash=%{$ref_geneB};

            foreach my $geneB_id(keys %geneB_hash){
                my $distance_value=$geneB_hash{$geneB_id};
                print PAIR_DISTANCE "$group_id\t$geneA_id\t$geneB_id\t$distance_value\n";
            }

        }          

    }
    close PAIR_DISTANCE;
}
=cut



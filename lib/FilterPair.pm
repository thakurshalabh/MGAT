##### Module to Filter BEST/PARTIAL/CHIMERIC hits for Query sequences based on threshold #####
##### Author: Shalabh Thakur ################
##### Date: 16-Aug-2013 #####################

#!/usr/bin/perl
package FilterPair;
use strict;
use Exporter;
use File::Basename;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(getHitForQuery);

sub getHitForQuery {

    my($hmm_program)=(shift);
    my($genome_name)=(shift);
    my(%Num_Reference)=%{(shift)};
    my(%PARSE_HMMER)=%{(shift)};
    my(%PARAM_HMMER)=%{(shift)};
    my($best_output_dir)=(shift);
    my($all_output_dir)=(shift);
    my($chimeric_output_dir)=(shift);
    my(%DomainTable)=%{(shift)};
    my(%SequenceAlignment)=%{(shift)};
    my(%SequenceID)=%{(shift)};
    my(%SequenceLen)=%{(shift)};
    
    my @prev_hsp=();
    my @domain=();    
    my @query_coordinate=();
    my @subject_coordinate=();
    my @start_coord=();
    my @end_coord=();
    my $counter=1;    
    my %PairRelation=();
    my %MatchedQuery=();
    my %QueryEvalue=();    
    my @chimera=();
    my %chimeratargetperquery=();
    my %chimeratarget=();
    my %chimeraquery=();
    my $totalidentity=0.0;
    my $totalsimilarity=0.0;
    my $total_alignment_length=0;
    my $file_size=keys(%DomainTable);

    my $best_output_file=$best_output_dir."/besthit_$genome_name.txt";
    my $all_output_file=$all_output_dir."/allhit_$genome_name.txt";
    my $chimeric_output_file=$chimeric_output_dir."/chimeric_$genome_name.txt";
    
    open(BESTHITFILE,">$best_output_file");
    open(ALLHITFILE,">$all_output_file");
    open(CHIMERAFILE,">$chimeric_output_file");

    foreach(sort{$a<=>$b}keys %DomainTable){
             
             my %feature=%{$DomainTable{$_}};           

             my $query_id=$feature{query_name};
             my $subject_id=$feature{target_name};
             my $accuracy=$feature{accuracy};
             my $seq_evalue=$feature{evalue};
             my $c_evalue=$feature{c_evalue};
             my $i_evalue=$feature{i_evalue};
             my $bit_score=$feature{bit_score};
             my $query_len=$feature{query_length};
             my $subject_len=$feature{target_length}; 
             my $qstart=$feature{query_start};
             my $qend=$feature{query_end};
             my $sstart=$feature{target_start};
             my $send=$feature{target_end};
             my $num_domain=$feature{num_domain};
             my $total_domain=$feature{total_domain};
             my $identical_position=$SequenceAlignment{$query_id}->{$subject_id}->{$num_domain}->{identical};
             my $similar_position=$SequenceAlignment{$query_id}->{$subject_id}->{$num_domain}->{similar};
             my $alignment_length=length($SequenceAlignment{$query_id}->{$subject_id}->{$num_domain}->{query_align});
             my $identity=0.0;
             my $similarity=0.0; 

             if($query_id eq $subject_id){
               next; 
             }                
             ##### Calculate Identtity and Similarity of each domain #####
             if($identical_position ne '' and $similar_position eq ''){
                $identical_position=0;
                $similar_position=0;               
             }                     

             $identity=sprintf "%.2f",(($identical_position/$alignment_length)*100);
             $similarity=sprintf "%.2f",(($similar_position/$alignment_length)*100);                  
          
             ##### CHECK ACCURACY, EVALUE IDENTITY AND SIMILARITY Threshold FOR HSP #####          
             if(($accuracy>=$PARSE_HMMER{ACCURACY}) and ($i_evalue <= $PARAM_HMMER{E}) and (($identity>=$PARSE_HMMER{IDENTITY}) or ($similarity>=$PARSE_HMMER{SIMILARITY}))){
                  
             ##### NEXT HSP MATCH FOR QUERY#####
                  if((defined($prev_hsp[1])) and (defined($prev_hsp[0])) and ($subject_id eq $prev_hsp[1]) and ($query_id eq  $prev_hsp[0])){                          
                      my @query_range=([$qstart,$qend]);
                      my @subject_range=([$sstart,$send]);
                      $totalidentity += $identity * $alignment_length;
                      $totalsimilarity += $similarity * $alignment_length;
                      $total_alignment_length += $alignment_length;

                      push(@domain,"$identity\t$similarity\t$num_domain\t$qstart\t$qend\t$sstart\t$send");          
                      push(@query_coordinate,@query_range);
                      push(@subject_coordinate,@subject_range);                                 
           
                   }else{                       
             ##### BEFORE CHECKING NEW QUERY SEQUENCE, PARSE PREVIOUS QUERY #####

                       if(($query_coordinate[0]) and ($subject_coordinate[0]) and ($counter>1)){                          
                          ##### Calculate Query and HMM Coverage #####
                          my $q_coverage=sprintf "%.2f",Coverage(\@query_coordinate,$prev_hsp[2]);
                          my $s_coverage=sprintf "%.2f",Coverage(\@subject_coordinate,$prev_hsp[3]); 

                          my $percent_identity=int($totalidentity/$total_alignment_length * 10 + .5)/10;
                          my $percent_similarity=int($totalsimilarity/$total_alignment_length * 10 + .5)/10;

                          push(@prev_hsp,$q_coverage);
                          push(@prev_hsp,$s_coverage); 
                          push(@prev_hsp,\@domain);
                          push(@prev_hsp,$percent_identity);
                          push(@prev_hsp,$percent_similarity);                         

                          #print "PER_IDENT: $percent_identity\tPER_SIM:$percent_similarity\n";

                          ##### Find Pair Relation between query and target #####
                          my ($relation)=find_QH_Relation(\@prev_hsp,\@query_coordinate,\@subject_coordinate,\%PARSE_HMMER);                           
                          push(@prev_hsp,$relation);
                          $PairRelation{$prev_hsp[0]}->{$prev_hsp[1]}=\@prev_hsp;                                                                  
                       }           
                      
                       #### Print Out Similarity information for query-target pair ######

                       my ($match_id,$best_evalue,$match,$chimera_match,$countchimeratargetperquery,$chimeratarget,$chimeraquery)=PrintOutPairSimilairty(\%PairRelation,\%MatchedQuery,\%QueryEvalue,$hmm_program,\%PARSE_HMMER,\%chimeratargetperquery,\%chimeraquery,\%chimeratarget); 
                          
                       if($match_id ne ''){
                          $MatchedQuery{$match_id}=$match;
                          $QueryEvalue{$match_id}=$best_evalue;

                       }elsif($chimera_match){
                          push(@chimera,$chimera_match);
                       } 
                        %chimeratargetperquery=%{$countchimeratargetperquery};
                        %chimeraquery=%{$chimeraquery};
                        %chimeratarget=%{$chimeratarget};                   
                          
                       #### REINITIALIZE VARIABLES FOR NEW QUERY SEQUENCE #####
                                         
                        my @query_range=([$qstart,$qend]);                                             
                        my @subject_range=([$sstart,$send]);

                        %PairRelation=();
                        @query_coordinate=();    
                        @subject_coordinate=();
                        @prev_hsp=();
                        @domain=();
                        $totalidentity=0.0;
                        $totalsimilarity=0.0;
                        $total_alignment_length=0;                                         
                        push(@query_coordinate,@query_range);
                        push(@subject_coordinate,@subject_range); 
                   
                        push(@prev_hsp,$query_id);
                        push(@prev_hsp,$subject_id);                      
                        push(@prev_hsp,$query_len); 
                        push(@prev_hsp,$subject_len);
                        push(@prev_hsp,$seq_evalue);
                        push(@prev_hsp,$total_domain);
                        push(@prev_hsp,$accuracy);
                        push(@prev_hsp,$bit_score);
                        $totalidentity += $identity * $alignment_length;
                        $totalsimilarity += $similarity * $alignment_length;
                        $total_alignment_length += $alignment_length;
                        push(@domain,"$identity\t$similarity\t$num_domain\t$qstart\t$qend\t$sstart\t$send");                                                                   
                   }
             } 

         ##### IF last line of file ######
         if($counter==$file_size){ 
           ##### Calculate Query and HMM Coverage #####
           if(($accuracy>=$PARSE_HMMER{ACCURACY}) and ($i_evalue<=$PARAM_HMMER{E}) and ($identity>=$PARSE_HMMER{IDENTITY}) or ($similarity>=$PARSE_HMMER{SIMILARITY})){                                       
              
              my $q_coverage=Coverage(\@query_coordinate,$prev_hsp[2]);        
              my $s_coverage=Coverage(\@subject_coordinate,$prev_hsp[3]);  

              my $percent_identity=int($totalidentity/$total_alignment_length * 10 + .5)/10;
              my $percent_similarity=int($totalsimilarity/$total_alignment_length * 10 + .5)/10; 
     
              push(@prev_hsp,$q_coverage);
              push(@prev_hsp,$s_coverage); 
              push(@prev_hsp,\@domain);
              push(@prev_hsp,$percent_identity);
              push(@prev_hsp,$percent_similarity);                        

              ##### Find Pair Relation #####
              my ($relation)=find_QH_Relation(\@prev_hsp,\@query_coordinate,\@subject_coordinate,\%PARSE_HMMER);             
              push(@prev_hsp,$relation);
              $PairRelation{$prev_hsp[0]}->{$prev_hsp[1]}=\@prev_hsp;               
              
              #### Print Out Similarity information for pair ######
              my ($match_id,$best_evalue,$match,$chimera_match,$countchimeratargetperquery,$chimeratarget,$chimeraquery)=PrintOutPairSimilairty(\%PairRelation,\%MatchedQuery,\%QueryEvalue,$hmm_program,\%PARSE_HMMER,\%chimeratargetperquery,\%chimeraquery,\%chimeratarget); 
                 
                  if($match_id ne ''){
                    $MatchedQuery{$match_id}=$match;
                    $QueryEvalue{$match_id}=$best_evalue;                 
                  }elsif($chimera_match){
                    push(@chimera,$chimera_match);
                  }
                  %chimeratargetperquery=%{$countchimeratargetperquery};
                  %chimeraquery=%{$chimeraquery};
                  %chimeratarget=%{$chimeratarget};                  
            }
         }
         $counter++; 
    } 

    foreach my $seq_id(keys %SequenceID){         
         unless(defined($MatchedQuery{$seq_id})){ 
             print BESTHITFILE "$seq_id\t$seq_id\t$SequenceLen{$seq_id}\t$SequenceLen{$seq_id}\t1\t1\t1\t$SequenceLen{$seq_id}\t1\t$SequenceLen{$seq_id}\t0\t100\t100.00\t100.00\t100.00\t100.00\t100.00\t100.00\tSINGLETON\t-\n"; 
             print ALLHITFILE "$seq_id\t$seq_id\t$SequenceLen{$seq_id}\t$SequenceLen{$seq_id}\t1\t1\t1\t$SequenceLen{$seq_id}\t1\t$SequenceLen{$seq_id}\t0\t100\t100.00\t100.00\t100.00\t100.00\t100.00\t100.00\tSINGLETON\t-\n";        
         }elsif(defined($MatchedQuery{$seq_id}) and ($MatchedQuery{$seq_id} ne '') and $hmm_program eq "hmmscan"){
             print BESTHITFILE "$MatchedQuery{$seq_id}";
             print ALLHITFILE "$MatchedQuery{$seq_id}";
         }        
    }

    ###### PRINT CHIMERA ####### 
    foreach my $chimeric_hit(@chimera){
            chomp($chimeric_hit);
            my @chimeric_hit=split("\t",$chimeric_hit);  

            if($chimeraquery{$chimeric_hit[0]} eq 1 and $chimeratarget{$chimeric_hit[1]}<=2){
              print CHIMERAFILE "$chimeric_hit\n";
            }            
    }

    close BESTHITFILE;
    close ALLHITFILE;
    close CHIMERAFILE;  
}


##### Finds percent coverage for the query or subject sequence ####
sub Coverage {

    my(@coordinate)=@{(shift)};
    my($len)=(shift);

    my $prev_start=0;
    my $prev_end=0;
    my $coverage_region=0;
    my $range=0;     

    foreach(sort{$a->[0] <=> $b->[0]} @coordinate){          
             
           my $start=$_->[0];
           my $end=$_->[1];              
         
            if($start>$prev_end){          
               $range=($end-$start)+1;
            }elsif($start<=$prev_end){              
               $range=(($end-$start)+1)-(($prev_end-$start)+1);
            }
            $coverage_region=$coverage_region + $range;
    }
    #### Calculate Percent Coverage #####   
    my $percent_coverage=($coverage_region/$len)*100;
    return $percent_coverage;
}

##### Finds Relation between Query and Hit Sequence #######

sub find_QH_Relation {
    
    my(@hsp_feature)=@{(shift)};
    my(@query_coordinate)=@{(shift)};
    my(@subject_coordinate)=@{(shift)};
    my(%PARSE_HMMER)=%{(shift)};  

    my($query_id)=$hsp_feature[0];
    my($subject_id)=$hsp_feature[1];
    my($query_coverage)=$hsp_feature[8];
    my($subject_coverage)=$hsp_feature[9];
    my($query_len)=$hsp_feature[2];
    my($subject_len)=$hsp_feature[3];
    my($evalue)=$hsp_feature[4];
    my($total_domain)=$hsp_feature[5];
    my($accuracy)=$hsp_feature[6]; 
    my($percent_identity)=$hsp_feature[11];
    my($percent_similarity)=$hsp_feature[12]; 
    my $relation='';
 
    ###### Check for Query-Subject Relation ######

    if($query_coverage>=$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage>=$PARSE_HMMER{HMM_COVERAGE} and (($percent_identity>=$PARSE_HMMER{IDENTITY}) and ($percent_similarity>=$PARSE_HMMER{SIMILARITY}))){
       $relation="BEST";
    }elsif((($query_coverage>=$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage<$PARSE_HMMER{HMM_COVERAGE}) or ($query_coverage<$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage>=$PARSE_HMMER{HMM_COVERAGE})) and (($percent_identity>=$PARSE_HMMER{IDENTITY}) and ($percent_similarity>=$PARSE_HMMER{SIMILARITY}))){
       $relation="TRUNCATED";    
    }elsif(($query_coverage<$PARSE_HMMER{QUERY_COVERAGE} and $subject_coverage<$PARSE_HMMER{HMM_COVERAGE}) and ($query_coverage>=$PARSE_HMMER{MIN_QUERY_COVERAGE} and $subject_coverage>=$PARSE_HMMER{MIN_HMM_COVERAGE}) and (($percent_similarity>=$PARSE_HMMER{SIMILARITY} and $percent_identity>=$PARSE_HMMER{MIN_CHIMERA_IDENTITY}) or ($percent_similarity>=$PARSE_HMMER{MIN_CHIMERA_SIMILARITY} and $percent_identity>=$PARSE_HMMER{IDENTITY}))){       
       ##### Check for Chimera #### 
       my $isChimera=checkForChimericGene($query_id,$subject_id,$evalue,\@query_coordinate,\@subject_coordinate,$query_len,$subject_len,$total_domain);      
       if($isChimera and $accuracy>=$PARSE_HMMER{CHIMERA_ACCURACY}){            
           $relation=$isChimera;        
       }#else{
       #    $relation="PARTIAL";
       #}      
    }    
    return($relation); 
}

##### Check for Chimera Gene ########

sub checkForChimericGene {

    my($query_id)=(shift);
    my($subject_id)=(shift);
    my($evalue)=(shift);   
    my(@query_coordinate)=@{(shift)};
    my(@subject_coordinate)=@{(shift)};
    my($query_len)=(shift);
    my($subject_len)=(shift);
    my($total_domain)=(shift);   
    my $isChimera='';

    ##### look for start and end position of query and subject #####
    my $start_query=$query_coordinate[0]->[0];
    my $end_query=$query_coordinate[0]->[1];

    my $start_subject=$subject_coordinate[0]->[0];
    my $end_subject=$subject_coordinate[0]->[1];
      
      if($total_domain==1){
	  if(($start_query<=30) and ($end_query<$query_len)){
             if(($start_subject<=30) and ($end_subject<$subject_len)){
                 $isChimera="N-Chimera";                 
             }
	  }elsif(($start_query>30) and  ((($query_len-$end_query)+1)<=15)){
             if(($start_subject>30) and ((($subject_len-$end_subject)+1)<=15)){
                $isChimera="C-Chimera";               
             }
	  }       
      } 
   return($isChimera);    
}

##### PrintOut Pairwise Similairty Table #####

sub PrintOutPairSimilairty {

    my(%pair_similarity)=%{(shift)}; 
    my(%MatchedQuery)=%{(shift)};
    my(%QueryEvalue)=%{(shift)};
    my($hmm_program)=(shift);
    my(%PARSE_HMMER)=%{(shift)};
    my(%chimeratargetperquery)=%{(shift)};
    my(%chimeraquery)=%{(shift)};
    my(%chimeratarget)=%{(shift)};
  
    my $match_id=''; 
    my $match=''; 
    my $best_evalue=0.0;
    my $chimera_match='';    

    foreach my $query_id(keys %pair_similarity){

          my %subject_key=%{$pair_similarity{$query_id}};
          $query_id=~/(\w+)(\|)(\w+)/;
          my $query_genome=$1;                

          foreach my $subject_id(keys %subject_key){

                  my @similarity=@{$subject_key{$subject_id}};                 
 
                  my $qlen=$similarity[2];
                  my $slen=$similarity[3];
                  my $evalue=$similarity[4];
                  my $bit_score=$similarity[7];
                  my $qcoverage=sprintf "%.2f",$similarity[8];
                  my $scoverage=sprintf "%.2f",$similarity[9];
                  my @coordinate=@{$similarity[10]};
                  my $percent_identity=$similarity[11];
                  my $percent_similarity=$similarity[12];
                  my $relation=$similarity[13];
                  my $total_domain=$similarity[5];
                  my $target_genome='';

                  if($subject_id=~/Group/){
                     $target_genome='';                     
                  }else{
                     $subject_id=~/(\w+)(\|)(.+)/;
                     $target_genome=$1;       
                  }                  

                  if($query_genome eq $target_genome){                   
                    $hmm_program="phmmer";
                  }elsif($subject_id=~/Group/){                    
                    $hmm_program="hmmscan";
                  }
                

                  foreach my $hsp_coord(@coordinate){

                      my @hsp_coord=split("\t",$hsp_coord);
                      my $identity=shift(@hsp_coord);
                      my $similarity=shift(@hsp_coord);
                         $hsp_coord=join("\t",@hsp_coord);

                      if($relation=~/(Chimera)/){
                         if($percent_similarity>=$PARSE_HMMER{SIMILARITY} or $percent_identity>=$PARSE_HMMER{IDENTITY}){
                            
                            if(!defined($chimeratarget{$subject_id})){
                               $chimeratarget{$subject_id}=0;
                            }
                            if(!defined($chimeraquery{$query_id})){
                               $chimeraquery{$query_id}=1;
                            }                             

                            if($chimeratargetperquery{$query_id}->{$target_genome}){ 
                                $chimeratargetperquery{$query_id}->{$target_genome}=$chimeratargetperquery{$query_id}->{$target_genome} + 1; 
                            }else{
                                $chimeratargetperquery{$query_id}->{$target_genome}=1;
                            }
                            if($chimeratarget{$subject_id}){
                                 $chimeratarget{$subject_id}=$chimeratarget{$subject_id} + 1;
                            }else{
                                 $chimeratarget{$subject_id}=1;
                            }                           
                            
                            if($chimeratargetperquery{$query_id}->{$target_genome}<=2 and $chimeraquery{$query_id} ne 0){
                               $chimeraquery{$query_id}=1;
                            }else{
                               $chimeraquery{$query_id}=0;
                            }           
                            $chimera_match="$query_id\t$subject_id\t$qlen\t$slen\t$total_domain\t$hsp_coord\t$evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$identity\t$similarity\t$qcoverage\t$scoverage\tPARTIAL\t$relation";                       
                         }                         
                      }                      
                      elsif($relation ne ''){
                         if($percent_similarity>=$PARSE_HMMER{SIMILARITY} and $percent_identity>=$PARSE_HMMER{IDENTITY} and $relation eq 'BEST'){

                              if($hmm_program eq "hmmscan" and (!defined($MatchedQuery{$query_id}) or  $evalue<$QueryEvalue{$query_id})){

                                  $match="$query_id\t$subject_id\t$qlen\t$slen\t$total_domain\t$hsp_coord\t$evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$identity\t$similarity\t$qcoverage\t$scoverage\t$relation\t-\n"; 
                                  $match_id=$query_id;
                                  $best_evalue=$evalue; 

                              }elsif($hmm_program eq "phmmer"){
                                  $match=''; 
                                  $match_id=$query_id;
                                  $best_evalue=$evalue;
                                  
                                  print BESTHITFILE "$query_id\t$subject_id\t$qlen\t$slen\t$total_domain\t$hsp_coord\t$evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$identity\t$similarity\t$qcoverage\t$scoverage\t$relation\t-\n";
                                  print ALLHITFILE "$query_id\t$subject_id\t$qlen\t$slen\t$total_domain\t$hsp_coord\t$evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$identity\t$similarity\t$qcoverage\t$scoverage\t$relation\t-\n"; 
                              }                                                                                                        
                         }else{ ### print truncated pairs ####                              
                              print ALLHITFILE "$query_id\t$subject_id\t$qlen\t$slen\t$total_domain\t$hsp_coord\t$evalue\t$bit_score\t$percent_identity\t$percent_similarity\t$identity\t$similarity\t$qcoverage\t$scoverage\t$relation\t-\n";
                              
                         }

                        # $chimeraquery{$query_id}=0;                        
                      }
                  }                                
          }
    }
 
    return($match_id,$best_evalue,$match,$chimera_match,\%chimeratargetperquery,\%chimeratarget,\%chimeraquery);
}

##### Module to access Hmmer package #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package Hmmer;
use strict;
use Exporter;
use File::Basename;

use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(Run Sort_table Read_domain_table Read_aligned_sequence);


sub Run {
    
    my($hmm_program)=(shift);
    my($hmm_file)=(shift);
    my($DB_file)=(shift);
    my($HMMFULL_OUT)=(shift);
    my($HMMTBL_OUT)=(shift);
    my($HMMDOM_OUT)=(shift);    
    my(%HMMER_Param)=%{(shift)};
    my $Param='';

    my $outfile=basename($hmm_file);
       $outfile=~s/(\.)(\w+)/\.out/g;   

    foreach(keys %HMMER_Param){

       if($_=~/^(o)$/){
          $Param=$Param." -".$_." ".$HMMFULL_OUT."/"."full_".$outfile;
       }elsif($_=~/^(tblout)$/){
          $Param=$Param." --".$_." ".$HMMTBL_OUT."/"."tbl_".$outfile;
       }elsif($_=~/^(domtblout)$/){
          $Param=$Param." --".$_." ".$HMMDOM_OUT."/"."dom_".$outfile;
       }elsif($_=~/^([E|A|T|Z])$/){
          $Param=$Param." -".$_." ".$HMMER_Param{$_};
       }else{
          $Param=$Param." --".$_." ".$HMMER_Param{$_};
       }
    }   
    if($hmm_program=~/phmmer/){
      system("$hmm_program $Param $hmm_file $DB_file");
    }elsif($hmm_program=~/hmmscan/){
      system("$hmm_program $Param $DB_file $hmm_file");
    }     
}

sub Sort_table {

    my($query_file)=(shift);
    my($HMMDOM_OUT)=(shift);
     
    my $file_line=1;
    my $file=basename($query_file);
       $file=~s/(\.)(\w+)/\.out/g;  
       $file="dom_".$file; 

    open(DOM_TAB,"$HMMDOM_OUT/$file") or die "Cannot open tabular output file\n";
    my @dom_result=<DOM_TAB>;
    close DOM_TAB;

    my $title_line='';
    my @result_row=();

    foreach(@dom_result){

         if($_=~/\#/ and $file_line<=3) {
             $title_line=$title_line.$_;
             $file_line++;         
             next;
          }
          elsif($_=~/\#/){
            next;
          }
          else{         
             my @dom_line=split(' ',$_);             
             push(@result_row,\@dom_line);
          }
    }  

    my @sorted_result=sort {$a->[3] cmp $b->[3] ||
                            $a->[6] <=> $b->[6]                            
                           } @result_row;

   open(DOM_TAB,">$HMMDOM_OUT/$file") or die "Cannot open tabular output file\n";

   print DOM_TAB "$title_line";

   foreach my $row(@sorted_result){
      foreach my $column(@{$row}){
         print DOM_TAB "$column\t";
      }      
      print DOM_TAB "\n";
  }
}

###### READ HMMER DOMAIN OUTPUT FILE ########

sub Read_domain_table {

    my($hmm_program)=(shift);
    my($query_file)=(shift);
    my($HMMDOM_OUT)=(shift); 
    my(%DomainTable)=();
    my(%HitFeature)=();
      
    my $file_line=1;
    my $file=basename($query_file);
       $file=~s/(\.)(\w+)/\.out/g;  
       $file="dom_".$file;

    #### Read Tabular HMM Dom Output ###

    open(DOM_TAB,"$HMMDOM_OUT/$file") or die "Cannot open tabular $hmm_program output file\n";    

     foreach(<DOM_TAB>){

         if($_=~/\#/) {$file_line++; next; }

         else{         
             my @dom_line=split(' ',$_);
         # target name(0) accession(1) tlen(2) query_name(3) accession(4) qlen(5) E-value(6) score(7) bias(8) #(9) of(10) c-Evalue(11) i-Evalue(12) score(13) bias(14) h_from(15) h_to(16) a_from(17) a_to(18) e_from(19) e_to(20) acc(21) description of target(22)
         #### Read Parameters in scalar variables #### 
            
             $DomainTable{$file_line}->{evalue}=$dom_line[6];
             $DomainTable{$file_line}->{bit_score}=$dom_line[7];
             $DomainTable{$file_line}->{bias}=$dom_line[8];
             $DomainTable{$file_line}->{num_domain}=$dom_line[9];
             $DomainTable{$file_line}->{total_domain}=$dom_line[10];
             $DomainTable{$file_line}->{c_evalue}=$dom_line[11];
             $DomainTable{$file_line}->{i_evalue}=$dom_line[12];
             $DomainTable{$file_line}->{domain_score}=$dom_line[13];
             $DomainTable{$file_line}->{domain_bias}=$dom_line[14];                   
             $DomainTable{$file_line}->{accuracy}=$dom_line[21];
             $DomainTable{$file_line}->{target_description}=$dom_line[22];
             $DomainTable{$file_line}->{envelope_start}=$dom_line[19];
             $DomainTable{$file_line}->{envelope_end}=$dom_line[20];
             
             if($hmm_program eq "phmmer"){
                 $DomainTable{$file_line}->{target_name}=$dom_line[0];
                 $DomainTable{$file_line}->{target_accession}=$dom_line[1];
                 $DomainTable{$file_line}->{target_length}=$dom_line[2]; 
                 $DomainTable{$file_line}->{query_name}=$dom_line[3];
                 $DomainTable{$file_line}->{query_accession}=$dom_line[4];
                 $DomainTable{$file_line}->{query_length}=$dom_line[5];
                 $DomainTable{$file_line}->{query_start}=$dom_line[15];
                 $DomainTable{$file_line}->{query_end}=$dom_line[16];
                 $DomainTable{$file_line}->{target_start}=$dom_line[17];
                 $DomainTable{$file_line}->{target_end}=$dom_line[18];                
             }
             elsif($hmm_program eq "hmmscan"){

                $DomainTable{$file_line}->{target_name}=$dom_line[0];
                $DomainTable{$file_line}->{target_accession}=$dom_line[1];
                $DomainTable{$file_line}->{target_length}=$dom_line[2]; 
                $DomainTable{$file_line}->{query_name}=$dom_line[3];
                $DomainTable{$file_line}->{query_accession}=$dom_line[4];
                $DomainTable{$file_line}->{query_length}=$dom_line[5];                
                $DomainTable{$file_line}->{query_start}=$dom_line[17];
                $DomainTable{$file_line}->{query_end}=$dom_line[18];
                $DomainTable{$file_line}->{target_start}=$dom_line[15];
                $DomainTable{$file_line}->{target_end}=$dom_line[16];
             }
             $HitFeature{$dom_line[3]}=$dom_line[5];
             $HitFeature{$dom_line[0]}=$dom_line[2];                
             $file_line++;
         }
      }      
      close DOM_TAB;
      return(\%DomainTable,\%HitFeature);
}

##### Read Full Output for HMMER Programs ####

sub Read_aligned_sequence {

    my($hmm_program)=(shift);
    my($query_file)=(shift);
    my($HMMFULL_OUT)=(shift);
    my(%HitFeature)=%{(shift)};
 
    my %SequenceAlignment=();

    my $file_line=1;
    my $start_aln=0;
    my $query_seq='';
    my $target_seq='';
    my $sim_line='';
    my $prob_line='';
    my $query_id='';
    my $target_id='';
    my $domain_num=0;

    my $file=basename($query_file);
       $file=~s/(\.)(\w+)/\.out/g;  
       $file="full_".$file;

    #### Read Tabular HMM Dom Output ###
    open(FULL_TAB,"$HMMFULL_OUT/$file") or die "Cannot open full $hmm_program output file\n";    

    foreach(<FULL_TAB>){

          if($_=~/phmmer/){
             $hmm_program="phmmer";
          }

          if($_=~/(\=\=)(\s+)(domain)(\s+)(\d+)/){

            if($query_seq ne '' and $target_seq ne '' and $sim_line ne ''){
                 my $len=0;
                 my $identity=0.0;
                 my $similarity=0.0;

                 if($HitFeature{$query_id}<=$HitFeature{$target_id}) 
                 { 
                    $len=$HitFeature{$query_id};
                 }elsif($HitFeature{$query_id}>$HitFeature{$target_id}){
                    $len=$HitFeature{$target_id};
                 }
 
                 my $identical_region=$sim_line;
                    $identical_region=~s/[\s|\+]//g;
                 my $num_identical=length($identical_region);
                    $identity=sprintf "%.2f",(($num_identical/$len)*100);

                 my $similar_region=$sim_line;
                    $similar_region=~s/\s//g;
                 my $num_similar=length($similar_region);
                    $similarity=sprintf "%.2f",(($num_similar/$len)*100);                

                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_align}=$query_seq; 
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_align}=$target_seq;
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{match_sequence}=$similar_region;          
            }

            $_=~/(\=\=)(\s+)(domain)(\s+)(\d+)/;

            $start_aln=1;
            $query_seq='';
            $target_seq='';
            $sim_line='';
            $prob_line='';
            $query_id='';
            $target_id='';
            $domain_num=$5;
            next;    
          }

          if($_=~/^(\s+)$/ and $start_aln!=2){
            next;
          }

          if($_=~/>>/ or $_=~/(Internal pipeline statistics summary:)/){
            $start_aln=0; next;
          }

          if($start_aln==1){
            my $aln_line=$_;
            $aln_line=~s/^(\s+)//g;
            my @aln_line=split(' ',$aln_line);
               if($hmm_program eq "phmmer"){
                   $query_id=shift(@aln_line);
                   my $seq=$aln_line[1];            
                   $query_seq=$query_seq.$seq;
                }
                elsif($hmm_program eq "hmmscan"){
                   $target_id=shift(@aln_line);
                   my $seq=$aln_line[1];            
                   $target_seq=$target_seq.$seq;
                }
            $start_aln=2;
            next;
          }

          if($start_aln==2){
            my $aln_line=$_;
            $aln_line=~s/^(\s+)//g;                    
            $sim_line=$sim_line.$aln_line;
            $start_aln=3;
            next;
          }

          if($start_aln==3){
            my $aln_line=$_;
            $aln_line=~s/^(\s+)//g;
            my @aln_line=split(' ',$aln_line);

                if($hmm_program eq "phmmer"){
                   $target_id=shift(@aln_line);
                   my $seq=$aln_line[1];            
                   $target_seq=$target_seq.$seq;
                }
                elsif($hmm_program eq "hmmscan"){
                   $query_id=shift(@aln_line);
                   my $seq=$aln_line[1];            
                   $query_seq=$query_seq.$seq;
                }
            $start_aln=4;
            next;
          }
 
          if($start_aln==4){
            my $aln_line=$_;
            $aln_line=~s/^(\s+)//g;                    
            $prob_line=$prob_line.$aln_line;
            $start_aln=1;
            next;
          }
    }

    ###### IF last hit in full alignment file ##########
    if(eof(FULL_TAB)){

         if($query_seq ne '' and $target_seq ne '' and $sim_line ne ''){
                 my $len=0;
                 my $identity=0.0;
                 my $similarity=0.0;

                 if($HitFeature{$query_id}<=$HitFeature{$target_id}) 
                 { 
                    $len=$HitFeature{$query_id};
                 }elsif($HitFeature{$query_id}>$HitFeature{$target_id}){
                    $len=$HitFeature{$target_id};
                 }
 
                 my $identical_region=$sim_line;
                    $identical_region=~s/[\s|\+]//g;                 
                 my $num_identical=length($identical_region); 
             
                 my $similar_region=$sim_line;
                    $similar_region=~s/\s//g;
                 my $num_similar=length($similar_region);                         

                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{identical}=$num_identical;
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{similar}=$num_similar;
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{query_align}=$query_seq; 
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{target_align}=$target_seq;
                 $SequenceAlignment{$query_id}->{$target_id}->{$domain_num}->{match_sequence}=$similar_region;          
           }
    }
    close FULL_TAB;
    return(\%SequenceAlignment);
}






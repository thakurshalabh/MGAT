##### Module to read protein family cluster#####
##### Author: Shalabh Thakur ################
##### Date: 6-AUG-2013 #####################

#!/usr/bin/perl
package Cluster;

use strict;
use warnings;
use Env;
use FindBin qw($Bin);
use lib "$Bin/../lib";


use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(read_cluster readGroupForGene getGroupLength);

####### Genes for each groups ######

sub read_cluster {

    my($cluster_file)=(shift);
    my($genome_name)=(shift);
    my(%homologue_group)=();
    my $genome_status=0;

    open(GROUP_OUT,"$cluster_file") or die "Cannot open the cluster file $cluster_file\n";

    foreach(<GROUP_OUT>){

          chomp($_);

          if($_=~/$genome_name/){

             $genome_status=1;
             last;
          }      
          $_=~/^(\w+)(\:)(\t)(.+)/;
          my $group_line=$4;
          my $group_id=$1;
          $group_line=~s/\t$//g;   
          $homologue_group{$group_id}=$group_line;            
    }

    return(\%homologue_group,$genome_status);
}

##### Group for each gene ####

sub readGroupForGene {

    my($cluster_file)=(shift);
    my %groupforgene=();
 
    open(GROUP_OUT,"$cluster_file") or die "Cannot open the cluster file $cluster_file\n";

    foreach(<GROUP_OUT>){
          chomp($_);
          $_=~/^(\w+)(\:)(\t)(.+)/;

          my $group_id=$1;
          my $group_line=$4;
             $group_line=~s/\t$//g;
          my @genes=split("\t",$group_line);

          foreach my $gene(@genes){
             chomp($gene);
              $groupforgene{$gene}=$group_id;                      
          }                            
    }

   return(\%groupforgene);
}

#### Group Length #####
sub getGroupLength {

    my(%homolog_group)=%{(shift)};
    my(%seq_len_table)=%{(shift)};
    my(%grouplength)=();    
   
    foreach my $group_id(keys %homolog_group){

          my @gene_list=split(' ',$homolog_group{$group_id});

          foreach(@gene_list){              
               $_=~/(\w+)(\|)(\w+)/;
               my $genome_name=$1;
               my $len=$seq_len_table{$1}->{$_};            
               
               if((!defined($grouplength{$group_id})) or ($len>=$grouplength{$group_id})){
                 $grouplength{$group_id}=$len;                               
               }
          }
                
    }
   return(\%grouplength);
}

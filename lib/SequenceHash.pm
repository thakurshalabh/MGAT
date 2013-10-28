##### Module to get sequence ID in hash table #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package SequenceHash;
use strict;
use Exporter;
use File::Basename;
use File::Copy;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(getSequenceFeature);


sub getSequenceFeature {

    my($input_seq_dir)=(shift);
   
    opendir(SEQ_DIR,$input_seq_dir);
    my(@SequenceFile)=readdir(SEQ_DIR); 

    my %SeqLenTable=();
    my %IDHashTable=();
    my %Sequence=();

    foreach(@SequenceFile){ 

       if($_=~/^\.+$/ or $_=~/\~/){next;}
 
       my $file_path=$input_seq_dir."/".$_;

       my $genome_name=$_;
          $genome_name=~s/(\.)(\w+)$//g;
       
       open(FILE,$file_path) or die "Error: File cannot be open\n";    
       my $seq_id='';
       my $seq='';

       # print $file_path,"\t",scalar(@seq_file),"\n";       

       while(<FILE>){                     
            if($_=~/^>/){
                if($seq_id ne ''){

                   $IDHashTable{$genome_name}->{$seq_id}=$seq_id;
                   $Sequence{$genome_name}->{$seq_id}=$seq;
                   $SeqLenTable{$genome_name}->{$seq_id}=length($seq);

                   $seq_id='';
                   $seq='';
                 }                                        
                 $_=~/(^>)(\w+)(\|)(\w+)/;                
                 $seq_id=$2.$3.$4;                
             }else{
               $_=~s/\s//g;              
               $seq=$seq.$_;
             } 

             if(eof(FILE)){                             
                $IDHashTable{$genome_name}->{$seq_id}=$seq_id;
                $Sequence{$genome_name}->{$seq_id}=$seq;
                $SeqLenTable{$genome_name}->{$seq_id}=length($seq);                
             }                                 
       }

      

       close(FILE); 

    }
    return(\%SeqLenTable,\%IDHashTable,\%Sequence);
}

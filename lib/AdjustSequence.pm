##### Module to adjust sequence file #####
##### Author: Shalabh Thakur ################
##### Date: 21-May-2013 #####################

#!/usr/bin/perl
package AdjustSequence;
use strict;
use Exporter;
use File::Basename;
use File::Copy;
use vars qw(@ISA @EXPORT @EXPORT_OK);

@ISA   = qw(Exporter);
@EXPORT= ();
@EXPORT_OK = qw(checkSequence);



sub checkSequence {

    my($Sequence_dir)=(shift);
    my($Sequence_format)=(shift);
    my($Adjust_sequence)=(shift);
    my($input_seq_dir)=(shift);
    my(@SequenceFile)=@{(shift)};    

              
    if($Adjust_sequence eq "YES"){              
             
           foreach(@SequenceFile){                                     
                    my $file_path1=$Sequence_dir."/".$_;
                    my $file_path2=$input_seq_dir."/".$_;
                    my($name)=basename($file_path1);

               if($name=~/\.fasta/ or $name=~/\.fa/ or $name=~/\.fas/){
                    
                    $name=~s/\.(\w+)//g;              
                    open(FILE1,$file_path1) or die "Error:File not found\n";
                    open(FILE2,">".$file_path2) or die "Error: File cannot be created\n";                  

                   while(<FILE1>){                     
                       if($_=~/^>/){                         
                           $_=~s/^>/>$name|/;
                           $_=~s/[\.\-]/\_/g;
                       }                     
                       print FILE2 "$_";                      
                   }
                   close FILE1;
                   close FILE2;
               }
          }
           
    }elsif($Adjust_sequence eq "NO"){

         foreach(@SequenceFile){
              my $file_path1=$Sequence_dir."/".$_;
              my $file_path2=$input_seq_dir."/".$_;
              my($name)=basename($file_path1);

              if($name=~/\.fasta/ or $name=~/\.fa/ or $name=~/\.fas/){
                 copy($file_path1,$file_path2) or die "cannot copy file\n";
              }
         }
    }
}

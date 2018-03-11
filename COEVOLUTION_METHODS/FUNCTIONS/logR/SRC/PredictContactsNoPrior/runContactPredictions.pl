use strict;
use warnings;


if(@ARGV!=1){
  print STDERR "usage: perl runContactPredictions.pl alignment_file(in multifasta format)\n";
  exit(1);
}

my $threshold=0.2; #columns with more than this percentage of gaps will be removed from the alignment

my ($file) = @ARGV;
my $aln_file=$file;
my $posterior_file=$file;
my $ungapped_aln_file=$file;
my $mapping_file=$file;
my $final_posterior_file=$file;
$aln_file=~s/fasta/aln/;
$posterior_file=~s/\.fasta/\_ungapped\.post/;
$ungapped_aln_file=~s/\.fasta/\_ungapped\.aln/;
$mapping_file=~s/fasta/mapping/;
$final_posterior_file=~s/\.fasta/\.post/;
system("perl change_alignment_format.pl $file > $aln_file");
system("perl remove_columns_with_many_gaps.pl $aln_file $ungapped_aln_file $mapping_file $threshold");
system("./calculate_posteriors $ungapped_aln_file > $posterior_file");
system("perl map_positions_back.pl $mapping_file $posterior_file > $final_posterior_file");
 

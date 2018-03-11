use strict;
use warnings;


if(@ARGV!=4){
    print STDERR "usage: alignment_infile alignment_outfile pos_mapping_file gap_threshold\n";
    exit(1);
}

my ($file, $outfile, $mapping_file, $threshold) = @ARGV;


open(IN, $file) || die;
my @m;
my $count=0;
my %frac_gaps;
my $line=<IN>;
my @l=split(/\s+/,$line);
while(<IN>){
    chomp;
    my @s=split(//, $_);
    @{$m[$count]}=@s;
    for(my $i=0;$i<@s;$i++){
	if(!defined $frac_gaps{$i}){
	    $frac_gaps{$i}=0;
	}
	if($s[$i] eq "-"){
	    $frac_gaps{$i}++;
	}
	
    }
    $count++;
}
my $num_cols=0;
foreach my $key (keys %frac_gaps){
    $frac_gaps{$key}/=$count;
    if($frac_gaps{$key}<$threshold){
	$num_cols++;
    }
}
close(IN);
open(OUT2, ">$mapping_file") || die;
open(OUT, ">$outfile") || die;
print OUT scalar(@m), " ", $num_cols, "\n";
for(my $i=0;$i<@m;$i++){
    my $count=0;
    for(my $j=0;$j<@{$m[$i]};$j++){
	if($frac_gaps{$j}<$threshold){
	    if($i==0){
		print OUT2 $j, "\t", $count, "\n";
		$count++;
	    }
	    print OUT $m[$i][$j];
	}
    }
    print OUT "\n";
}
close(OUT);



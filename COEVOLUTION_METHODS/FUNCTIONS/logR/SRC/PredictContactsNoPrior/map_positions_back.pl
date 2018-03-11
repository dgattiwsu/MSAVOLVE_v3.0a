use strict;
use warnings;

if(@ARGV!=2){
    print STDERR "usage: positional_mapping_file posterior_file\n";
    exit(1);
}

my ($pos_file, $post_file) = @ARGV;

open(IN, $pos_file) || die;
my %mapping;
while(<IN>){
    my @s=split;
    $mapping{$s[1]}=$s[0];
}
close(IN);


my %post;
open(IN, $post_file) || die;
while(<IN>){
    my @s=split;
    print $mapping{$s[0]}, " ", $mapping{$s[1]}, " ", $s[2], "\n";
}
close(IN);



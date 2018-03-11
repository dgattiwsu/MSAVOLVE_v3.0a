use strict;
use warnings;

if(@ARGV!=1){
    print STDERR "usage: file\n";
    exit(1);
}

my ($file) = @ARGV;


open(IN, $file) || die;
my %seq;
my $prot;
while(<IN>){
    if($_=~/>(\S+)/){
	$prot=$1;
    }
    else{
	chomp;
	$seq{$prot}.=$_;
    }
}
close(IN);

my @p=keys %seq;
print scalar(@p), "\t", length($seq{$p[0]}), "\n";
foreach my $protein (keys %seq){
    $seq{$protein}=~s/\./\-/g;
    print "$seq{$protein}\n";
}


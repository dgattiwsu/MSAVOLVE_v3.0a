function [ sumCOVdif ] = testrec( glob_COV,glob_COV_orig,recomb_COV,recomb_COV_orig )
 
%   This function checks the coevolution matrices after a recombination cycle.

COVdif = (glob_COV - glob_COV_orig) - (recomb_COV - recomb_COV_orig);
sumCOVdif = sum(COVdif(1:end));


end


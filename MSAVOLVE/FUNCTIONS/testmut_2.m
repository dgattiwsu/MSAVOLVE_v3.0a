function [ sumCOVdif ] = testmut_2( glob_COV,glob_COV_orig,mutcov_COV,...
    mutcov_COV_orig )
 
%  This function checks the coevolution matrices after a mutation cycle.

COVdif = (glob_COV - glob_COV_orig) - (mutcov_COV - mutcov_COV_orig);
sumCOVdif = sum(COVdif(1:end));


end


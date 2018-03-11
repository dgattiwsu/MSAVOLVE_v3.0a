function [ sumCOVdif ] = testmut( COV,COV_orig,mut_COV,mut_COV_orig,...
    cov_COV,cov_COV_orig )

% This function checks the coevolution matrices after a mutation cycle.

COVdif = (COV - COV_orig) - (mut_COV - mut_COV_orig ... 
    + cov_COV - cov_COV_orig);
sumCOVdif = sum(COVdif(1:end));


end


function [ Q_zeros_entropy,Q_ones_entropy,G_entropy ] = binmsa_to_Q_entropy( binmsa )
% This function calculates the global entropy of a binary msa from the
% partition function.

[nrows,ncols] = size(binmsa);

% Here we calculate the global entropy from the partition function using 
% the 0's.
% Either on the columns...

%    E = (nrows-sum(binmsa,1))./nrows;

% Or the rows. If we use the rows, we can apply weights derived from the
% similarity between sequences, for example as it would be calculated 
% using the function 'nmsa_to_weights. 
% W = nmsa_to_weights(MSA_test,threshold,nrows,ncols);

    E = (ncols-sum(binmsa,2))./ncols;
            
%    Q_zeros_entropy = log(sum(exp(-E)./W));
% Where W = 1/m and m is the number of similar sequences within the
% threshold of similarity chosen

    Q_zeros_entropy = log(sum(exp(-E)));
    
    
% Here we calculate the global entropy from the partition function using 
% the 1's.
% Either on the columns...

%    E = sum(binmsa,1)./nrows;

% Or the rows. If we use the rows, we can apply weights.
    
    E = sum(binmsa,2)./ncols;
            
%    Q_ones_entropy = log(sum(exp(-E).*W));

    Q_ones_entropy = log(sum(exp(-E)));
    
% Here we calculate the global entropy in the traditional way either from
% Matlab function:

%    G_entropy = entropy(binmsa);

% or from our own:

    G_entropy = BNMSA_to_gH(binmsa);
    
end


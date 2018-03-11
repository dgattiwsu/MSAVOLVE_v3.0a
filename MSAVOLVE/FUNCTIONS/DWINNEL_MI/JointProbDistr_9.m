 function [ alphabet ] = JointProbDistr_9(nmsa)

% This function calculates the joint probability 
% distribution for the 20 aa's in multiple columns of an msa. 

alphabet = unique(nmsa,'rows');
nrows_alphabet = size(alphabet,1);
[nrows_nmsa,ncols_nmsa] = size(nmsa);

% Loop through the alphabet.

for i = 1:nrows_alphabet
    for j = 1:nrows_nmsa
    ind(j,1) = sum(nmsa(j,:) == alphabet(i,:));
    end    
    ind2 = sum(ind(:,1) == ncols_nmsa);
    frequencies(i,1)= ind2/nrows_nmsa;
end
alphabet = [alphabet frequencies];
end





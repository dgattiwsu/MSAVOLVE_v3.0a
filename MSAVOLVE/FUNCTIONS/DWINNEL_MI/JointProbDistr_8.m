 function [ alphabet ] = JointProbDistr_8(X)

% This function calculates the joint probability 
% distribution for the 20 aa's in two columns of a msa. 

alphabet = unique(X,'rows');
ncols_D = size(alphabet,1);

% Loop through the alphabet.
% Here use "size" instead of "length". "length" is ambiguous if there is
% only one row.

nrows_X = size(X,1);

% Preallocate the logical array D for speed.

D = false(nrows_X,ncols_D);

for i = 1:ncols_D
    
% The following is only for testing.
% for i = 43:43;

% A = X == alphabet(i,1);
% B = X == alphabet(i,2);
% C = [A(:,1) B(:,2)];
A = X(:,1) == alphabet(i,1);
B = X(:,2) == alphabet(i,2);
D(:,i) = A == 1 & B == 1;

% Calculate now the length of the nmsa_D array.
nrows_D = sum(D(:,i));

% Here we calculate the frequency of each pair.
alphabet(i,3) = nrows_D/nrows_X;

end





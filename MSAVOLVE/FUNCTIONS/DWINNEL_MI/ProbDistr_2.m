function [ alphabet ] = ProbDistr_2(X,alphabet)

% Here we calculate the probability distribution
% for the all the symbols (including gaps) in a single column of a msa.
% This script supersedes ProbDistr.m.

unique_alphabet = unique(X);

% Loop through the unique observed alphabet

for i = 1:size(unique_alphabet,1)
    
A = X == unique_alphabet(i,1);
    
unique_alphabet(i,2) = sum(A(:,1))/size(X,1);

% Now replace the unique values of frequencies in the imported complete 
% alphabet.

alphabet(unique_alphabet(i,1),2) = unique_alphabet(i,2);

end









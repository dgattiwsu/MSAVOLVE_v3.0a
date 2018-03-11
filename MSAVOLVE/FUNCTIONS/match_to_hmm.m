function [ hmm_vec,hmm_vec_unif_distr ] = match_to_hmm( input_mat,null_vec )
% This function derives the HMM vector (=the height of a HMM logo) from
% the emission probabilities calculated by the HMMER 3.0 program. We 
% calculate the vector both based on the average distribution of aa and 
% on a uniform distribution (p=1/20 for every aa).
 
[rows,cols]=size(input_mat);
hmm_vec=zeros(rows,1);
hmm_row=zeros(1,cols);
hmm_vec_unif_distr=zeros(rows,1);
hmm_row_unif_distr=zeros(1,cols);

for i=1:rows
for j=1:cols
   
     % Calculate the relative entropies.
     % For the average distribution of all 20 amino acids.
     
     hmm_row(j)=input_mat(i,j)*log2(input_mat(i,j)/null_vec(j));
     
     % For a uniform distribution of all 20 amino acids
     
     hmm_row_unif_distr(j)=input_mat(i,j)*log2(20*input_mat(i,j));  

end

% Here we sum all the contributions.

hmm_vec(i)=sum(hmm_row(:));
hmm_vec_unif_distr(i)=sum(hmm_row_unif_distr(:));

end

% We subtract the minimum value to get all positive numbers as in a HMM 
% logo.

hmm_vec=hmm_vec-min(hmm_vec);
hmm_vec_unif_distr=hmm_vec_unif_distr-min(hmm_vec_unif_distr);

end


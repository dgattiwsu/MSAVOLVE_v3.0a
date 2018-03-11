% This script reorder the covarions and returns their probability 
% distribution. No need to run as a function.

covar_vec_2 = sortrows(covar_vec_rand,1);
covar_vec_2 = sort(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec_2 = sortrows(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec = covar_vec_2;

[ncov,ncovcols] = size(covar_vec);

% Covarying PD of each pair or multiple in the entire MSA.
REF_cov_prob_PD = cell(ncov,1);
REF_cov_alphabet = zeros(nseq,ncovcols+1,ncov);

for i = 1:ncov
nmsa = REF_nmsa(:,nonzeros(covar_vec(i,:)));
nmsa25 = nmsa == 25;
nmsa25 = sum(nmsa25,2);
nogaps = nmsa25 == 0;
nmsa = nmsa(nogaps,:);
alphabet = JointProbDistr_9(nmsa);
    [n_alphabet,alphabet_cols] = size(alphabet);
    REF_cov_alphabet(1:n_alphabet,1:alphabet_cols,i) = alphabet;
    REF_cov_prob_PD{i} = ...
        fitdist((1:1:n_alphabet)','kernel',...
        'width',0.01,'frequency',...
        REF_cov_alphabet(1:n_alphabet,alphabet_cols,i));    
    clear nmsa nmsa25 nogaps alphabet n_alphabet
end

% Covarying PD of each pair or multiple in the main branches.
REF_branch_cov_prob_PD = cell(ncov,nbranches);
REF_branch_cov_alphabet = zeros(nseq,ncovcols+1,ncov,nbranches);

for k = 1:nbranches
REF_branch_ind = find(REF_Clusters == k);
REF_nmsa_branch = REF_nmsa(REF_branch_ind,:);    

for i = 1:ncov
nmsa = REF_nmsa_branch(:,nonzeros(covar_vec(i,:)));
nmsa25 = nmsa == 25;
nmsa25 = sum(nmsa25,2);
nogaps = nmsa25 == 0;
nmsa = nmsa(nogaps,:);
alphabet = JointProbDistr_9(nmsa);
    [n_alphabet,alphabet_cols] = size(alphabet);
    REF_branch_cov_alphabet(1:n_alphabet,1:alphabet_cols,i,k) = alphabet;
    REF_branch_cov_prob_PD{i,k} = ...
        fitdist((1:1:n_alphabet)','kernel',...
        'width',0.01,'frequency',...
        REF_branch_cov_alphabet(1:n_alphabet,alphabet_cols,i,k));    
    clear nmsa nmsa25 nogaps alphabet n_alphabet
end
clear REF_nmsa_branch
end

clear covar_vec_2 covar_vec_rand

% The following line is critical: don't remove it.
covar_vec_bk = covar_vec;
covar_vec_large = covar_vec;
cov_coord = covar_vec;
ncov_coord = size(cov_coord,1);


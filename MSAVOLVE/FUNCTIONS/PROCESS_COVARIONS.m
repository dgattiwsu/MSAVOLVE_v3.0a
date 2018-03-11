% This script reorder the covarions and returns their probability 
% distribution. No need to run as a function.

covar_vec_2 = sortrows(covar_vec_rand,1);
covar_vec_2 = sort(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec_2 = sortrows(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec = covar_vec_2;

cov_prob_PD = cell(ncov,1);
cov_alphabet = zeros(400,3,ncov);

for i = 1:ncov
    A = covar_vec(i,1);
    B = covar_vec(i,2);
    gap=find((REF_nmsa(:,A) ~= 25) & (REF_nmsa(:,B) ~= 25));    
    alphabet = JointProbDistr_8([REF_nmsa(gap,A) REF_nmsa(gap,B)]); 
    n_alphabet = size(alphabet,1);
    cov_alphabet(1:n_alphabet,1:3,i) = alphabet;
    cov_prob_PD{i} = ...
        fitdist((1:1:n_alphabet)','kernel',...
        'width',0.01,'frequency',...
        cov_alphabet(1:n_alphabet,3,i));    
end

REF_cov_prob_PD = cov_prob_PD;
REF_cov_alphabet = cov_alphabet;

clear gap alphabet n_alphabet

% Here we determine the covarying PD of each pair in the main branches.

REF_branch_cov_prob_PD = cell(ncov,nbranches);
REF_branch_cov_alphabet = zeros(400,3,ncov,nbranches);

for k = 1:nbranches
REF_branch_ind = REF_Clusters == k;
REF_nmsa_branch = REF_nmsa(REF_branch_ind,:);    
    for i = 1:ncov
    A = covar_vec(i,1);
    B = covar_vec(i,2);
    gap=find((REF_nmsa_branch(:,A) ~= 25) ...
        & (REF_nmsa_branch(:,B) ~= 25));    
    alphabet = JointProbDistr_8([REF_nmsa_branch(gap,A) ...
        REF_nmsa_branch(gap,B)]); 
    n_alphabet = size(alphabet,1);
    REF_branch_cov_alphabet(1:n_alphabet,1:3,i,k) = alphabet;
    REF_branch_cov_prob_PD{i,k} = ...
        fitdist((1:1:n_alphabet)','kernel',...
        'width',0.01,'frequency',...
        REF_branch_cov_alphabet(1:n_alphabet,3,i,k));    
    clear gap alphabet n_alphabet
    end
clear REF_branch_ind REF_nmsa_branch
end

% clear covar_vec_2


% % Here we determine the covarying PD of each pair in the main branches.
% 
% REF_branch_cov_prob_PD = cell(ncov,nbranches);
% REF_branch_cov_alphabet = zeros(400,3,ncov,nbranches);
% 
% for k = 1:nbranches
% REF_branch_ind = find(REF_LeafClusters == k);
% REF_o_nmsa_branch = REF_o_nmsa(REF_branch_ind,:);    
%     for i = 1:ncov
%     A = covar_vec(i,1);
%     B = covar_vec(i,2);
%     % Here we use only the branch of the reordered msa.
%     gap=find((REF_o_nmsa_branch(:,A) ~= 25) ...
%         & (REF_o_nmsa_branch(:,B) ~= 25));    
%     alphabet = JointProbDistr_8([REF_o_nmsa_branch(gap,A) ...
%         REF_o_nmsa_branch(gap,B)]); 
%     n_alphabet = size(alphabet,1);
%     REF_branch_cov_alphabet(1:n_alphabet,1:3,i,k) = alphabet;
%     REF_branch_cov_prob_PD{i,k} = ...
%         fitdist((1:1:n_alphabet)','kernel',...
%         'width',0.01,'frequency',...
%         REF_branch_cov_alphabet(1:n_alphabet,3,i,k));    
%     end
% clear REF_o_nmsa_branch
% end
% 
% clear gap alphabet n_alphabet

clear covar_vec_2 covar_vec_rand

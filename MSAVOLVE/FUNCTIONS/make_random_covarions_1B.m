function [ nmsa ] = ...
    make_random_covarions_1B( nmsa,branch_start,branch_end,...
    covar_vec,cov_prob_PD,cov_alphabet,ncycles )

% Here we evolve the covarying position of a protein inside a single matrix.

[nseq,npos] = size(nmsa);

% First we create an arrays for the covarying positions.

msa_cov_pos_mat = false(nseq,npos);

% Now we loop through each covarying pair and we select the indices for the 
% covarying positions. We find the occurrence of a member of the covarion list  
% in the logical vector of positions flagged for covariation. Since the 
% covar_vec list contains multiple entries for all possible multiplets we 
% select random from these multiple entries. This is achieved in a somewhat
% complicated way by randomizing the covar_vec list and picking only the
% last occurrence of a unique covarion in the randomized list.
perm_ind = randperm(size(covar_vec,1))';
perm_covar_vec = covar_vec(perm_ind,:);
[unique_perm_covar_vec,last_perm_ind,~] = unique(perm_covar_vec(:,1));
last_perm_covar_vec =[unique_perm_covar_vec last_perm_ind];

for k = 1:size(last_perm_ind)
    i = perm_ind(last_perm_ind(k));
    ind = nonzeros(covar_vec(i,:));
    for j = branch_start:branch_end
        msa_cov_pos_mat(j,ind) = true;
    end
end
    
% Now we evolve only the covarying positions.
% Here we are going to scan along the sequences rather than along the
% positions. If there are multiplets, the last pair in a multiplets gets
% assigned based on its joint PD. This produces a bias in the jointPD's
% based on the order of the pairs in the covar_vec. For this reason, each
% time we scramble the indices of the covar_vec. This is indeed equivalent
% to assuming that if A and B are coevolving and B and C are coevolving,
% sometimes A and B go first and sometimes B and C, before all three become
% segregated.

[ncols,nrows] = size(nmsa);
nmsa1 = nmsa;
nmsa_s = zeros(ncols,nrows);
nmsa1_s = zeros(ncols,nrows);

perm_vec = randperm(size(covar_vec,1));

% Here we take the indices from the permuted vector. 
for k = 1:size(covar_vec,1)
    i = perm_vec(k);
    % Here we loop over the number of cycles requested to maximize the
    % frequency of the covarions
    for cycle = 1:ncycles
    % Here we loop over all the sequences in the msa
    for j = branch_start:branch_end
        % Here for every sequence we check if the covarying positions
        % represented in the loop are flagged.
        covars = nonzeros(covar_vec(i,:));
        ncovars = numel(covars);
        if msa_cov_pos_mat(j,covars)
            % Here we get the covarying residues.
            ind = round(random(cov_prob_PD{i}));
            nmsa1(j,covars) = cov_alphabet(ind,1:ncovars,i);
            nmsa1_s(j,covars) = cov_alphabet(ind,end,i);
        end
    end
    evolved = nmsa1_s >= nmsa_s;
    nmsa(evolved) = nmsa1(evolved);
    nmsa_s = nmsa1_s;
    end
end

end




function [ evolved_MSA,history,mut_history,cov_history,...
    COV,mut_COV,cov_COV,glob_COV,mutcov_COV,...
    msa_mut_pos_mat,msa_cov_pos_mat] = ...
    evolve_6s_3D_cov_hist_multiplets( MSA,history,mut_history,cov_history,...
    COV,mut_COV,cov_COV,glob_COV,mutcov_COV,mut_rate,covar_vec,...
    mut_prob_PD,cov_prob_PD,cov_alphabet,~ )

% Here we evolve a protein inside a single matrix containing n identical
% rows corresponding to the ancestral protein.
% For example we will generate a random number from a Poisson distribution 
% assuming that 10% of all the sequence positions may have changed during 
% evolution in an undefined amount of time, which is covered by the 
% evolution cycle. Thus our target number 'nmut' (lambda of 
% the Posson distribution) is 300*.1 = 30. We will draw a number 
% once for each sequence in the msa.

% symbols = (1:1:20)';

% Please, REMEMBER! the 'nseq' inside this function is not the same value
% as 'nseq' in the main program. Here 'nseq' has a different value in each 
% of the three levels.
[nseq,npos] = size(MSA);

% MSAs_mat = zeros(nseq,npos);
% for i = 1:nseq
%     for j = 1:npos  
%     MSAn = MSA(i,j);        
%     MSAs_mat(i,j) = mut_freq_PD(MSAn,j);
%     end
% end

msa_mut_num_vec = zeros(nseq,1);
nmut = round(npos*mut_rate/100);
for i = 1:nseq
msa_mut_num_vec(i,1) = random('Poisson',nmut);
end

% Now for each sequence we will draw the positions that must change from a
% uniform distribution of npos positions. 
% First we create two arrays, one for the normal positions and one for the
% covarying positions.

msa_mut_pos_mat = false(nseq,npos);
msa_cov_pos_mat = false(nseq,npos);
for i = 1:nseq
    %   For a uniform distribution
    ind = randi(npos,1,msa_mut_num_vec(i));
    %   For a distribution previously specified (like conserv_PD).
    %   ind = round(random(conserv_PD,1,msa_mut_num_vec(i)));
    msa_mut_pos_mat(i,ind) = true;
end

% Now we loop through each covarying pair and we split the indices for all
% the normal positions from the indices for the covarying positions. We
% need to do this in tree passes.
% 1st pass: we find the occurrence of a member of the covarion list in the 
% logical vector of positions flagged for mutation. Since the covar_vec
% list contains multiple entries for all possible multiplets we select
% random from these multiple entries. This is achieved in a somewhat
% complicated way by randomizing the covar_vec list and picking only the
% last occurrence of a unique covarion in the randomized list.
perm_ind = randperm(size(covar_vec,1))';
perm_covar_vec = covar_vec(perm_ind,:);
[unique_perm_covar_vec,last_perm_ind,~] = unique(perm_covar_vec(:,1));
last_perm_covar_vec =[unique_perm_covar_vec last_perm_ind];

for k = 1:size(last_perm_ind)
    i = perm_ind(last_perm_ind(k));
    ind = nonzeros(covar_vec(i,:));
    for j = 1:nseq
        if sum((msa_mut_pos_mat(j,ind))) > 0
        msa_cov_pos_mat(j,ind) = true;
        end
    end
 end

% 2nd pass: we remove from the mutation list any position that was put in
% the covarion list.
 for j = 1:nseq
        msa_mut_pos_mat(j,msa_cov_pos_mat(j,:)) = false;
 end

% First we evolve the standard positions and we accept all mutations.
MSA1 = MSA;
for j = 1:npos
    ind = find(msa_mut_pos_mat(:,j));
    for n = 1:length(ind)
        MSA1(ind(n),j) = round(random(mut_prob_PD{j}));
    end
end

% We can make the evolution more "selective" by accepting only changes that
% correspond to higher values in the probability distribution. This will
% lead to significantly more positions becoming fully conserved. Use 3C for
% this purpose.
% MSA1s_mat = zeros(nseq,npos);
% for i = 1:nseq
%    for j = 1:npos  
%    MSA1n = MSA1(i,j);        
%    MSA1s_mat(i,j) = mut_freq_PD(MSA1n,j);
%    end
% end
%
% evolved = MSA1s_mat >= MSAs_mat;
% MSA(evolved) = MSA1(evolved);
% MSA = MSA1;

history_part_1 = MSA1 ~= MSA;
% IMPORTANT: convert to numeric or unpredictable results may occur.
history_part_1 = + history_part_1;

    for k = 1:nseq
        COV_part = zeros(npos,npos,nseq);
        for m = 1:npos
            for j = m:npos    
                COV_part(m,j,k) = history_part_1(k,m)*history_part_1(k,j);
            end
        end
        COV = COV + COV_part;    
        mut_COV = mut_COV + COV_part;    
    end
    
% Now we evolve only the covarying positions and we accept all mutations.
% Here we are going to scan along the sequences rather than along the
% positions. If there are multiplets, the last pair in a multiplets gets
% assigned based on its joint PD. This produces a bias in the jointPD's
% based on the order of the pairs in the covar_vec. For this reason, each
% time we scramble the indices of the covar_vec. This is indeed equivalent
% to assuming that if A and B are coevolving and B and C are coevolving,
% sometimes A and B go first and sometimes B and C, before all three become
% segregated.

MSA2 = MSA1;
perm_vec = randperm(size(covar_vec,1));

% While only the last pair in a multiplet is implemented it is
% actually different to go through all the pairs rather than selecting only 
% the last one.
% perm_covar_vec = covar_vec(perm_vec,:);
% [b,m,n] = unique(perm_covar_vec(:,1))

% Here we loop over the covarying positions, but we take the indices from
% the permuted vector.

for k = 1:size(covar_vec,1)
    i = perm_vec(k);
    % Here we loop over all the sequences in the msa
    for j = 1:nseq
        % Here for every sequence we check if the covarying positions
        % represented in the loop are flagged
        if msa_cov_pos_mat(j,nonzeros(covar_vec(i,:)))
            % Here we get the covarying residues.
            ind = round(random(cov_prob_PD{i}));
            MSA2(j,nonzeros(covar_vec(i,1))) = cov_alphabet(ind,1,i);
        end
    end
end
%--------------------------------------------------------------------------    

history_part_2 = MSA2 ~= MSA1;
% IMPORTANT: convert to numeric or unpredictable results may occur.
history_part_2 = + history_part_2;

    for k = 1:nseq
        COV_part = zeros(npos,npos,nseq);
        for m = 1:npos
            for j = m:npos    
                COV_part(m,j,k) = history_part_2(k,m)*history_part_2(k,j);
            end
        end
        COV = COV + COV_part;    
        cov_COV = cov_COV + COV_part;    
    end

history_part_3 = MSA2 ~= MSA;
% IMPORTANT: convert to numeric or unpredictable results may occur.
history_part_3 = + history_part_3;

    for k = 1:nseq
        COV_part = zeros(npos,npos,nseq);
        for m = 1:npos
            for j = m:npos    
                COV_part(m,j,k) = history_part_3(k,m)*history_part_3(k,j);
            end
        end
        glob_COV = glob_COV + COV_part;    
        mutcov_COV = mutcov_COV + COV_part;    
    end

    
evolved_MSA = MSA2;
% IMPORTANT:    make sure the history arrays are all numeric or 
%               unpredictable results may occur.
history = + history;
mut_history = + mut_history;
cov_history = + cov_history;
history = history + history_part_1 + history_part_2;
mut_history = mut_history + history_part_1;
cov_history = cov_history + history_part_2;

end




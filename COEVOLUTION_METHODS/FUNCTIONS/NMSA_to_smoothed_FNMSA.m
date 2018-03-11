function [fnmsa,A_hat] = NMSA_to_smoothed_FNMSA(nmsa,dist_method,threshold,...
    psc_method,psc_lambda,psc_blosum,nsymbols,recover_method,recover_level)
% This function converts every symbol in the msa to the frequency of that
% symbol in each column. Weights and pseudocounts are applied to the
% frequencies to correct for the similarities between sequences.
% dist_method = 'GAPS'|'NOGAPS' (include or not gaps in the distance matrix).
% threshold = (similarity threshold: e.g. 0.9).
% psc_method = 'DCA'|'SDP' (pseudocount method DCA or SDP style as described 
% in "SDPpred ..." Kalinina, OV et al. Nucleic Acid Research 2004, Vol32, 
% W424-W428, and in "H2r ..." Merkl, R. and Zwick, M., BMC Bioinformatics 
% 2008, 9,151).
% psc_lambda = pseudocount scale (e.g. 0.5 for DCA; 1.0 for SDP).
% psc_blosum = 'AUTO'|30|35|40|45|50|55|60|62|65|70|75|80|85|90|100. Blosum  
% matrix used to determine the pseudocounts in the SDP method. If set to 
% 'AUTO' the best blosum matrix for the mean value of the distance matrix
% is chosen automatically.
% nsymbols = 20|21 (if set to 20 gaps frequencies are zeroed; 20 is 
% recommended if using psc_method = 'DCA'.
% recover_method = 'ALM'|'iALM'|'APG'|'pAPG'|'NONE'; it is the robust PCA
% method used to recover the low rank matrix corresponding to the frequency
% converted msa with noise removed. It uses one of the methods developed by
% the Perception and Decision Lab at the University of Illinois
% (http://perception.csl.uiuc.edu/). The options are: 'exact ALM' (ALM)[can
% be very slow], 'inexact ALM' (iALM)[strongly recommended], 'accelerated
% proximal gradient' (APG) , 'partial accelerated proximal gradient' (pAPG)
% [recommended], or no recovering (NONE). 
% recover_level = specifies the degree of noise removed: useful values are 
% in the range 0.2-0.05. If using recovering (e.g. recover_method = iALM),
% the recommended value is recover_level = 0.05.

% Possible usage:
% [fnmsa] = NMSA_to_smoothed_FNMSA(nmsa,'NOGAPS',0.90,'SDP',0.5,'AUTO',21,...
%                                  'NONE',0.05);
% [fnmsa,noise_red_fnmsa] = NMSA_to_smoothed_FNMSA(nmsa,'NOGAPS',0.90,...
%                                  'SDP',0.5,'AUTO',21,'iALM',0.05);

[nrows,ncols] = size(nmsa);

% Replace unusual symbols
ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

% Calculate the distance matrix.

[dist] = get_distance_matrix(nmsa,dist_method);

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff=round(sum(W));
end

fprintf('Meff = %d \n', Meff);

% Here we create a 3d matrix in which every layer has the dimensions of the
% nmsa and represents one of 20 symbols (if gaps are excluded). Then every row
% of each layer is scaled by the weight of that sequence.

nmsa_3 = zeros(nrows,ncols,21);

for i = 1:nsymbols
    nmsa_3(:,:,i) = nmsa == i;
    for j = 1:nrows
        nmsa_3(j,:,i) = nmsa_3(j,:,i)*W(j);
    end
end

% We have two options to determine pseudocount corrected frequencies:

switch psc_method
    case 'SDP'
    % SDPpred style (much more complicated).
        
    Fi = zeros(21,ncols);
    Si = zeros(21,ncols);

    % First we determine the profile as total counts.

    for i = 1:21
        layer = squeeze(nmsa_3(:,:,i));
        Si(i,:) = sum(layer);
    end

% Here we retrieve the blosum substitution matrix and we convert into
% probabilities. In order to complete this operation we will also need to 
% obtain the background probabilities  either 'in general' or in the 
% specific msa that is being studied. We recall that the entries Sij in the  
% blosum matrix are Sij=lambda*log2(Mij/pj), where Mij is the probability of 
% i mutating to j and pj is the frequency of j.

if psc_blosum == 'AUTO'
    mean_dist = mean(sum(dist)-ones(1,nrows))*100/nrows;
    blosum_vec = [30 35 40 45 50 55 60 62 65 70 75 80 85 90 100];
    blosum_vec_dif = abs(blosum_vec - mean_dist);
    [~,blosum_ind] = min(blosum_vec_dif);
    psc_blosum = blosum_vec(blosum_ind);
    fprintf('Using the Blosum%d matrix \n', psc_blosum);
end
    
    [blosum_mat,info] = blosum(psc_blosum);
    blosum_mat = blosum_mat(1:20,1:20);
    scale = info.Scale;
    blosum_mat = pow2(blosum_mat*scale);
% 
    Si_sum = sum(Si(1:20,:));
    fSi = zeros(20,ncols);

% Here we calculate the background probabilities for this msa.
    for i = 1:20
        fSi(i,:) = Si(i,:)./Si_sum;
    end
    bg_prob = mean(fSi,2);

    Si_stack = zeros(21,ncols,21);
    % Si_sum_stack = zeros(20,ncols);

    for i = 1:20
        for l = 1:20
            % Si_stack(l,:,i) = Si(l,:)*blosum_mat(l,i);
            % Si_stack(l,:,i) = fSi(l,:).*Si(l,:)*blosum_mat(l,i);
            % This is the probability that each aa 'l' will mutate to aa
            % 'i' at every position.
            Si_stack(l,:,i) = bg_prob(l).*Si(l,:)*blosum_mat(l,i);
        end
        Si_stack(i,:,i) = 0;
        % Si_sum_stack(i,:) = Si_sum - Si(i,:);
    end

    psc_Si_sum = Si_sum + psc_lambda.*sqrt(Si_sum);

    for i = 1:20
        % psc_sum = (sum(Si_stack(:,:,i)))./Si_sum_stack(i,:);
        % psc_sum = (sum(Si_stack(:,:,i)))./Si_sum;
        % This is the total probability that any aa would convert to aa 'i'
        % at each position.
        psc_sum = (sum(Si_stack(:,:,i)))./sqrt(Si_sum);
        Fi(i,:) = Si(i,:) + psc_lambda*psc_sum;
        Fi(i,:) = Fi(i,:)./psc_Si_sum;
    end

% Here we make the assumption that the frequency of the gaps is whatever is
% left of the unobserved frequency.

    Fi_sum = sum(Fi);
    Fi(21,:) = 1 - Fi_sum;

    case 'DCA'
    % DCA style
    loq = psc_lambda/21;
    Fi = zeros(21,ncols);
    for i = 1:21
        layer = squeeze(nmsa_3(:,:,i));
        % Pseudocount.
        Fi(i,:) = (loq + sum(layer))/(loq + Meff);
    end

end

% Here we remove the frequencies of the gaps and we scale the remaining
% frequencies.

if nsymbols == 20
    Fi(21,:) = 0;
end
    Fi_sum = sum(Fi);
    scaleFi = repmat(Fi_sum,21,1);
    Fi = Fi./scaleFi;

% Here we create a new msa in which every symbol is replaced by the
% weigthed frequence of that symbol in that column of the msa.

fnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        fnmsa(i,j) = Fi(row(j),j);
        end
end

% MI = NMSA_to_MI(fnmsa);

% First we recover the msa matrix with noise removed using a technique of
% robust PCA.
switch recover_method
    case 'ALM'
    [A_hat E_hat iter] = exact_alm_rpca(fnmsa,recover_level);
    case 'iALM'        
    [A_hat E_hat iter] = inexact_alm_rpca(fnmsa,recover_level);
    case 'APG'
    [A_hat E_hat iter] = proximal_gradient_rpca(fnmsa,recover_level);
    case 'pAPG'
    [A_hat E_hat iter] = partial_proximal_gradient_rpca(fnmsa,recover_level);
    case 'NONE'
    A_hat = fnmsa;
end

end

function [dist] = get_distance_matrix(nmsa,dist_method)

[nrows,ncols] = size(nmsa);

switch dist_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(nmsa);        
        dist = bin_ordered*bin_ordered';

        sdist = zeros(nrows,nrows);      
        for i = 1:nrows
            sdist(i,:) = dist(i,:)/dist(i,i);
        end

        udist = triu(sdist,1);
        ldist = tril(sdist,-1)';
        lind = udist < ldist;
        udist(lind) = 0;
        ldist(~lind) = 0;
        mdist = udist + ldist;
        mdist = mdist + mdist' + eye(nrows);

        dist = mdist;
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;

end

end

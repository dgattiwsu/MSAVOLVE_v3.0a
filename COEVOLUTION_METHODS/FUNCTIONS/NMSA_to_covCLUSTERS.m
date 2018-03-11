function [ fnmsa,Fi,A_hat,E_hat,...
    o_COV,u_COV,f_COV,...
    o_data1,o_center1,o_dist1,...
    u_data1,u_center1,u_dist1,...
    f_data1,f_center1,f_dist1,...
    o_data2,o_center2,o_dist2,...
    u_data2,u_center2,u_dist2,...
    f_data2,f_center2,f_dist2,...    
    o_data3,o_center3,o_dist3,o_icaA,...
    u_data3,u_center3,u_dist3,u_icaA,...
    f_data3,f_center3,f_dist3,f_icaA ] = ...
    NMSA_to_covCLUSTERS(nmsa,threshold,psc_lambda,nsymbols,...
    recover_method,noise_corr,ica_distance_method,edge_scale,plotoption,...
    stabil_1,stabil_2,stabil_3)

% This function calculates a covariance matrix without first converting the
% MSA to a long binary format. This is achieved by converting every symbol
% in the msa to the frequency of that symbol in each column. Weights and
% pseudocount are applied to the frequencies to correct for the
% similarities between sequences. While this type of covariance matrix is
% not very good for the purpose of finding 'strictly' covarying pairs, its
% spectral analysis leads to the identification of clusters of position
% that have similar patterns of covariation. 'recover_method' is the robust
% PCA method used to recover the low rank matrix corresponding to the
% frequency converted msa with noise removed. 'noise_corr' is a value that
% define the expected amount of noise correction: a value of 1 produces no
% changes in the msa; smaller values produce a progressive reduction of the
% noise [recommended: 0.05] of the msa. It uses one of the methods
% developed by the Perception and Decision Lab at the University of
% Illinois (http://perception.csl.uiuc.edu/). The options are: 'exact ALM'
% (ALM)[can be very slow], 'inexact ALM' (iALM)[recommended], 'accelerated
% proximal gradient' (APG)[recommended] ,'partial accelerated proximal
% gradient' (pAPG)[recommended], or no recovering (NONE). The most stable
% method is APG, but iALM may give better clustering. A noise correction
% value of 0.05 gives almost always good results. Smaller values tend to
% progressively remove differences between the columns of the msa, and thus
% ultimately group 'all' residues in 2-3 clusters. 'ica_distance_method'
% defines whether the ICA distance between points are calculated on the
% orthogonal 'ORTHO' or skewed 'SKEW' representation. Independent Component
% Analysis is carried out with the 'fastICA' package by Aapo Hyvarinen at
% the University of Helsinki.

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

% First we calculate a traditional positional covariance matrix on the
% 'nmsa' . For this purpose we find convenient to use a reduced output from
% the function NMSA_to_posCOV:
 
[o_COV] = ...
     NMSA_to_posCOV(nmsa,threshold,psc_lambda,'FRO',0.0,0.001,nsymbols);
 
[o_data1,o_dist1,o_center1,o_data2,o_dist2,o_center2,...
    o_data3,o_dist3,o_center3,o_icaA] = ...
    find_3d_dist(o_COV,ica_distance_method,edge_scale,plotoption,stabil_1);

% Now we will calculate the positional covariance matrix in a completely
% new way. First we calculate the sequence distance matrix.

bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';

dist_threshold = dist > threshold*ncols; 

if threshold == 1
    W = ones(nrows,1);
else
W = 1./sum(dist_threshold,2);
end
Meff=sum(W);

fprintf('Meff = %f \n', Meff);

% lambda for pseudocount.
loq = psc_lambda/21;

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

% Here we calculate the profile (Fi) of the nmsa taking into consideration 
% the pseudocount:

Fi = zeros(21,ncols);
for i = 1:21
    layer = squeeze(nmsa_3(:,:,i));
    % Pseudocount.
    Fi(i,:) = (loq + sum(layer))/(loq + Meff);
end

% If we removed the frequencies of the gaps here we scale the remaining
% frequencies.

if nsymbols == 20
sumFi = sum(Fi);
scaleFi = repmat(sumFi,21,1);
Fi = Fi./scaleFi;
end

% Here we create a new msa in which every symbol is replaced by the
% weigthed frequence of that symbol in that column of the msa.
fnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        fnmsa(i,j) = Fi(row(j),j);
        end
end

% Here we can show that the MI matrix is unchanged: we have not lost the
% information on covariation by converting symbols to frequencies.
% MI_1 = NMSA_to_MI(nmsa);
% MI_2 = NMSA_to_MI(fnmsa);

% Here we calculate the positional covariance matrix for this
% frequency-based fnmsa

u_COV = cov(fnmsa);
[u_data1,u_dist1,u_center1,u_data2,u_dist2,u_center2,...
    u_data3,u_dist3,u_center3,u_icaA] = ...
    find_3d_dist(u_COV,ica_distance_method,edge_scale,plotoption,stabil_2);

% Finally we calculate an fnmsa matrix with noise removed using a technique of
% robust PCA.

switch recover_method
    case 'ALM'
    [A_hat E_hat iter] = exact_alm_rpca(fnmsa,noise_corr);
    case 'iALM'        
    [A_hat E_hat iter] = inexact_alm_rpca(fnmsa,noise_corr);
    case 'APG'
    [A_hat E_hat iter] = proximal_gradient_rpca(fnmsa,noise_corr);
    case 'pAPG'
    [A_hat E_hat iter] = partial_proximal_gradient_rpca(fnmsa,noise_corr);
    case 'NONE'
    A_hat = fnmsa;
end

f_COV = cov(A_hat);
[f_data1,f_dist1,f_center1,f_data2,f_dist2,f_center2,...
    f_data3,f_dist3,f_center3,f_icaA] = ...
    find_3d_dist(f_COV,ica_distance_method,edge_scale,plotoption,stabil_3); 
  
end


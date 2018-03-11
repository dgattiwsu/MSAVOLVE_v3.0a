function [o_cov_mat] = sort_cov_mat_2(cov_mat)
% First we calculate the distance 'd' between all the sequences using the
% p-distance method, where p is the proportion of sites at which the two 
% sequences differ. p is close to 1 for poorly related sequences, and is 
% close to 0 for similar sequences.
[nseq,~] = size(cov_mat);

% Here we calculate the tree and we cluster the sequences. 

 tree = seqlinkage(1-cov_mat);

% We set a limit to the number of clusters detected to avoid errors.

 [LeafClusters,~,~] = ...
     cluster(tree,[],'criterion','gain','MaxClust',6);
  
% Here we reorder the msa and the distance matrix based on the phylogenetic 
% tree.

 leaf_name_order = get(tree,'LeafNames');
 leaf_number_order = zeros(nseq,1);

 for i = 1:nseq
    name = leaf_name_order{i,1}; 
    leaf_number_order(i,1) = str2double(name(6:end));
 end

% However, clustering can give results different from the tree ranking. For
% this reason we combine the two informations to achieve a more careful
% reordering of the sequences that reflects both the phylogenetic tree and
% the clustering. This is useful to display the heat map of the covariance
% matrix (see below).

 cluster_ind = [LeafClusters leaf_number_order];
 cluster_ind = sortrows(cluster_ind,1);
 cluster_ind = cluster_ind(:,2);

 [V,D] = spectral(cov_mat);
 V = V(:,cluster_ind);
 D = diag(D);
 D = D(cluster_ind);
 D_mat = zeros(nseq);

    for i = 1:nseq
    D_mat(i,i) = D(i,1);
    end
     
 o_cov_mat = V*D_mat*V';

 

     
  

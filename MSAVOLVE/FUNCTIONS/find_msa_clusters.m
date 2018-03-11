function [ cluster,cluster_ind,cluster_COV ] = find_msa_clusters(nmsa,...
    COV_3D,DISTANCE,DIMENSIONS,MAX_CLUSTERS,cluster_number)

% This functions identifies the indices of a requested cluster in a MSA and
% returns the cluster indices, the cluster MSA and the cluster totCOV
% matrix.
npos = size(nmsa,2);
cmsa = int2aa(nmsa);
 if DISTANCE == 1
 dist = seqpdist(cmsa,'SquareForm',true,'Method','p-distance');
 cov_mat = 1-dist;
 else 
 binmsa = nmsa_to_binmsa(nmsa);
 cov_mat = cov(binmsa',1);  
 end
[pc,~] = spectral(cov_mat);
clusters = clusterdata(pc(:,1:DIMENSIONS),'maxclust',MAX_CLUSTERS);

cluster_ind = find(clusters == cluster_number);
cluster = nmsa(cluster_ind,:);

COV = sum(COV_3D(:,:,cluster_ind),3);
COV = COV+COV';
for i = 1:npos
COV(i,i) = NaN;
end
cluster_COV = COV;
% cluster_COV = COV_3D;
end


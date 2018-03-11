function [ cluster,cluster_ind,cluster_COV ] = find_msa_branches(nmsa,...
    COV_3D,b_ind2,branch_number)

% This functions identifies the indices of a requested branch in a MSA and
% returns the branch indices, the branch MSA and the branch totCOV matrix.
npos = size(nmsa,2);
cmsa = int2aa(nmsa);

cluster_ind = (b_ind2(branch_number,1):b_ind2(branch_number,2));
cluster = nmsa(cluster_ind,:);

COV = sum(COV_3D(:,:,cluster_ind),3);
COV = COV+COV';
for i = 1:npos
COV(i,i) = NaN;
end
cluster_COV = COV;
% cluster_COV = COV_3D;
end


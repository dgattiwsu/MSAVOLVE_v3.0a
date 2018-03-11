function[all_sca_mat,nv_sca_mat,nm_sca_mat]=get_sca(DDG_mat)
% Here the imput is the DDGation matrix containing the trajectory of
% all the possible aa at the different positions in all the trials. 
% From it we derive the larger matrix containing all the possible dot
% products of these vectors (all_sca_mat) and the npos*npos matrix of sca
% values. Either the vector (Frobenius) norm is used (nv_sca_mat) or the 
% matrix norm (largest singular value, nm_sca_mat).
[naa,npos,ntrial] = size(DDG_mat);
all_sca_mat = zeros(npos,npos,naa,naa);
nv_sca_mat = zeros(npos,npos);
nm_sca_mat = zeros(npos,npos);

for i = 1:npos
    for j = i:npos
        sca_mat_i = squeeze(DDG_mat(:,i,:));
        sca_mat_j = squeeze(DDG_mat(:,j,:));
        sca_mat_ij = sca_mat_i * sca_mat_j';
        all_sca_mat(i,j,:,:) = sca_mat_ij;
        nv_sca_mat(i,j) = (norm(sca_mat_ij(1:end)))/ntrial;
        nm_sca_mat(i,j) = (norm(sca_mat_ij))/ntrial;        
    end
end

nv_sca_mat = nv_sca_mat + nv_sca_mat';
nm_sca_mat = nm_sca_mat + nm_sca_mat';

for i = 1:npos
    nv_sca_mat(i,i) = NaN;
    nm_sca_mat(i,i) = NaN;
end

end

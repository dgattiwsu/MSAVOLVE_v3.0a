% This script exports the fasta alignment and the cov_COV matrix from a 
% MSAvolve run both as such and as a covar_vec.
for i = 1:end_cycle
    s1 = ['kdo8ps_msa_' int2str(i) '.fasta'];
    nmsa_to_faln(MSA_select_ALL(:,:,i),s1);
    s2 = ['kdo8ps_msa_' int2str(i) '_coev.matrix'];
    dlmwrite(s2,cov_COV_ALL(:,:,i));    
    sorted_COV = sort_matrix_descend(cov_COV_ALL(:,:,1));
    ind = find(sorted_COV(:,1));
    sorted_COV = sorted_COV(ind,:);
    s3 = ['kdo8ps_msa_' int2str(i) '_covar_vec.txt'];
    dlmwrite(s3,sorted_COV);
end

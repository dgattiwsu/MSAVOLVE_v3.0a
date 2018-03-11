function [s_coev_mat,s_inter_mat] = cum_dist_coev(coev_mat,inter_mat,near)
% coev_mat: coevolution matrix derived from a method
% inter_mat: sparse interaction matrix derived from 'coev_distance_matrix'
% near: threshold to cut out residues close in sequence. A
% value of 1 means residues are accepted if separated by 1 (e.g. 31 and
% 32). A value of 2 means residues are accepted if separated by 2 (e.g. 31 and
% 33). And so on ...

% s_coev_mat = sort_matrix_descend_2(coev_mat,near);
s_coev_mat = sort_matrix_descend_2(coev_mat,1);

% cutoff = s_coev_mat(npreds,1);
% thr = coev_mat>=cutoff;
% inter_mat = full(inter_mat);
% inter_mat(~thr) = 0;
% inter_u = triu(inter_mat,1);
% inter_l = tril(inter_mat,-1);
% inter_mat = inter_u + inter_l;

s_inter_mat = sort_matrix_descend_2(inter_mat,near);
nnz_ind = nnz(s_inter_mat(:,1));
[c,ia,ib] = intersect(s_coev_mat(:,2:3),s_inter_mat(1:nnz_ind,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    s_coev_mat(ia(i),4) = s_inter_mat(ib(i),1);
end
s_coev_mat(:,5)= cumsum(s_coev_mat(:,4));
s_inter_mat(:,4) = 0;
s_inter_mat(:,5) = cumsum(s_inter_mat(:,1));


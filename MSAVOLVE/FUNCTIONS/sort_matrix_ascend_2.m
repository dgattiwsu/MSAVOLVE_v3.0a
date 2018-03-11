function [ sorted_mat ] = sort_matrix_ascend_2( mat,order )
% This function sorts a symmetric matrix in descending order converting all 
% the NaN to zero before the sorting. Then the matrix subscripts of each 
% matrix value are identified and written as two adjacent columns (first 
% rows then cols) next to the sorted values column. Only the unique part of
% the matrix, based on 'order' is retained, the rest is equal to zero.

mat = nantozero(mat);
mat = triu(mat,order);
[rank,ind]=sort(mat(:),'ascend');
s=size(mat);
[sub_row,sub_col]=ind2sub(s,ind);

% The original linear indices of the matrix are no longer needed.

sorted_mat=[rank sub_row sub_col];

end


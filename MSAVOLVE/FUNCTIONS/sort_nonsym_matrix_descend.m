function [ sorted_mat ] = sort_nonsym_matrix_descend( mat )
% This function sorts a symmetric matrix in descending order converting all 
% the NaN to zero before the sorting. Then the matrix subscripts of each 
% matrix value are identified and written as two adjacent columns (first 
% rows then cols) next to the sorted values column.

% [rank,ind]=sort(nantozero(mat(1:end)),'descend');
% s=size(mat);
% [sub_row,sub_col]=ind2sub(s,ind);
% 
% % The original linear indices of the matrix are no longer needed.
% % sorted_mat=[rank' ind' sub_row' sub_col'];
% 
% sorted_mat=[rank' sub_row' sub_col'];
% 
% % Since the matrix is symmetric save every other result.
% sorted_mat=sorted_mat(2:2:end,:);

% Modified from the original above on 03/26/12 to make compatible with
% sort_matrix_descend_2. By default the code below uses an order of 1 in
% the triu function: this means the diagonal is not included in the
% sorting.

mat = nantozero(mat);
[rank,ind]=sort(mat(:),'descend');
s=size(mat);
[sub_row,sub_col]=ind2sub(s,ind);

% The original linear indices of the matrix are no longer needed.

sorted_mat=[rank sub_row sub_col];

end


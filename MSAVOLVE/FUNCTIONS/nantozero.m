function [ nantozero_matrix ] = nantozero( matrix )
% NANTOZERO: converts all the NaN in a vector or matrix to zero.
nantozero_matrix=matrix;
nan_matrix=isnan(matrix);
nan_matrix_index= nan_matrix==1;
nantozero_matrix(nan_matrix_index)=0;
end


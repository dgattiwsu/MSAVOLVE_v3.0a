function [ mat ] = diag_to_nan( mat )
% This function converts the diagonal of a square matrix to all NaN's.
[~,ncols] = size(mat);
diag_ind = logical(eye(ncols));
mat(diag_ind) = NaN;
end


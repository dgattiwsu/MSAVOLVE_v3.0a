function [ zerotonan_matrix ] = zerotonan( matrix )
% NANTOZERO: converts all zeros in a vector or matrix to NaN.
zerotonan_matrix = matrix;
zero_matrix_index = matrix == 0;
zerotonan_matrix(zero_matrix_index) = NaN;
end


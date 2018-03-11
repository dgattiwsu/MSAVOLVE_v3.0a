function [scale] = scale_matrices(mat1,mat2)
% This function returns the best scale such that mat1*scale = mat2.
% Notice that this is the same as projecting the linearized mat2 vector
% onto the linearized mat1 vector.
[rows,cols] = size(mat1);
template = ones(rows,cols);
upper = triu(template,1);
triu_ind = template == upper;
mat1 = mat1(triu_ind);
mat2 = mat2(triu_ind);
scale = mat1\mat2;
% Possible scales:
% scale1 = mat1\mat2;
% scale1 = (mat1'*mat2)/(mat1'*mat1);
% scale2 = mat2\mat1;
% scale2 = (mat2'*mat1)/(mat2'*mat2);
end

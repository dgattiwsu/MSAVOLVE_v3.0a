function [ZNMI_mat] = MI_to_ZNMI(input_mat)
% This function calculates a ZNMI matrix from an input MI matrix according
% to the algorithm of Brown and Brown.

mat=input_mat;

[rows,cols]=size(mat);
mean_row=zeros(rows,1);
var_row=zeros(rows,1);
mu_nd=zeros(rows,cols);
var_nd=zeros(rows,cols);
sig_nd=zeros(rows,cols);
ZNMI_mat=zeros(rows,cols);

for i=1:rows
    mean_row(i)=nanmean(mat(i,:));
    var_row(i)=nanvar(mat(i,:));   
end

% Only the mean and variance of the column NMI must be calculated according
% to equation 4 in PLoS ONE,5(6),p 13, 2010. When referred to this mean and
% variance each value of the input matrix becomes a value of the ZNMI
% matrix.

for i=1:rows
    for j=i:rows
    mu_nd(i,j)=(mean_row(i)*var_row(j)+mean_row(j)*var_row(i))/...
        (var_row(i)+var_row(j));
    var_nd(i,j)=(var_row(j)*var_row(i))/(var_row(i)+var_row(j));
    sig_nd(i,j)=sqrt(var_nd(i,j));
    ZNMI_mat(i,j)=(mat(i,j)-mu_nd(i,j))/sig_nd(i,j);
    ZNMI_mat(j,i)=ZNMI_mat(i,j);
    end
    
ZNMI_mat(i,i)=NaN;

end

end

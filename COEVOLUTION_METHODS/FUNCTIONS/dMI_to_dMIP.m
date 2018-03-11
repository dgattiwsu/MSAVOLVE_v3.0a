function [ MIP_mat ] = dMI_to_dMIP( input_mat )
% This function calculates the MIp matrix from an input MI matrix according
% to the algorithm of Dunn, Wahl, and Gloor (2008). Differs from MI_to_MIP
% because it keeps the diagonal values.

mat=input_mat;

[rows,cols]=size(mat);
mean_mat=nanmean(mat(1:numel(mat)));
mean_row=zeros(rows,1);
var_row=zeros(rows,1);
MCA_mat=zeros(rows,cols);

% Here  we calculate the MCA matrix.

for i=1:rows
    mean_row(i)=nanmean(mat(i,:));
    var_row(i)=nanvar(mat(i,:));   
end

for i=1:rows
    for j=i:rows

    MCA_mat(i,j)=(mean_row(i)*mean_row(j))/mean_mat;
    MCA_mat(j,i)=MCA_mat(i,j);
    
    end
% MCA_mat(i,i)=NaN;

end

%MCA_mat=MCA_mat+MCA_mat';    

%for n=1:cols
    
%MCA_mat(n,n)=MCA_mat(n,n)/2;

%end

% Finally we subtract the MCA matrix from the MI matrix

MIP_mat=mat-MCA_mat;

end


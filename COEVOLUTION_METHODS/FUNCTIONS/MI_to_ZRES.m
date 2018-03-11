function [ ZRES_mat ] = MI_to_ZRES( input_mat )
% This function calculates the ZRes matrix from a simple MI matrix 
% according to the algorithm by Little and Chen.

mat=input_mat;

[rows,cols]=size(mat);
mean_row=zeros(rows,1);
var_row=zeros(rows,1);
std_row=zeros(rows,1);
CP_mat=zeros(rows,cols);
FIT_mat=zeros(rows,cols);
ZRES_mat=zeros(rows,cols);

% Calculate the column product (CP) matrix 

for i=1:rows
    mean_row(i)=nanmean(mat(i,:));
    var_row(i)=nanvar(mat(i,:));   
end
for i=1:rows
    for j=i:rows

    CP_mat(i,j)=(mean_row(i)*mean_row(j));
    CP_mat(j,i)=CP_mat(i,j);
    end

CP_mat(i,i)=NaN;

end

% Least-square fit the CP matrix to the input matrix. Ignore the warning of
% rank deficiency.

for i=1:rows
    betahat=CP_mat(i,:)'\mat(i,:)';
    FIT_mat(i,:)=CP_mat(i,:)*betahat;
end

% Calculate the residual matrix and make sure it is symmetric.

RES_mat=mat-FIT_mat;
RES_mat=(RES_mat+RES_mat')/2;

% Calculate the ZRES matrix

for i=1:rows
    mean_row(i)=nanmean(RES_mat(i,:));
    std_row(i)=nanstd(RES_mat(i,:));   
end
for i=1:rows
    for j=i:rows
    zscore_1=(RES_mat(i,j)-mean_row(i))/std_row(i);
    zscore_2=(RES_mat(i,j)-mean_row(j))/std_row(j);
    ZRES_mat(i,j)=zscore_1*zscore_2;
    if (zscore_1<0&&zscore_2<0)
        ZRES_mat(i,j)=-ZRES_mat(i,j);
    end
    
    ZRES_mat(j,i)=ZRES_mat(i,j);
        
    end

ZRES_mat(i,i)=NaN;

end

end


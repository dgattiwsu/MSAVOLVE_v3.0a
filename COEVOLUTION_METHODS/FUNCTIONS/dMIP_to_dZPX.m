function [ ZPX_mat,ZPX2_mat ] = dMIP_to_dZPX( input_mat )
% Same as MI_to_ZPX. This function calculates a Zpx and a squared Zpx
% matrix from an input MIP matrix according to the algorith of Gloor et al.
% (2010). The squared Zpx matrix is identical to the ZRes matrix calculated
% with the algorithm of Little and Chen, and implemented in the function
% MI_to_ZRES. Differs from MIP_to_ZPX because it keeps the diagonal values.

mat=input_mat;

[rows,cols]=size(mat);
mean_row=zeros(rows,1);
std_row=zeros(rows,1);
ZPX2_mat=zeros(rows,cols);

for i=1:rows
    mean_row(i)=nanmean(mat(i,:));
    std_row(i)=nanstd(mat(i,:));   
end
for i=1:rows
    for j=i:rows

    ZPX2_i=(mat(i,j)-mean_row(i))/std_row(i);
    ZPX2_j=(mat(i,j)-mean_row(j))/std_row(j);
    
    ZPX2_mat(i,j)=(ZPX2_i*ZPX2_j);

% Here we correct for the product of two negative ZPX2_i and ZPX2_J, which
% would give the wrong MI.
    
    if (ZPX2_i<0&&ZPX2_j<0)
        ZPX2_mat(i,j)=-ZPX2_mat(i,j);
    end

% Symmetrize.

    ZPX2_mat(j,i)=ZPX2_mat(i,j);
    end
    
    % ZPX2_mat(i,i)=NaN;

% Here we take the square root and then only the real part.

ZPX_mat=real(sqrt(ZPX2_mat));

end
end


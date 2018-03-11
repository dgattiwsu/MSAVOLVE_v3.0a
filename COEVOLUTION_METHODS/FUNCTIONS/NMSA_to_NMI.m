function [ MI_mat,NMI_mat ] = NMSA_to_NMI( nmsa )
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (.nmsa). The function uses  the "NormMutualInformation" function 
% from Dwinnel package.

[nrows,ncols]=size(nmsa);
MI_mat=zeros(ncols,ncols);
NMI_mat=zeros(ncols,ncols);
gap_mat=zeros(ncols,ncols);

% Here we start looping over all the column and calculate the MI. 
% For every pair of columns only the indices that do not have gaps 
% (represented by the number 25 in the nmsa matrix) in both columns  are 
% considered. Note that the gap matrix is initially all zero, but as we
% loop, it becomes populated with the number of rows in which both columns
% do not have gaps.

for i=1:ncols
for j=i:ncols
    gap=find((nmsa(:,i)~=25) & (nmsa(:,j)~=25));
    gap_mat(i,j)=numel(gap);    

    [MI_mat(i,j),NMI_mat(i,j)]=(NormMutualInformation(nmsa(gap,i),nmsa(gap,j)));

end
end

% Here we divide the gap matrix by the total number of sequences in the
% alignment. This means the scale matrix contains for every ij position a
% scaling factors corresponding to the fraction of rows in any ij pair that
% did not contained gaps.

scale_mat=gap_mat/nrows;

% Symmetrize both scaling and MI matrices

scale_mat=scale_mat+scale_mat';
MI_mat=MI_mat+MI_mat';
NMI_mat=NMI_mat+NMI_mat';

% Set the diagonal to the self values.

for i=1:ncols
    
MI_mat(i,i)=MI_mat(i,i)/2;
NMI_mat(i,i)=NMI_mat(i,i)/2;
scale_mat(i,i)=scale_mat(i,i)/2;

end


% Multiply point by point

MI_mat=MI_mat.*scale_mat;
NMI_mat=NMI_mat.*scale_mat;

% Better set the diagonal to NaN than to 0.

for i=1:ncols

MI_mat(i,i)=NaN;
NMI_mat(i,i)=NaN;

end

end


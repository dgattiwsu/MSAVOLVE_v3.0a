function [ SU_mat ] = NMSA_to_SU( nmsa )
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (.nmsa). The function uses  the "mutualinformation" or
% "symmetricuncertainty" functions of the Information Theory Toolbox v.1.0
% for Matlab (http://www.mathworks.com/matlabcentral/fileexchange/
% 17993-information-theory-toolbox-v1-0) developed by Joaquin Goni 
% (Dept. of Physics and Applied Mathematics, University of Navarra, 
% Pamplona, Spain), which can be downloaded from Matlabcentral.
%   
[nrows,ncols]=size(nmsa);
SU_mat=zeros(ncols,ncols);
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
    SU_mat(i,j)=(SymmetricUncertainty(nmsa(gap,i),nmsa(gap,j)));

% Use the following instead for unnormalized MI instead of SU.    
% SU_mat(i,j)=(MutualInformation(nmsa(gap,i),nmsa(gap,j)));

end
end

% Here we divide the gap matrix by the total number of sequences in the
% alignment. This means the scale matrix contains for every ij position a
% scaling factors corresponding to the fraction of rows in any ij pair that
% did not contained gaps.

scale_mat=gap_mat/nrows;

% Symmetrize both scaling and MI matrices

scale_mat=scale_mat+scale_mat';
SU_mat=SU_mat+SU_mat';

% Multiply point by point

SU_mat=SU_mat.*scale_mat;

% Better set the diagonal to NaN than to 0.

for i=1:ncols

SU_mat(i,i)=NaN;

end

end


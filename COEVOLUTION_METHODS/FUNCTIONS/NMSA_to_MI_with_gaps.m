function [ MI_mat ] = NMSA_to_MI_with_gaps( nmsa )
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (.nmsa). The function uses  the "MutualInformation" function from 
% Dwinnel package.

[nrows,ncols]=size(nmsa);
MI_mat=zeros(ncols,ncols);

% Here we start looping over all the column and calculate the MI. 

for i=1:ncols
for j=i:ncols

    MI_mat(i,j) = (MutualInformation(nmsa(:,i),nmsa(:,j)));
    % Alternatively we can use the Mi function from Goni's toolbox.
    % MI_mat(i,j) = (mutualinformation_goni(nmsa(:,i),nmsa(:,j)));
    MI_mat(j,i) = MI_mat(i,j);
    
end
end

% Better set the diagonal to NaN than to 0.

for i=1:ncols

MI_mat(i,i)=NaN;

end

end


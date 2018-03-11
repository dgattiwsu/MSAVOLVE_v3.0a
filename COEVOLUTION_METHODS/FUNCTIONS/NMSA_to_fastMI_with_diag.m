function [ MI_mat ] = NMSA_to_fastMI_with_diag( nmsa )
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (.nmsa). The function uses  the "fastMI" function. It does not
% apply a correction for gaps.

[~,ncols] = size(nmsa);
MI_mat = zeros(ncols,ncols);

for i = 1:ncols
    for j = i:ncols
    MI_mat(i,j) = fastMI(nmsa(:,i),nmsa(:,j));
    MI_mat(j,i) = MI_mat(i,j);    
    end
end

end


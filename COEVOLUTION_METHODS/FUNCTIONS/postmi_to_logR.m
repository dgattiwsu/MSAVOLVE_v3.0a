function [ bayesMI,bayesMI_covar_vec ] = postmi_to_bayesMI( postmi_mat,npos,fcov )
% This function returns the bayesMI matrices and a vector of the top
% fraction of covarying residues based on the value of fcov (eg., fcov = 10
% means 10% of all positions).

ncov = round(npos*fcov/100);

bayesMI = zeros(npos,npos);
[nrows,~] = size(postmi_mat);

for k = 1:nrows
    i = postmi_mat(k,1) + 1;
    j = postmi_mat(k,2) + 1;
    logP = postmi_mat(k,3);
    
    bayesMI(i,j) = logP;
end

bayesMI = bayesMI + bayesMI';
bayesMI = exp(bayesMI);

for i = 1:npos
    bayesMI(i,i) = 0;
end

sorted_bayesMI = sort_matrix_descend(bayesMI);
bayesMI_covar_vec = sorted_bayesMI(1:ncov,2:3); 

for i = 1:npos
    bayesMI(i,i) = NaN;
end

% bayesMI_hits = numel(intersect(covar_vec,bayesMI_covar_vec,'rows'))/2

end


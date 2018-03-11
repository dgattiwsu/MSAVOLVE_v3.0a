function [ DCA,DCA_covar_vec ] = plm_dca_to_wMI_DCA( dca,fcov )
% This function returns the MI and DCA matrices and a vector of the top
% fraction of covarying residues based on the value of fcov (eg., fcov = 10
% means 10% of all positions).
npos = max(dca(:,2));
ncov = round(npos*fcov/100);

wMI = zeros(npos,npos);
DCA = zeros(npos,npos);
[nrows,~] = size(dca);

for k = 1:nrows
    i = dca(k,1);
    j = dca(k,2);
    di = dca(k,3);
    
    DCA(i,j) = di;
end

clear wmi di

DCA = DCA + DCA';

sorted_DCA = sort_matrix_descend(DCA);
DCA_covar_vec = sorted_DCA(1:ncov,2:3); 

for i = 1:npos
    DCA(i,i) = NaN;
end

end


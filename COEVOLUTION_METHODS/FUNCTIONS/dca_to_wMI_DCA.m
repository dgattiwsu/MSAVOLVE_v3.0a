function [ wMI,DCA,wMI_covar_vec,DCA_covar_vec ] = dca_to_wMI_DCA( dca,fcov )
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
    wmi = dca(k,3);
    di = dca(k,4);
    
    wMI(i,j) = wmi;
    DCA(i,j) = di;
end

clear wmi di

wMI = wMI + wMI';
DCA = DCA + DCA';

sorted_wMI = sort_matrix_descend(wMI);
wMI_covar_vec = sorted_wMI(1:ncov,2:3);            

sorted_DCA = sort_matrix_descend(DCA);
DCA_covar_vec = sorted_DCA(1:ncov,2:3); 

for i = 1:npos
    wMI(i,i) = NaN;
    DCA(i,i) = NaN;
end

end


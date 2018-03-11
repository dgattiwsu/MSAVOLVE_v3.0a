function [ mrH ] = bin_rel_entropy( bnmsa )
%--------------------------------------------------------------------------
% relative MI calculation; the following are the rules used in generating 
% the arrays.
% pI(1) == 1 ; pI(2) == 0 

bnmsa = + bnmsa;
const1 = log2(nrows);
all = nrows*ncols;
sum_of_1 = sum(bnmsa(:));
sum_of_0 = all - sum_of_1;
const2(1,1) = log2(sum_of_1/all);
const2(2,1) = log2(sum_of_0/all);
rH = zeros(1,ncols);

pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        H(I) = (sum(pI(ind,I).*(log2(pI(ind,I) - const1 .- const2)))/nrows;
end

mrH = mean(H);

end


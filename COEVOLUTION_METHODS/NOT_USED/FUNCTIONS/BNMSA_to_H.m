function [pH] = BNMSA_to_H(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% Binary Entropy calculation.

bnmsa = + bnmsa;
const = log2(nrows);
H = zeros(1,ncols);
pH = zeros(1,2);
pMI = zeros(ncols,ncols);
pIJ = zeros(4,1);
epIJ = zeros(4,1);

pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
end

pH = mean(H);

end



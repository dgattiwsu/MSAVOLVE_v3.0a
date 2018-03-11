function [gH] = BNMSA_to_gH(bnmsa)
%--------------------------------------------------------------------------
% Binary Entropy calculation.
[nrows,ncols] = size(bnmsa);
bnmsa = + bnmsa(:);
nrows = nrows*ncols;
const = log2(nrows);

% Even faster
pI = zeros(2,1);
        pI(1) = sum(bnmsa);
        pI(2) = nrows - pI(1);
        ind = find(pI);        
        H = -(sum(pI(ind).*(log2(pI(ind)) - const)))/nrows;
gH = H;

end


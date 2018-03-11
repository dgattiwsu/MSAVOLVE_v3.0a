function [ mrH ] = BNMSA_to_mrH_2( bnmsa,nrows,ncols )
%--------------------------------------------------------------------------
% relative MI calculation; the following are the rules used in generating 
% the arrays.
% pI(1) == 1 ; pI(2) == 0 
% pJ(1) == 1 ; pJ(2) == 0

bnmsa = + bnmsa;
const1 = log2(nrows);
%all = nrows*ncols;
sum_of_1 = mean(sum(bnmsa))/nrows;
sum_of_0 = mean(sum(~bnmsa))/nrows;
const2(1,1) = log2(sum_of_1);
const2(2,1) = log2(sum_of_0);
rH = zeros(1,ncols);

pI = zeros(2,1);
for I = 1:ncols
        pI(1) = sum(bnmsa(:,I));
        pI(2) = nrows - pI(1);
        ind = find(pI);
        rH(I) = (sum(pI(ind).*(log2(pI(ind)) - const1 - const2(ind))))/nrows;
end
mrH = mean(rH);

% The previous code is equivalent to the following, but faster:

% bnmsa = + bnmsa;
% all = nrows*ncols;
% sum_of_1 = sum(bnmsa(:));
% sum_of_0 = all - sum_of_1;
% bg_pI(1,1) = sum_of_1/all;
% bg_pI(2,1) = sum_of_0/all;
% rH = zeros(1,ncols);
% 
% pI = zeros(2,1);
% for I = 1:ncols
%         pI(1) = sum(bnmsa(:,I));
%         pI(2) = nrows - pI(1);
%         pI(1) = pI(1)/nrows;
%         pI(2) = pI(2)/nrows;
%         
%         ind = find(pI(:));
%         rH(I) = (sum(pI(ind).*(log2(pI(ind)./bg_pI(ind)))));
% end
% mrH = mean(rH);


end


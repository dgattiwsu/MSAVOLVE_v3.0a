function [binmsa] = nmsa_to_v_binmsa_21q(nmsa)
% Returns each sequence of length L as 21 row vectors of size L with 0 and 1. 
% Number 21 represents gaps (which would be # 25 in the original
% Matlab numeric representation of an MSA.

[nseq,npos]=size(nmsa);
ind25 = nmsa == 25;
nmsa(ind25) = 21;
binmsa=zeros(21*nseq,npos);
for i=1:nseq 
    for aa=1:21 
        binmsa(21*(i-1)+aa,:)=(nmsa(i,:)==aa); 
    end; 
end;

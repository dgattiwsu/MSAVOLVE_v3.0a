function [binmsa]=nmsa_to_v_binmsa_20q(nmsa)
% Returns each sequence of length L as 20 row vectors of size L with 0 and 1. 
% Gaps (which would be # 25 in the original Matlab numeric representation 
% of an MSA are ignored. Thus, if at a certain position there is a gap that
% position will be converted into a column vector of 20 0s.

[nseq,npos]=size(nmsa);
binmsa=zeros(20*nseq,npos);
for i=1:nseq 
    for aa=1:20 
        binmsa(20*(i-1)+aa,:)=(nmsa(i,:)==aa); 
    end; 
end;

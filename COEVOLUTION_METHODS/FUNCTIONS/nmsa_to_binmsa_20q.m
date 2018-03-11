function [binmsa]=nmsa_to_binmsa_20q(nmsa)
% Returns each sequence of length L as a vector of size 20L with 0 and 1. 
% Gaps (which would be # 25 in the original Matlab numeric representation 
% of an MSA are ignored. Thus, if at a certain position there is a gap that
% position will be converted into a vector of 20 0s.

[nseq,npos]=size(nmsa);
binmsa=zeros(nseq,20*npos);
for i=1:npos 
    for aa=1:20 
        binmsa(:,20*(i-1)+aa)=(nmsa(:,i)==aa); 
    end; 
end;

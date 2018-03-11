function [binmsa]=nmsa_to_binmsa(nmsa)
% Returns each sequence of length L as a vector of size 20L with 0 and 1. 
% Numbers 21-25 represent various possible symbols (eg. gaps) in the 
% Matlab numeric representation of an MSA.

[nseq,npos]=size(nmsa); 
binmsa=zeros(nseq,25*npos);
for i=1:npos 
    for aa=1:25 
        binmsa(:,25*(i-1)+aa)=(nmsa(:,i)==aa); 
    end; 
end;

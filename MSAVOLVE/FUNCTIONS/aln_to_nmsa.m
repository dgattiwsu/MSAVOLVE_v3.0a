function [ msa, nmsa ] = aln_to_nmsa(aln)

% This function converts a multiple alignment in clustal format into a
% multiple alignment in matlab numeric format (nmsa). In invoking this
% function remember to place the name of the clustal file in quotes. For
% example to convert the msa "eco.aln" we would write:
%
% eco = aln_to_nmsa('eco.aln')
%

msa = multialignread(aln);

[~,nseq] = size(msa);
[~,nres] = size(msa(1,1).Sequence);

nmsa = zeros(nseq,nres);

for i=1:nseq
    
   nmsa(i,:)= aa2int(msa(1,i).Sequence);
  
end


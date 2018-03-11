function [ msa, nmsa ] = faln_to_nmsa(aln)

% This function converts a multiple alignment in fasta format into a
% multiple alignment in matlab numeric format (nmsa). In invoking this
% function remember to place the name of the fasta file in quotes. For
% example to convert the msa "eco.faln" we would write:
%
% eco = faln_to_nmsa('eco.faln')
%

msa = fastaread(aln);

[nseq,~] = size(msa);
[~,nres] = size(msa(1,1).Sequence);

nmsa = zeros(nseq,nres);

for i=1:nseq
    
   nmsa(i,:)= aa2int(msa(i,1).Sequence);
  
end


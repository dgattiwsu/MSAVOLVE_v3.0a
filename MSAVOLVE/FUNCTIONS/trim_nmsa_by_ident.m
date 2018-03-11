function [ trimmed_nmsa,nkeep,keep ] = trim_nmsa_by_ident( nmsa,threshold )
% This function removes from an nmsa all sequences whose identity is higher
% than a defined threshold, but keeps 1 example of each group.
% Usage: [trimmed_nmsa,nkeep,keep_index] = trim_nmsa_by_ident(nmsa,0.85)

[nseq,~] = size(nmsa);
keep = ones(nseq,1);
id_mat = zeros(nseq);
for i = 1:nseq
    for j = i:nseq;
    seq_1 = nmsa(i,:);
    seq_2 = nmsa(j,:);
    shared = (seq_1 ~= 25) | (seq_2 ~= 25); 
    ident = sum((seq_1(shared) == seq_2(shared))/sum(shared));
    id_mat(i,j) = ident;
    end
end
for i = 1:nseq
    out = find(id_mat(:,i)>threshold);
    out = setdiff(out,i);
    keep(out) = 0;
end
nkeep = sum(keep);
keep = logical(keep);
trimmed_nmsa = nmsa(keep,:);
end


function [trimmed_msa] = trim_msa_by_threshold(msa,threshold)
% This function trims an msa by keeping only the sequences that are less
% than threshold (e.g. 0.8 = 80%) identical to any other sequence. For
% every threshold used it provides the Meff of the original msa
% corresponding to that threshold.
[nrows,~] = size(msa);
ind1 = ones(nrows,1);

dist = seqpdist(msa,'UseParallel','true','SquareForm','true',...
                'ScoringMatrix','BLOSUM62');
dist = max(dist(:)) - dist;
dist_threshold = dist > threshold; 

W = 1./sum(dist_threshold)';
Meff=sum(W);
fprintf('Meff = %f \n', Meff);

for i = 1:nrows
    dist_threshold(i,i) = 0;
end

for i = 1:nrows
    [~,ind1(i)] = max(dist_threshold(i,:));
end
accept = sum(dist_threshold)';
ind2 = find(~accept);
discard = sum(triu(dist_threshold));
ind3 = find(discard);
ind = unique([ind1;ind2]);
ind = setdiff(ind,ind3);
trimmed_msa = msa(ind,:);

end

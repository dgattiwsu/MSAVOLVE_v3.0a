function [trimmed_nmsa,trimmed_smsa] = trim_nmsa_by_threshold(nmsa,smsa,threshold)
% This function trims an msa by keeping only the sequences that are less
% than threshold (e.g. 0.8 = 80%) identical to any other sequence. For
% every threshold used it provides the Meff of the original msa
% corresponding to that threshold.
[nrows,ncols] = size(nmsa);
W = ones(nrows,1);
ind1 = ones(nrows,1);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:))/ncols;
        dist(j,i) = dist(i,j);
    end
end
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
trimmed_nmsa = nmsa(ind,:);
trimmed_smsa = smsa(ind,1);

end

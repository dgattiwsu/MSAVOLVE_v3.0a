function [ fnmsa ] = nmsa_to_fnmsa(nmsa,threshold)
%
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

ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
% ind_23 = nmsa == 23;
% nmsa(ind_23) = 21;

loq = 0.5/21;
loq2 = 0.5/(21^2);
Fi = zeros(21,ncols);
nmsa_3 = zeros(nrows,ncols,21);

% Here we create a 3d matrix in which every layer has the dimensions of the 
% nmsa and represents one of 21 symbols. Then every row of each layer is scaled
% by the weight of that sequence.

for i = 1:21
    nmsa_3(:,:,i) = nmsa == i;
    for j = 1:nrows
        nmsa_3(j,:,i) = nmsa_3(j,:,i)*W(j);
    end
end

for i = 1:21
    layer = squeeze(nmsa_3(:,:,i));
    Fi(i,:) = sum(layer)/Meff;
    % With pseudocount.
    % Fi(i,:) = (sum(layer)+loq)/(loq+Meff);
end

% Here we create a new msa in which every symbol is replaced by the
% weigthed frequence of that symbol in that column of the msa.
fnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        fnmsa(i,j) = Fi(row(j),j);
        end
end

end


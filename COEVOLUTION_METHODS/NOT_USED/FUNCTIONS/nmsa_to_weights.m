function [W] = nmsa_to_weights(nmsa,threshold,nrows,ncols)
%--------------------------------------------------------------------------
% Weights calculation for the partition function
W = ones(nrows,1);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:))/ncols;
        dist(j,i) = dist(i,j);
    end
end
dist_threshold = dist > threshold; 
    W = 1./sum(dist_threshold)';
end

function [dist] = distance_matrix(nmsa)
[nrows,ncols] = size(nmsa);

% for i = 1:nrows
%     for j = i:nrows
%         dist(i,j) = sum(nmsa(i,:)~=nmsa(j,:));
%         dist(j,i) = dist(i,j);
%     end
% end

% The following method is much faster than the previous 'for' loop.
bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';
dist = ncols-dist;


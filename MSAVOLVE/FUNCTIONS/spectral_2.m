function [top_V,top_D_mat] = spectral_2(mat,nvec)
% Spectral analysis of a square matrix.
% This function returns the top 'nvec' eigenvectors and eigenvalues sorted
% in descending order. Similar to [V,D] = eigs(mat,nvec,'la') except for
% the signs of the vectors, which are fixed such that the mean of each 
% vector is positive. Eigenvalues are returned as a matrix.

[V,D]=eig(mat);
[sorted_D,ind]=sort(diag(D),'descend');
top_V = V(:,ind(1:nvec));
top_D = sorted_D(1:nvec);

top_D_mat = zeros(nvec);

for i = 1:nvec
top_D_mat(i,i) = top_D(i,1);
end

for i=1:nvec
    switch_sign = sign(mean(top_V(:,i)));
    top_V(:,i) = switch_sign*top_V(:,i);
end

end
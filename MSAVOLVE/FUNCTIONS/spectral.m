function [ sV,sD_mat ] = spectral( mat )

% Spectral analysis of a square matrix.
% This function returns all the eigenvectors and eigenvalues sorted
% in descending order. Similar to [V,D] = eig(mat) except for the ordering,
% which is descending instead of ascending, and the signs of the vectors,
% which are fixed such that the mean of each vector is positive.  
 
size_mat = size(mat,1);

[V,D]=eig(mat);

% We sort the eigenvalues in descending order, and then the eigenvectors
% accordingly.

[sD,ind]=sort(diag(D),'descend');

sD_mat = zeros(size_mat);

for i = 1:size_mat
sD_mat(i,i) = sD(i,1);
end

sV = V(:,ind);

% We fix the sign of each component of the vectors, such that the mean of  
% each vector is positive.

for i=1:size(sV,2)
    switch_sign = sign(mean(sV(:,i)));
    sV(:,i) = switch_sign*sV(:,i);
end

end
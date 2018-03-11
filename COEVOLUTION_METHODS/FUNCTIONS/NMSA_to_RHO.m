function [rho] = NMSA_to_RHO(nmsa)

[nrows,~] = size(nmsa);
binmsa = nmsa_to_binmsa_21q(nmsa);
cov_mat = cov(binmsa',1);
cov_inv = inv(cov_mat);
rho = zeros(nrows,nrows);

for i = 1:nrows
    for j = i:nrows
        rho(i,j) = -cov_inv(i,j)/sqrt(cov_inv(i,i)*cov_inv(j,j));
        rho(j,i) = rho(i,j);
    end
end

end


function [ corr1,corr2 ] = how_similar( A,B )
% This function tries to determine how similar are two matrices.
Afin = A;
Afin(isnan(A))=0;
Bfin = B;
Bfin(isnan(B))=0;
GE = eig(Afin,Bfin);
fin = isfinite(GE);
% corr1 = mean(GE(fin));
corr1 = norm(GE(fin));
EA = eig(Afin);
EB = eig(Bfin);
corr2 = corr(EA,EB);
end


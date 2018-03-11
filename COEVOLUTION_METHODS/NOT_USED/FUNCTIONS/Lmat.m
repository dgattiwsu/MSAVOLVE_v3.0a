ncols = 50;
Lmat = ones(ncols);
Lmat_scale = linspace(1,0,ncols);
Lmat = ones(ncols);
for i = 1:ncols
    vsize = ncols +1 -i;
    v = ones(vsize,1)*Lmat_scale(i);
    Lmat = Lmat + diag(v,i-1);
end
Lmat_u = triu(Lmat,1);
Lmat_l = tril(Lmat',0);
Lmat = Lmat_u + Lmat_l;
Lmat = Lmat -1;
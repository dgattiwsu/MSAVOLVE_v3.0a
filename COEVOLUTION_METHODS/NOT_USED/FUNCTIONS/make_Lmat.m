function [Lmat] = make_Lmat(ncols,method)

Lmat = ones(ncols);

switch method
    case 'LINEAR'
        Lmat_scale = linspace(0,1,ncols);
    case 'LOG'
        Lmat_scale = logspace(0,1,ncols);
end

for i = 1:ncols
    vsize = ncols +1 -i;
    v = ones(vsize,1)*Lmat_scale(i);
    Lmat = Lmat + diag(v,i-1);
end
Lmat_u = triu(Lmat,1);
Lmat_l = tril(Lmat',0);
Lmat = Lmat_u + Lmat_l;
Lmat = Lmat -1;
end
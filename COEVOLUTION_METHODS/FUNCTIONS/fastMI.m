function MI = fastMI(X,Y)
% Mutual information between two column vectors, X and Y, having
% integer values.
% 
% Modified from the original MI_GG function written by Giangregorio 
% Generoso (ggiangre@unisannio.it).

N = size(X,1);

const = log2(N);

h = accumarray([X Y], 1); 

xy_prod = sum(h,2)*sum(h,1);

xy_ind = h~=0;

MI = const + (sum(h(xy_ind) .* log2(h(xy_ind)./xy_prod(xy_ind))))/N;

end





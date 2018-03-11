function MI = fastMI_2(X,Y)
% Compute the mutual information between two column vectors, X and Y, having
% integer values. This function is a little slower than 'fastMI' but is 
% useful as a template for the fast multidimensional MI calculation in 
% 'NMSA_to_md2_ZPX2'.
% 
% Modified from an original implementation by GIANGREGORIO Generoso. 
% E-MAIL: ggiangre@unisannio.it
%__________________________________________________________________________

N = size(X,1);
const = log2(N);

% Joint histogram
h = accumarray([X Y], 1); 

x_marg = sum(h,2);
y_marg = sum(h,1);

x_marg = x_marg(x_marg~=0);
y_marg = y_marg(y_marg~=0);
xy_marg = h(h~=0);

% Entropy of X
Hx = const - sum(x_marg.*log2(x_marg))/N; 

% Entropy of Y
Hy = const - sum(y_marg.*log2(y_marg))/N;

% Joint Entropy of X and Y
Hxy = const - sum(xy_marg.*log2(xy_marg))/N; 

% Mutual Information
MI = Hx + Hy - Hxy; 

end


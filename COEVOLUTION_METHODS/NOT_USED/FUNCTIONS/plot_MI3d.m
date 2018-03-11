function [ X,Y,Z,W ] = plot_MI3d(mat)

W = mat(:);
[a b c] = size(mat);
lin_num = numel(W);
indices = 1:lin_num;
dims = [a b c];
[X Y Z] = ind2sub(dims,indices);
X = X';Y = Y';Z = Z';
wind = find(W);
W = W(wind);
X = X(wind);
Y = Y(wind);
Z = Z(wind);
scatter3(X,Y,Z,1000*W,W,'filled')

end


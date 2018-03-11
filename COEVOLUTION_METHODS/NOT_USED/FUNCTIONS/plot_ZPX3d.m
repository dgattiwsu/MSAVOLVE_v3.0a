function [ X,Y,Z,W ] = plot_ZPX3d(mat,flag,scale)
% If the flag is set to 0 every negative value is set to 0. If flag => 1,
% then all values are scaled in order to have a completely positive matrix.
% We recommend flag = 0, chosing appropriately the scale value to represent
% the dots.
if flag
min_mat = nanmin(mat(:));
mat = mat - min_mat;
else
neg_ind = mat < 0;
mat(neg_ind) = 0;
end
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
scatter3(X,Y,Z,scale*W,W,'filled')
box('on');

end


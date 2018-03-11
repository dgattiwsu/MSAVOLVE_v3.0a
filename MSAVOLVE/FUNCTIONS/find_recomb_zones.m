function [fitted_mat,model] = find_recomb_zones(mat,span,hrf)
% This function take a coevolution matrix and returns a model of a surface
% fitted (using the 'loess' algorithm) to the values of the matrix 
% (interpreted as Z values). It also returns a new matrix with the values 
% of the surface fit. 'span' is the window of matrix elements used in the
% loess fit. 'hrf' is the vector of recombination zones. Finally, the
% function plot for you the fitted 2D contour plot with the recombination
% zones shown as yellow dotted lines.
npos = size(mat,1);
X_vec = zeros((npos*npos),1);
Y_vec = zeros((npos*npos),1);
Z_vec = zeros((npos*npos),1);
ind = 0;
for i = 1:npos
for j = 1:npos
ind = ind+1;
X_vec(ind,1) = i;
Y_vec(ind,1) = j;
Z_vec(ind,1) = mat(i,j);
end
end

nanind = isnan(Z_vec);
X_vec = X_vec(~nanind);
Y_vec = Y_vec(~nanind);
Z_vec = Z_vec(~nanind);

% options = fitoptions('lowess'); % Lowess fit uses linear polynomial
options = fitoptions('loess'); % Loess fit uses quadratic polynomial
options.Span = span;
% options.Robust = 'Bisquare';
options.Robust = 'Off';
model = fit([X_vec,Y_vec],Z_vec,'loess',options);

[XI,YI] = meshgrid((1:npos),(1:npos));
fitted_mat = model(XI,YI);
imagesc(fitted_mat);
hold on

title('Loess fit of coevolution matrix','FontSize',14,'FontWeight','n');
        hrfy = [1 npos];
        for i = 1:length(hrf)
        hrfx = [hrf(i) hrf(i)];
        plot(hrfx,hrfy,':y','LineWidth',1.5)
        hold on
        end

        hrfx = [1 npos];
        for i = 1:length(hrf)
        hrfy = [hrf(i) hrf(i)];
        plot(hrfx,hrfy,':y','LineWidth',1.5)
        hold on
        end    
hold off
set(gca,'YDir','normal','YTick',zeros(1,0));    


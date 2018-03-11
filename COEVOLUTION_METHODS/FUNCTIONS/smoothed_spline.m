function [spline1] = smoothed_spline(x,y,tol_scale)
% LEAST_SQUARE_SPLINE  
%   Make sure the data are in rows ...
x = x(:).'; y = y(:).';
% ... and start by plotting the data specific to the highlighted spline fit.
tol = mean(diff(y))*tol_scale;
names={'data'};
% tol = (.5)^2*(2*pi);
spline1 = spaps(x, y, tol);

names{end+1} = 'spline1'; 
% fnplt(spline1,'-',2)




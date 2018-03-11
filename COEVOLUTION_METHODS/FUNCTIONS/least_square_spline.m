function [spline1] = least_square_spline(x,y)
% LEAST_SQUARE_SPLINE  
%   Make sure the data are in rows ...
x = x(:).'; y = y(:).';
% ... and start by plotting the data specific to the highlighted spline fit.

names={'data'};
% weights = ones(1,length(x));

% we are starting off with the least squares polynomial approximation of 
% order 4:
spline1 = spap2(1,4,x,y); 
% extract knots from current approximation:
knots = fnbrk(spline1,'knots'); 
% you changed the order:
knots = augknt(knots,14); 
% least-squares approximation:
spline1 = spap2(knots,14,x,y);
names{end+1} = 'spline1'; 
% fnplt(spline1,'-',2)




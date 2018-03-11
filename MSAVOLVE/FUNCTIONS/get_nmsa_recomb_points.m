function [cf_,x,nmsa_rel_entropy,y,crossover_points] = ...
    get_nmsa_recomb_points(nmsa,smoothing_param,delta)
% This function provides a hint of the position of the recombination
% points, based on a plot of the nmsa relative entropy.
% Recommended smoothing param between 0.995 and 1.0. Values closer to 1
% will give more crossover points.
% Recommended delta is between 0.05 and 0.01. Smaller values will give more
% crossover points

% First we get the relative entropy
[ nmsa_rel_entropy ] = get_nmsa_rel_entropy( nmsa );

% Then we create a fit
x = [1:numel(nmsa_rel_entropy)];

fo_ = fitoptions('method','SmoothingSpline','Normalize','on',...
    'SmoothingParam',smoothing_param);
ft_ = fittype('smoothingspline');

cf_ = fit(x',nmsa_rel_entropy',ft_,fo_);

% Here we get the y values of the fitted curve
y = cf_(x);

% Here we plot the fit.
plot(x,y);
set(gca,'Xlim',[1,max(x)]);

% Here we get an automatic assignment of the crossover points
[maxtab, mintab]=peakdet(y, delta, x);
crossover_points = maxtab(:,1)';

end


function [fit_model] = shape_pres_fit(X1,Y1)
ok_ = isfinite(X1) & isfinite(Y1);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs',...
        'Ignoring NaNs and Infs in data.' );
end
ft_ = fittype('pchipinterp');

% Fit this model using new data
fit_model = fit(X1(ok_),Y1(ok_),ft_);



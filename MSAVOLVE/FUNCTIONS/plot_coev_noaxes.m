function [ coev_plot ] = plot_coev_noaxes( coev_all,covar_vec )
%This function plots the coevolution matrix
    coev = nanmean(coev_all,3);
    finite_ind = isfinite(coev);
%    figure
    imagesc(coev,[0 max(coev(finite_ind))]);figure(gcf);
    set(gca,'YDir','normal')
    hold on
    scatter(covar_vec(:,1),covar_vec(:,2),...
    'YDataSource','covar_vec(:,2)',...
    'MarkerEdgeColor','yellow');figure(gcf)
    scatter(covar_vec(:,2),covar_vec(:,1),...
    'YDataSource','covar_vec(:,1)',...
    'MarkerEdgeColor','yellow');figure(gcf)
%    axis equal tight

end


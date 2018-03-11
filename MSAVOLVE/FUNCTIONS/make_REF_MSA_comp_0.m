%% Plots of various comparisons for the Reference MSA.
 
REF_cov_mat_lin = nonzeros(triu(REF_cov_mat,1));
    
REF_features = figure; clf; 
set(REF_features,'Units','normalized','Position',[0 0.2 0.6 0.8],'Name',...
    'REF Sequence Correlations');

figure(REF_features); 
subplot(2,2,1); 
hist(REF_cov_mat_lin,npos/2);
xlabel('Covariance','FontSize',14,'FontWeight','n'); 
ylabel('Number of Sequence Pairs','FontSize',14,'FontWeight','n'); 
title('REF covariance matrix histogram','FontSize',14,'FontWeight','n');
grid on

subplot(2,2,2);
o_REF_cov_mat = sort_cov_mat(REF_nmsa,REF_cov_mat);
imagesc(o_REF_cov_mat);figure(gcf);
xlabel('Seq No','FontSize',14,'FontWeight','n'); 
ylabel('Seq No','FontSize',14,'FontWeight','n'); 
title('REF covariance matrix heat map','FontSize',14,'FontWeight','n');
set(gca,'YDir','normal')

subplot(2,2,3); 
scatter3(REF_pc(:,1),REF_pc(:,2),REF_pc(:,3),'ok','SizeData',30,...
    'MarkerFaceColor','none');
xlabel('PC1','FontSize',14,'Fontweight','normal');
ylabel('PC2','FontSize',14,'Fontweight','normal');
zlabel('PC3','FontSize',14,'Fontweight','normal');
title('Covariance matrix spectral analysis',...
    'FontSize',14,'FontWeight','n');
% We recolor.
subplot(2,2,3); 
scatter3(REF_pc(REF_clust1,1),REF_pc(REF_clust1,2),REF_pc(REF_clust1,3),'ok',...
    'SizeData',30,'MarkerFaceColor','b');
hold on
scatter3(REF_pc(REF_clust2,1),REF_pc(REF_clust2,2),REF_pc(REF_clust2,3),'ok',...
    'SizeData',30,'MarkerFaceColor','r');
scatter3(REF_pc(REF_clust3,1),REF_pc(REF_clust3,2),REF_pc(REF_clust3,3),'ok',...
    'SizeData',30,'MarkerFaceColor','g');
scatter3(REF_pc(REF_clust4,1),REF_pc(REF_clust4,2),REF_pc(REF_clust4,3),'ok',...
    'SizeData',30,'MarkerFaceColor','c');
scatter3(REF_pc(REF_clust5,1),REF_pc(REF_clust5,2),REF_pc(REF_clust5,3),'ok',...
    'SizeData',30,'MarkerFaceColor','y');
scatter3(REF_pc(REF_clust6,1),REF_pc(REF_clust6,2),REF_pc(REF_clust6,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[0.0,0.5,1.0]);
scatter3(REF_pc(REF_clust7,1),REF_pc(REF_clust7,2),REF_pc(REF_clust7,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[1.0,0.5,0.0]);
scatter3(REF_pc(REF_clust8,1),REF_pc(REF_clust8,2),REF_pc(REF_clust8,3),'ok',...
    'SizeData',30,'MarkerFaceColor','m');
scatter3(REF_pc(REF_clust9,1),REF_pc(REF_clust9,2),REF_pc(REF_clust9,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[0.5,0.5,0.5]);
scatter3(REF_pc(REF_clust10,1),REF_pc(REF_clust10,2),REF_pc(REF_clust10,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[0.5,0.0,1.0]);
hold off
grid on; box on
xlabel('PC1','FontSize',14,'Fontweight','normal');
ylabel('PC2','FontSize',14,'Fontweight','normal');
zlabel('PC3','FontSize',14,'Fontweight','normal');
title('Covariance matrix spectral analysis',...
    'FontSize',14,'FontWeight','n');
grid on; box on



if isempty(covar_vec)
subplot(2,2,4); 
stairs(REF_rel_entropy,'DisplayName','REF_rel_entropy');
title('REF Rel. Entropy',...
    'FontSize',14,'FontWeight','n');
else    
covar_vec_size = size(covar_vec,1);
covar_pos_y1 = REF_rel_entropy(covar_vec(:,1));
covar_pos_y2 = REF_rel_entropy(covar_vec(:,2));
subplot(2,2,4); 
stairs(REF_rel_entropy,'DisplayName','REF_rel_entropy');
hold on
scatter(covar_vec(:,1),covar_pos_y1,30,'filled',...
    'MarkerEdgeColor','green');
scatter(covar_vec(:,2),covar_pos_y2,30,'filled',...
    'MarkerEdgeColor','green');
title('REF Rel. Entropy and Covarions Position',...
    'FontSize',14,'FontWeight','n');
end
xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
ylabel('Relative Entropy','FontSize',14,'FontWeight','n');
set(gca,'XLim',[1 npos]);
hold off

% end


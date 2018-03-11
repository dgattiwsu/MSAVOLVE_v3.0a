function [ MSA_Clusters ] = ...
    make_REF_MSA_comp_3(REF_nmsa,MSA_select,DISTANCE,DIMENSIONS,...
    MAX_CLUSTERS,REF_rel_entropy,MSA_rel_entropy,...
    covar_vec)

% This functions plots various comparisons between the Reference and
% Evolved MSAs.
REF_cmsa = int2aa(REF_nmsa);
c_MSA_select = int2aa(MSA_select);
 if DISTANCE
 REF_dist = seqpdist(REF_cmsa,'SquareForm',true,'Method','p-distance');
 REF_cov_mat = 1-REF_dist;
 MSA_dist = seqpdist(c_MSA_select,'SquareForm',true,'Method','p-distance');
 MSA_cov_mat = 1-MSA_dist;
 else
 REF_binmsa = nmsa_to_binmsa(REF_nmsa);
 REF_cov_mat = cov(REF_binmsa',1);  
 MSA_binmsa = nmsa_to_binmsa(MSA_select);
 MSA_cov_mat = cov(MSA_binmsa',1);  
 end
[REF_pc,REF_ev] = spectral(REF_cov_mat);
REF_Clusters = clusterdata(REF_pc(:,1:DIMENSIONS),'maxclust',MAX_CLUSTERS);
[MSA_pc,MSA_ev] = spectral(MSA_cov_mat);
MSA_Clusters = clusterdata(MSA_pc(:,1:DIMENSIONS),'maxclust',MAX_CLUSTERS);

% First we identify all the clusters.
REF_clust1 = find(REF_Clusters == 1);
REF_clust2 = find(REF_Clusters == 2);
REF_clust3 = find(REF_Clusters == 3);
REF_clust4 = find(REF_Clusters == 4);
REF_clust5 = find(REF_Clusters == 5);
REF_clust6 = find(REF_Clusters == 6);
MSA_clust1 = find(MSA_Clusters == 1);
MSA_clust2 = find(MSA_Clusters == 2);
MSA_clust3 = find(MSA_Clusters == 3);
MSA_clust4 = find(MSA_Clusters == 4);
MSA_clust5 = find(MSA_Clusters == 5);
MSA_clust6 = find(MSA_Clusters == 6);

npos = size(REF_rel_entropy,2);
 
REF_cov_mat_lin = nonzeros(triu(REF_cov_mat,1));
MSA_cov_mat_lin = nonzeros(triu(MSA_cov_mat,1));
    
REF_MSA_seq_sim = figure; clf; 
set(REF_MSA_seq_sim,'Units','normalized','Position',[0 0.3 1.0 0.7],'Name',...
    'REF Sequence Correlations');

figure(REF_MSA_seq_sim); 
subplot(2,4,1); 
hist(REF_cov_mat_lin,npos/2);
xlabel('Covariance','FontSize',14,'FontWeight','n'); 
ylabel('Number of Sequence Pairs','FontSize',14,'FontWeight','n'); 
title('REF covariance matrix histogram','FontSize',14,'FontWeight','n');
grid on

figure(REF_MSA_seq_sim); 
subplot(2,4,2);

o_REF_cov_mat = sort_cov_mat(REF_nmsa,REF_cov_mat);
imagesc(o_REF_cov_mat);figure(gcf);
xlabel('Seq No','FontSize',14,'FontWeight','n'); 
ylabel('Seq No','FontSize',14,'FontWeight','n'); 
title('REF covariance matrix heat map','FontSize',14,'FontWeight','n');
set(gca,'YDir','normal')

figure(REF_MSA_seq_sim);
subplot(2,4,3); 
scatter3(REF_pc(:,1),REF_pc(:,2),REF_pc(:,3),'ok','SizeData',30,...
    'MarkerFaceColor','none');
xlabel('PC1','FontSize',14,'Fontweight','normal');
ylabel('PC2','FontSize',14,'Fontweight','normal');
zlabel('PC3','FontSize',14,'Fontweight','normal');
title('Covariance matrix spectral analysis',...
    'FontSize',14,'FontWeight','n');
% We recolor.
figure(REF_MSA_seq_sim);
subplot(2,4,3); 
scatter3(REF_pc(REF_clust1,1),REF_pc(REF_clust1,2),REF_pc(REF_clust1,3),'ok',...
    'SizeData',30,'MarkerFaceColor','g');
hold on
scatter3(REF_pc(REF_clust2,1),REF_pc(REF_clust2,2),REF_pc(REF_clust2,3),'ok',...
    'SizeData',30,'MarkerFaceColor','m');
scatter3(REF_pc(REF_clust3,1),REF_pc(REF_clust3,2),REF_pc(REF_clust3,3),'ok',...
    'SizeData',30,'MarkerFaceColor','c');
scatter3(REF_pc(REF_clust4,1),REF_pc(REF_clust4,2),REF_pc(REF_clust4,3),'ok',...
    'SizeData',30,'MarkerFaceColor','r');
scatter3(REF_pc(REF_clust5,1),REF_pc(REF_clust5,2),REF_pc(REF_clust5,3),'ok',...
    'SizeData',30,'MarkerFaceColor','y');
scatter3(REF_pc(REF_clust6,1),REF_pc(REF_clust6,2),REF_pc(REF_clust6,3),'ok',...
    'SizeData',30,'MarkerFaceColor','b');
hold off
xlabel('PC1','FontSize',14,'Fontweight','normal');
ylabel('PC2','FontSize',14,'Fontweight','normal');
zlabel('PC3','FontSize',14,'Fontweight','normal');
title('Covariance matrix spectral analysis',...
    'FontSize',14,'FontWeight','n');
grid on; box on

covar_vec_size = size(covar_vec,1);
covar_pos_y1 = REF_rel_entropy(covar_vec(:,1));
covar_pos_y2 = REF_rel_entropy(covar_vec(:,2));

figure(REF_MSA_seq_sim); 
subplot(2,4,4); 
stairs(REF_rel_entropy,'DisplayName','REF_rel_entropy');
hold on
scatter(covar_vec(:,1),covar_pos_y1,30,'filled',...
    'MarkerEdgeColor','green');
scatter(covar_vec(:,2),covar_pos_y2,30,'filled',...
    'MarkerEdgeColor','green');
xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
ylabel('Relative Entropy','FontSize',14,'FontWeight','n');
title('REF Rel. Entropy and Covarions Position',...
    'FontSize',14,'FontWeight','n');
set(gca,'XLim',[1 npos]);
hold off

figure(REF_MSA_seq_sim); 
subplot(2,4,5); 
hist(MSA_cov_mat_lin,npos/2);
xlabel('Covariance','FontSize',14,'FontWeight','n'); 
ylabel('Number of Sequence Pairs','FontSize',14,'FontWeight','n'); 
title('MSA covariance matrix histogram','FontSize',14,'FontWeight','n');
grid on

figure(REF_MSA_seq_sim); 
subplot(2,4,6); 
o_MSA_cov_mat = sort_cov_mat(MSA_select,MSA_cov_mat);
imagesc(o_MSA_cov_mat);figure(gcf);
xlabel('Seq No','FontSize',14,'FontWeight','n'); 
ylabel('Seq No','FontSize',14,'FontWeight','n'); 
title('MSA covariance matrix heat map','FontSize',14,'FontWeight','n');
set(gca,'YDir','normal')

figure(REF_MSA_seq_sim);
subplot(2,4,7); 
scatter3(MSA_pc(:,1),MSA_pc(:,2),MSA_pc(:,3),'ok','SizeData',30,...
    'MarkerFaceColor','none');
xlabel('PC1','FontSize',14,'Fontweight','normal');
ylabel('PC2','FontSize',14,'Fontweight','normal');
zlabel('PC3','FontSize',14,'Fontweight','normal');
title('Covariance matrix spectral analysis',...
    'FontSize',14,'FontWeight','n');
% We recolor.
figure(REF_MSA_seq_sim);
subplot(2,4,7); 
scatter3(MSA_pc(MSA_clust1,1),MSA_pc(MSA_clust1,2),MSA_pc(MSA_clust1,3),'ok',...
    'SizeData',30,'MarkerFaceColor','g');
hold on
scatter3(MSA_pc(MSA_clust2,1),MSA_pc(MSA_clust2,2),MSA_pc(MSA_clust2,3),'ok',...
    'SizeData',30,'MarkerFaceColor','m');
scatter3(MSA_pc(MSA_clust3,1),MSA_pc(MSA_clust3,2),MSA_pc(MSA_clust3,3),'ok',...
    'SizeData',30,'MarkerFaceColor','c');
scatter3(MSA_pc(MSA_clust4,1),MSA_pc(MSA_clust4,2),MSA_pc(MSA_clust4,3),'ok',...
    'SizeData',30,'MarkerFaceColor','r');
scatter3(MSA_pc(MSA_clust5,1),MSA_pc(MSA_clust5,2),MSA_pc(MSA_clust5,3),'ok',...
    'SizeData',30,'MarkerFaceColor','y');
scatter3(MSA_pc(MSA_clust6,1),MSA_pc(MSA_clust6,2),MSA_pc(MSA_clust6,3),'ok',...
    'SizeData',30,'MarkerFaceColor','b');
hold off
xlabel('PC1','FontSize',14,'Fontweight','normal');
ylabel('PC2','FontSize',14,'Fontweight','normal');
zlabel('PC3','FontSize',14,'Fontweight','normal');
title('Covariance matrix spectral analysis',...
    'FontSize',14,'FontWeight','n');
grid on; box on

figure(REF_MSA_seq_sim);
subplot(2,4,8); 
stairs(REF_rel_entropy,'-b','DisplayName','REF_rel_entropy');
hold on
stairs(MSA_rel_entropy,'-r','DisplayName','msa_rel_entropy');
hold off
xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
ylabel('Relative Entropy','FontSize',14,'FontWeight','n');
title('REF and MSA Relative Entropy',...
    'FontSize',14,'FontWeight','n');
set(gca,'XLim',[1 npos]);

end


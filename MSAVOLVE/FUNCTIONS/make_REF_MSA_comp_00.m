%% Plots of various comparisons between Reference and Evolved MSAs.
[~,~,MSA_select_z] = size(MSA_select_ALL);
if n > MSA_select_z
    MSA_select = nanmean(MSA_select_ALL,3);
else
MSA_select = MSA_select_ALL(:,:,n);
c_MSA_select = int2aa(MSA_select_ALL(:,:,n));
end

 if DISTANCE == 1
 MSA_dist = seqpdist(c_MSA_select,'SquareForm',true,'Method','p-distance');
 MSA_cov_mat = 1-MSA_dist;
 elseif DISTANCE == 2
 MSA_dist = seqpdist(c_MSA_select,'SquareForm',true,'Method','alignment-score');
 MSA_cov_mat = 1-MSA_dist;
 elseif DISTANCE == 3
 MSA_dist = seqpdist(c_MSA_select,'SquareForm',true,'Method','Jukes-Cantor');
 MSA_cov_mat = (max(MSA_dist(:))-MSA_dist)/max(MSA_dist(:));
 else
 MSA_binmsa = nmsa_to_binmsa(MSA_select);
 MSA_cov_mat = cov(MSA_binmsa',1);  
 end

 [MSA_pc,~] = pcacov(MSA_cov_mat);
 
% We fix the sign of each component of the vectors, such that the mean of  
% each vector is positive.

for i=1:size(MSA_pc,2)
    switch_sign = sign(mean(MSA_pc(:,i)));
    MSA_pc(:,i) = switch_sign*MSA_pc(:,i);
end
 
% use fastICA if package is installed.
if ICA
    MSA_pc_t = MSA_pc(:,1:DIMENSIONS)';
    [MSA_ic,~,~] = fastica(MSA_pc_t);
    MSA_pc = MSA_ic';
end

MSA_Clusters = clusterdata(MSA_pc(:,1:DIMENSIONS),...
    'distance',cluster_distance,'linkage',cluster_linkage,'maxclust',MAX_CLUSTERS);

for i = 1:MAX_CLUSTERS
    v = ['MSA_clust' int2str(i)];
    eval([v ' = (b_ind2(i,1):b_ind2(i,2));']); 
end
for i = MAX_CLUSTERS+1:10
    v = ['MSA_clust' int2str(i)];
    eval([v ' = [];']); 
end

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
    'SizeData',30,'MarkerFaceColor','b');
hold on
scatter3(MSA_pc(MSA_clust2,1),MSA_pc(MSA_clust2,2),MSA_pc(MSA_clust2,3),'ok',...
    'SizeData',30,'MarkerFaceColor','r');
scatter3(MSA_pc(MSA_clust3,1),MSA_pc(MSA_clust3,2),MSA_pc(MSA_clust3,3),'ok',...
    'SizeData',30,'MarkerFaceColor','g');
scatter3(MSA_pc(MSA_clust4,1),MSA_pc(MSA_clust4,2),MSA_pc(MSA_clust4,3),'ok',...
    'SizeData',30,'MarkerFaceColor','c');
scatter3(MSA_pc(MSA_clust5,1),MSA_pc(MSA_clust5,2),MSA_pc(MSA_clust5,3),'ok',...
    'SizeData',30,'MarkerFaceColor','y');
scatter3(MSA_pc(MSA_clust6,1),MSA_pc(MSA_clust6,2),MSA_pc(MSA_clust6,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[0.0,0.5,1.0]);
scatter3(MSA_pc(MSA_clust7,1),MSA_pc(MSA_clust7,2),MSA_pc(MSA_clust7,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[1.0,0.5,0.0]);
scatter3(MSA_pc(MSA_clust8,1),MSA_pc(MSA_clust8,2),MSA_pc(MSA_clust8,3),'ok',...
    'SizeData',30,'MarkerFaceColor','m');
scatter3(MSA_pc(MSA_clust9,1),MSA_pc(MSA_clust9,2),MSA_pc(MSA_clust9,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[0.5,0.5,0.5]);
scatter3(MSA_pc(MSA_clust10,1),MSA_pc(MSA_clust10,2),MSA_pc(MSA_clust10,3),'ok',...
    'SizeData',30,'MarkerFaceColor',[0.5,0.0,1.0]);

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
stairs(MSA_rel_entropy(n,:),'-r','DisplayName','msa_rel_entropy');
hold off
xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
ylabel('Relative Entropy','FontSize',14,'FontWeight','n');
title('REF and MSA Relative Entropy',...
    'FontSize',14,'FontWeight','n');
set(gca,'XLim',[1 npos]);


function [ covar_vec ] = ...
    make_REF_MSA_comp_1(REF_nmsa,DISTANCE,DIMENSIONS,...
    MAX_CLUSTERS,REF_rel_entropy,covar_vec,TREE_ONLY,...
    distance_method,tree_method,neighjoin_method)

% This functions plots various comparisons between the Reference and
% Evolved MSAs.
REF_cmsa = int2aa(REF_nmsa);
 if DISTANCE
 REF_dist = seqpdist(REF_cmsa,'SquareForm',true,'Method','p-distance');
 REF_cov_mat = 1-REF_dist;
 else
 REF_binmsa = nmsa_to_binmsa(REF_nmsa);
 REF_cov_mat = cov(REF_binmsa',1);  
 end
[REF_pc,REF_ev] = spectral(REF_cov_mat);
REF_Clusters = clusterdata(REF_pc(:,1:DIMENSIONS),'maxclust',MAX_CLUSTERS);

if TREE_ONLY
    
    cmsa_dist = seqpdist(REF_cmsa,'Method',distance_method);
    if tree_method == 'seqneighjoin'
        cmsa_tree = seqneighjoin(cmsa_dist,neighjoin_method);
    else
        cmsa_tree = seqlinkage(cmsa_dist,tree_method);
    end
    [REF_Clusters,NODE_Clusters,~] = ...
     cluster(cmsa_tree,[],'criterion','gain','MaxClust',MAX_CLUSTERS);

%     h = plot(smsa_tree,'Orientation','left','TerminalLabels','false'); 
%     set(h.BranchLines(NODE_Clusters==1),'Color','b')
%     set(h.BranchLines(NODE_Clusters==2),'Color','r')
%     set(h.BranchLines(NODE_Clusters==3),'Color','g')
%     set(h.BranchLines(NODE_Clusters==4),'Color','c')
%     set(h.BranchLines(NODE_Clusters==5),'Color','y')
%     set(h.BranchLines(NODE_Clusters==6),'Color','k')
%     set(h.BranchLines(NODE_Clusters>6),'Color','m')
end

% First we identify all the clusters.
REF_clust1 = find(REF_Clusters == 1);
REF_clust2 = find(REF_Clusters == 2);
REF_clust3 = find(REF_Clusters == 3);
REF_clust4 = find(REF_Clusters == 4);
REF_clust5 = find(REF_Clusters == 5);
REF_clust6 = find(REF_Clusters == 6);

npos = size(REF_rel_entropy,2);
 
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

end


%% 'clustering' section

cluster_COV = 0;
cluster_MI = 0;
cluster_MIP = 0;
cluster_ZPX = 0;
cluster_ZPX2 = 0;
cluster_ZNMI = 0;
cluster_OMES = 0;
cluster_SCA = 0;
cluster_RAMA_SCA = 0;

if cluster_COV
    
% If we want to plot the clustergram with row normalization prior  
% to clustering we can standardize either the rows or the columns.
% The following loop standardizes the rows.
% for i = 1:npos
%    row_mean = nanmean(COV(i,:));
%    row_std =nanstd(COV(i,:));
%    zCOV(i,:) = (COV(i,:) - row_mean)/row_std;
% end
% zCOV_max = max(nantozero(zCOV(1:end)));
% zclCOV2 = nantozero(zCOV(clCOV_row_ind,clCOV_col_ind));
% imagesc(zclCOV2);figure(gcf);
% set(gca,'YDir','normal')
% However, there should not be a need to normalize the rows.

% To cluster only along the columns, producing clustered rows:
 COV_max = max(nantozero(COV(1:end)));
 cl_r_COV = clustergram(nantozero(COV),'Standardize',3,...
 'Cluster',1,'Symmetric','false','Colormap','jet','DisplayRange',COV_max);
% However, the input to clustergram appears to have a bug that forces row
% normalization. We remove manually this normalization.
 set(cl_r_COV,'Standardize',3);

 cl_r_COV_row_ind = str2num(cell2mat(cl_r_COV.RowLabels));
 cl_r_COV_col_ind = str2num(cell2mat(cl_r_COV.ColumnLabels'));
 cl_r_COV2 = nantozero(COV(cl_r_COV_row_ind,cl_r_COV_col_ind));
 COV_clustered_rows = figure; 
 set(COV_clustered_rows,'Units','normalized','Position',[0 0.3 0.5 0.7],...
    'Name','COV_clustered_rows'); clf;
 imagesc(cl_r_COV2);figure(gcf);
 axis equal tight
 set(gca,'YDir','normal')
 
% We can superimpose on the clustering figure the positions of the zones of
% recombination. In this case the zones in yellow will separate various 
% groups of columns.
 figure(gcf);
 hold on
 hrfy = [1 npos];
  for i = 1:length(hrf)
    hrf_col_ind(i) = find(cl_r_COV_col_ind == hrf(i));
    hrfx = [hrf_col_ind(i) hrf_col_ind(i)];
    plot(hrfx,hrfy,'--y','LineWidth',1)
    hold on
  end
  
% If we also superimpose the recombination zones in order to separate
% various groups of rows the result (red lines) is not obvious because in 
% this case the zones are not contiguous.  
 figure(gcf);
 hrfx = [1 npos];
 for i = 1:length(hrf)
    hrf_row_ind(i) = find(cl_r_COV_row_ind == hrf(i));
    hrfy = [hrf_row_ind(i) hrf_row_ind(i)];
    plot(hrfx,hrfy,'--r','LineWidth',1)
    hold on
 end    
 hold off
 
% To cluster first along the columns and then along the rows of the row
% clustered data, producing also clustered columns use instead:
 cl_rc_COV = clustergram(nantozero(COV),'Standardize',3,...
 'Cluster',3,'Symmetric','false','Colormap','jet','DisplayRange',COV_max);
 set(cl_rc_COV,'Standardize',3);

 cl_rc_COV_row_ind = str2num(cell2mat(cl_rc_COV.RowLabels));
 cl_rc_COV_col_ind = str2num(cell2mat(cl_rc_COV.ColumnLabels'));
 cl_rc_COV2 = nantozero(COV(cl_rc_COV_row_ind,cl_rc_COV_col_ind));
 COV_clustered_rows_and_columns = figure;
 set(COV_clustered_rows_and_columns,'Units','normalized',...
    'Position',[0 0.3 0.5 0.7],...
    'Name','COV_clustered_rows_and_columns'); clf;
 imagesc(cl_rc_COV2);figure(gcf);
 set(gca,'YDir','normal')

% We can superimpose on the clustering figure the positions of the zones of
% recombination. In this case the zones in yellow will separate various 
% groups of columns. The zones in red separate various groups of rows;  
% however, because of the previous clusering the resulting red lines do not
% represent contiguous zones.
 figure(gcf);
 hold on
 hrfy = [1 npos];
  for i = 1:length(hrf)
    hrf_col_ind(i) = find(cl_rc_COV_col_ind == hrf(i));
    hrfx = [hrf_col_ind(i) hrf_col_ind(i)];
    plot(hrfx,hrfy,'--y','LineWidth',1)
    hold on
  end
  
 figure(gcf);
 hrfx = [1 npos];
 for i = 1:length(hrf)
    hrf_row_ind(i) = find(cl_rc_COV_row_ind == hrf(i));
    hrfy = [hrf_row_ind(i) hrf_row_ind(i)];
    plot(hrfx,hrfy,'--r','LineWidth',1)
    hold on
 end    
 hold off
% We can also try a clustering as implemented in the SCA 4.5 package.
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/Dwinnel_MI_Package');
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME');
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/PROTEIN_EVOLUTION');
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA');
% addpath(genpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA4.5/sca45_release'));     

% COV_ats_mat = 1:1:npos;
% COV_ats = mat2cell(COV_ats_mat,[1],ones(1,npos));
 
% For the following you must have a license to SCA v4.5. 
% usage: [p,l,sort_order,sorted]=SCAcluster(matrix,pos,1.0,jet,1);
%****************SCAcluster.m*****************
%*********************************************
% Author: Rama Ranganathan (rama.ranganathan@UTSouthwestern.edu)
%
% Two dimensional hierarchical clustering of SCA correlation matrix using
% city-block distance metric and complete linkage. The function takes in the
% correlation matrix, the position labels, a max_scale for linear mapping
% of the color map to correlation values, the colormap, and a flag
% (raw_or_not) that determines whether an unclustered version of the matrix
% is kicked out as well. 
% Returns:
% 1. the distance for positions (pdist output, p) 
% 2. the clusters for positions(linkage output, l) 
% 3. the sorted indices for positions (sort_order)
% 4. figures of the clustered matrix
% 5. the position dendrogram, and 
% 6. the unclustered matrix (if desired)
%
% 07/2003 - initial
% 01/2005 - modified for SCA2.0 and R14S1
% 08/2008 - modified for SCA3.0/4.0
%
% Copyright R.Ranganathan 1999-2010
%*********************************************
%*********************************************
% [COV_p,COV_l,COV_sort_order,COV_sorted] = SCAcluster(COV_ntz,...
%     COV_ats,COV_max);
% If we want to use Bill Lane scacursor program to query the clustered
% matrix we must first scale it to a 0-2 range
% scaled_COV_sorted = 2*COV_sorted/COV_max;
% COV_sort_ats = COV_ats(COV_sort_order);
% COV_clustered_sca_style = figure;
% set(COV_clustered_sca_style,'Units','normalized',...
%    'Position',[0 0.3 0.5 0.7],...
%    'Name','COV_clustered_sca_style'); clf; 
% imagesc(scaled_COV_sorted,[0 2]);
% figure(gcf);
% scacursor(COV_sort_ats,COV_sort_ats,scaled_COV_sorted);
 
% Here we carry out a spectral analysis of the coevolution matrix.
% [COV_evec,COV_eval]=eigenvect(COV_ntz);

% 3-D plots of the top three eigenvectors
% COV_top_3_modes=figure; set(COV_top_3_modes,'Units','normalized',...
%    'Position',[0 0.3 0.5 0.8],'Name','Top 3 Eigenmodes'); clf; 
% scatter3(COV_evec(:,1),COV_evec(:,2),COV_evec(:,3),'ko','SizeData',...
%    50, 'MarkerFaceColor','b');
% xlabel('ev 1','FontSize',12,'FontWeight','b');
% ylabel('ev 2','FontSize',12,'FontWeight','b');
% zlabel('ev 3','FontSize',12,'FontWeight','b');
% hold on;
 % Next we add labels to identify the corresponding position in the msa 
 % for each coefficient.
% for i=1:npos;
%    text(COV_evec(i,1)+.005,COV_evec(i,2)+.005,COV_evec(i,3)+.005,COV_ats(i));
% end;
% hold off;
% grid on;
% box on;

% The first eigenmode is simply related to the average correlation value at
% each position:

% COV_evec1_corr=figure; 
% set(COV_evec1_corr,'Units','normalized','Position',[0.6 0.7 0.3 0.4],...
%    'Name','First COV Eigenvector'); clf; 
% scatter(COV_evec(:,1),mean(COV_ntz),'ko','SizeData',50,...
%    'MarkerFaceColor','b');
% grid on
% xlabel('ev 1','FontSize',12,'FontWeight','b');
% ylabel('mean COV position value','FontSize',12,'FontWeight','b');
 
% rmpath(genpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA4.5/sca45_release'));     
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/Dwinnel_MI_Package');
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME');
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/PROTEIN_EVOLUTION');
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if cluster_ZPX
% Here we cluster our style, first along the columns and then along the 
% rows of the row clustered data: 
 ZPX2_max = max(nantozero(ZPX2(1:end)));
 ZPX2_ntz = nantozero(ZPX2);
 cl_rc_ZPX2 = clustergram(nantozero(ZPX2),'Standardize',3,...
 'Cluster',3,'Symmetric','false','Colormap','jet','DisplayRange',ZPX2_max);
 set(cl_rc_ZPX2,'Standardize',3);

 cl_rc_ZPX2_row_ind = str2num(cell2mat(cl_rc_ZPX2.RowLabels));
 cl_rc_ZPX2_col_ind = str2num(cell2mat(cl_rc_ZPX2.ColumnLabels'));
 cl_rc_ZPX22 = nantozero(ZPX2(cl_rc_ZPX2_row_ind,cl_rc_ZPX2_col_ind));
 ZPX2_clustered_rows_and_columns = figure;
 set(ZPX2_clustered_rows_and_columns,'Units','normalized',...
    'Position',[0 0.3 0.5 0.7],...
    'Name','ZPX2_clustered_rows_and_columns'); clf;
 imagesc(cl_rc_ZPX22);figure(gcf);
 set(gca,'YDir','normal')
 axis equal tight

% Here we cluster SCA style.
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/Dwinnel_MI_Package');
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME');
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/PROTEIN_EVOLUTION');
% rmpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA');
% addpath(genpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA4.5/sca45_release'));     

% ZPX2_ats_mat = 1:1:npos;
% ZPX2_ats = mat2cell(ZPX2_ats_mat,[1],ones(1,npos));
 
% [ZPX2_p,ZPX2_l,ZPX2_sort_order,ZPX2_sorted] = ...
%     SCAcluster(ZPX2_ntz,ZPX2_ats,ZPX2_max);
% To use scacursor we scale the clustered to a 0-2 range
% scaled_ZPX2_sorted = 2*ZPX2_sorted/ZPX2_max;
% ZPX2_sort_ats = ZPX2_ats(ZPX2_sort_order);
% ZPX2_clustered_sca_style = figure;
% set(ZPX2_clustered_sca_style,'Units','normalized',...
%    'Position',[0 0.3 0.5 0.7],...
%    'Name','ZPX2_clustered_sca_style'); clf; 
% imagesc(scaled_ZPX2_sorted,[0 2]);
% figure(gcf);
% scacursor(ZPX2_sort_ats,ZPX2_sort_ats,scaled_ZPX2_sorted); 

% rmpath(genpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA4.5/sca45_release'));     
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/Dwinnel_MI_Package');
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME');
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/PROTEIN_EVOLUTION');
% addpath('/Users/mimo/Documents/KDO/Papers/KDO_PLoS_ONE_13/SECOND_SUBMISSION/MATLAB_NEIME/SCA');
end

if cluster_RAMA_SCA
% Here we cluster our style, first along the columns and then along the 
% rows of the row clustered data. The clustering is done for both the
% original and the noise filtered sca matrix.

 cl_rc_RAMA_simpleSCA = clustergram(nantozero(RAMA_simpleSCA),'Standardize',3,...
 'Cluster',3,'Symmetric','false','Colormap','jet','DisplayRange',...
 RAMA_simpleSCA_max);
 set(cl_rc_RAMA_simpleSCA,'Standardize',3);

 cl_rc_RAMA_simpleSCA_row_ind = str2num(cell2mat(cl_rc_RAMA_simpleSCA.RowLabels));
 cl_rc_RAMA_simpleSCA_col_ind = str2num(cell2mat(cl_rc_RAMA_simpleSCA.ColumnLabels'));
 cl_rc_RAMA_simpleSCA2 = nantozero(RAMA_simpleSCA(cl_rc_RAMA_simpleSCA_row_ind,...
     cl_rc_RAMA_simpleSCA_col_ind));
 RAMA_simpleSCA_clustered_rows_and_columns = figure;
 set(RAMA_simpleSCA_clustered_rows_and_columns,'Units','normalized',...
    'Position',[0 0.3 0.5 0.7],...
    'Name','RAMA_simpleSCA_clustered_rows_and_columns'); clf;
 imagesc(cl_rc_RAMA_simpleSCA2);figure(gcf);
 set(gca,'YDir','normal')
 axis equal tight


 cl_rc_RAMA_SCA = clustergram(nantozero(RAMA_SCA),'Standardize',3,...
 'Cluster',3,'Symmetric','false','Colormap','jet','DisplayRange',...
 RAMA_SCA_max);
 set(cl_rc_RAMA_SCA,'Standardize',3);

 cl_rc_RAMA_SCA_row_ind = str2num(cell2mat(cl_rc_RAMA_SCA.RowLabels));
 cl_rc_RAMA_SCA_col_ind = str2num(cell2mat(cl_rc_RAMA_SCA.ColumnLabels'));
 cl_rc_RAMA_SCA2 = nantozero(RAMA_SCA(cl_rc_RAMA_SCA_row_ind,...
     cl_rc_RAMA_SCA_col_ind));
 RAMA_SCA_clustered_rows_and_columns = figure;
 set(RAMA_SCA_clustered_rows_and_columns,'Units','normalized',...
    'Position',[0 0.3 0.5 0.7],...
    'Name','RAMA_SCA_clustered_rows_and_columns'); clf;
 imagesc(cl_rc_RAMA_SCA2);figure(gcf);
 set(gca,'YDir','normal')
 axis equal tight

% Here we cluster SCA style:
% RAMA_SCA_ats_mat = 1:1:npos;
% RAMA_SCA_ats = mat2cell(RAMA_SCA_ats_mat,[1],ones(1,npos));
 
% [RAMA_SCA_p,RAMA_SCA_l,RAMA_SCA_sort_order,RAMA_SCA_sorted] = ...
%     SCAcluster(RAMA_SCA,RAMA_SCA_ats,RAMA_SCA_max);
% To use scacursor we scale the clustered to a 0-2 range
% scaled_RAMA_SCA_sorted = 2*RAMA_SCA_sorted/RAMA_SCA_max;
% RAMA_SCA_sort_ats = RAMA_SCA_ats(RAMA_SCA_sort_order);
% RAMA_SCA_clustered_sca_style = figure;
% set(RAMA_SCA_clustered_sca_style,'Units','normalized',...
%    'Position',[0 0.3 0.5 0.7],...
%    'Name','RAMA_SCA_clustered_sca_style'); clf; 
% imagesc(scaled_RAMA_SCA_sorted,[0 2]);
% figure(gcf);
% scacursor(RAMA_SCA_sort_ats,RAMA_SCA_sort_ats,scaled_RAMA_SCA_sorted);     
end

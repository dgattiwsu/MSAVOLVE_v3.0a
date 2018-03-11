% This script can be run to correct the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW,gapW2,gapW3) can be applied. bayesMI, and SCA do
% not need a correction.
%%
[nrows,ncols] = size(MSA_select_ALL(:,:,1));

gW = zeros(ncols,ncols,end_cycle); 
gapW0 = ones(ncols,ncols);

for n = 1:end_cycle
    gapW = correct_coevmat_forgaps(MSA_select_ALL(:,:,n));
    gapW2 = gapW.^2;
    gapW3 = gapW.^3;
    gW(:,:,n) = gapW;
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 MI = MI_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_MI(n,:),cov_zscore_MI(n,:),...
     covar_vec_zscore_MI(:,:,n)] = ... 
     coev_stats_2(MI,ncov,covar_vec); 
  
end

 stat_fcov_MI = [nanmean(fcov_MI);nanstd(fcov_MI)];
 mean_cov_zscore_MI = nanmean(cov_zscore_MI,2);
 
for i = 1:end_cycle
 corr_cov_zscore_MI(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_MI(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_MI(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_MI(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_MI(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_MI(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 MIP = MIP_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_MIP(n,:),cov_zscore_MIP(n,:),...
     covar_vec_zscore_MIP(:,:,n)] = ... 
     coev_stats_2(MIP,ncov,covar_vec); 
  
end

 stat_fcov_MIP = [nanmean(fcov_MIP);nanstd(fcov_MIP)];
 mean_cov_zscore_MIP = nanmean(cov_zscore_MIP,2);
 
for i = 1:end_cycle
 corr_cov_zscore_MIP(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_MIP(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_MIP(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_MIP(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_MIP(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_MIP(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 ZPX2 = ZPX2_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_ZPX2(n,:),cov_zscore_ZPX2(n,:),...
     covar_vec_zscore_ZPX2(:,:,n)] = ... 
     coev_stats_2(ZPX2,ncov,covar_vec); 
  
end

 stat_fcov_ZPX2 = [nanmean(fcov_ZPX2);nanstd(fcov_ZPX2)];
 mean_cov_zscore_ZPX2 = nanmean(cov_zscore_ZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX2(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_ZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX2(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_ZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX2(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 nbZPX2 = nbZPX2_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_nbZPX2(n,:),cov_zscore_nbZPX2(n,:),...
     covar_vec_zscore_nbZPX2(:,:,n)] = ... 
     coev_stats_2(nbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_nbZPX2 = [nanmean(fcov_nbZPX2);nanstd(fcov_nbZPX2)];
 mean_cov_zscore_nbZPX2 = nanmean(cov_zscore_nbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_nbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_nbZPX2(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_nbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_nbZPX2(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_nbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_nbZPX2(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 gbZPX2 = gbZPX2_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_gbZPX2(n,:),cov_zscore_gbZPX2(n,:),...
     covar_vec_zscore_gbZPX2(:,:,n)] = ... 
     coev_stats_2(gbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_gbZPX2 = [nanmean(fcov_gbZPX2);nanstd(fcov_gbZPX2)];
 mean_cov_zscore_gbZPX2 = nanmean(cov_zscore_gbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_gbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_gbZPX2(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_gbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_gbZPX2(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_gbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_gbZPX2(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 dbZPX2 = dbZPX2_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_dbZPX2(n,:),cov_zscore_dbZPX2(n,:),...
     covar_vec_zscore_dbZPX2(:,:,n)] = ... 
     coev_stats_2(dbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_dbZPX2 = [nanmean(fcov_dbZPX2);nanstd(fcov_dbZPX2)];
 mean_cov_zscore_dbZPX2 = nanmean(cov_zscore_dbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_dbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_dbZPX2(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_dbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_dbZPX2(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_dbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_dbZPX2(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 fgbZPX2 = fgbZPX2_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_fgbZPX2(n,:),cov_zscore_fgbZPX2(n,:),...
     covar_vec_zscore_fgbZPX2(:,:,n)] = ... 
     coev_stats_2(fgbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_fgbZPX2 = [nanmean(fcov_fgbZPX2);nanstd(fcov_fgbZPX2)];
 mean_cov_zscore_fgbZPX2 = nanmean(cov_zscore_fgbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fgbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fgbZPX2(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_fgbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fgbZPX2(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_fgbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fgbZPX2(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 dgbZPX2 = dgbZPX2_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_dgbZPX2(n,:),cov_zscore_dgbZPX2(n,:),...
     covar_vec_zscore_dgbZPX2(:,:,n)] = ... 
     coev_stats_2(dgbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_dgbZPX2 = [nanmean(fcov_dgbZPX2);nanstd(fcov_dgbZPX2)];
 mean_cov_zscore_dgbZPX2 = nanmean(cov_zscore_dgbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_dgbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_dgbZPX2(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_dgbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_dgbZPX2(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_dgbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_dgbZPX2(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 DCA = DCA_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_DCA(n,:),cov_zscore_DCA(n,:),...
     covar_vec_zscore_DCA(:,:,n)] = ... 
     coev_stats_2(DCA,ncov,covar_vec); 
  
end

 stat_fcov_DCA = [nanmean(fcov_DCA);nanstd(fcov_DCA)];
 mean_cov_zscore_DCA = nanmean(cov_zscore_DCA,2);
 
for i = 1:end_cycle
 corr_cov_zscore_DCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_DCA(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_DCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_DCA(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_DCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_DCA(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 fodorMI = fodorMI_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_fodorMI(n,:),cov_zscore_fodorMI(n,:),...
     covar_vec_zscore_fodorMI(:,:,n)] = ... 
     coev_stats_2(fodorMI,ncov,covar_vec); 
  
end

 stat_fcov_fodorMI = [nanmean(fcov_fodorMI);nanstd(fcov_fodorMI)];
 mean_cov_zscore_fodorMI = nanmean(cov_zscore_fodorMI,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fodorMI(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorMI(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_fodorMI(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorMI(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_fodorMI(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorMI(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 fodorSCA = fodorSCA_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_fodorSCA(n,:),cov_zscore_fodorSCA(n,:),...
     covar_vec_zscore_fodorSCA(:,:,n)] = ... 
     coev_stats_2(fodorSCA,ncov,covar_vec); 
  
end

 stat_fcov_fodorSCA = [nanmean(fcov_fodorSCA);nanstd(fcov_fodorSCA)];
 mean_cov_zscore_fodorSCA = nanmean(cov_zscore_fodorSCA,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fodorSCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorSCA(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_fodorSCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorSCA(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_fodorSCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorSCA(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 OMES = OMES_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_OMES(n,:),cov_zscore_OMES(n,:),...
     covar_vec_zscore_OMES(:,:,n)] = ... 
     coev_stats_2(OMES,ncov,covar_vec); 
  
end

 stat_fcov_OMES = [nanmean(fcov_OMES);nanstd(fcov_OMES)];
 mean_cov_zscore_OMES = nanmean(cov_zscore_OMES,2);
 
for i = 1:end_cycle
 corr_cov_zscore_OMES(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_OMES(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_OMES(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_OMES(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_OMES(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_OMES(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 McBASC = McBASC_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_McBASC(n,:),cov_zscore_McBASC(n,:),...
     covar_vec_zscore_McBASC(:,:,n)] = ... 
     coev_stats_2(McBASC,ncov,covar_vec); 
  
end

 stat_fcov_McBASC = [nanmean(fcov_McBASC);nanstd(fcov_McBASC)];
 mean_cov_zscore_McBASC = nanmean(cov_zscore_McBASC,2);
 
for i = 1:end_cycle
 corr_cov_zscore_McBASC(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_McBASC(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_McBASC(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_McBASC(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_McBASC(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_McBASC(i,:)','type','Kendall','rows','pairwise');
end

%%
for n = 1:end_cycle

 % Here we calculate some statistics after correcting for gaps.
 ELSC = ELSC_ALL(:,:,n).*gW(:,:,n);
 
 [fcov_ELSC(n,:),cov_zscore_ELSC(n,:),...
     covar_vec_zscore_ELSC(:,:,n)] = ... 
     coev_stats_2(ELSC,ncov,covar_vec); 
  
end

 stat_fcov_ELSC = [nanmean(fcov_ELSC);nanstd(fcov_ELSC)];
 mean_cov_zscore_ELSC = nanmean(cov_zscore_ELSC,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ELSC(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ELSC(i,:)','type','Pearson','rows','pairwise');
 corr_cov_zscore_ELSC(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ELSC(i,:)','type','Spearman','rows','pairwise');
 corr_cov_zscore_ELSC(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ELSC(i,:)','type','Kendall','rows','pairwise');
end


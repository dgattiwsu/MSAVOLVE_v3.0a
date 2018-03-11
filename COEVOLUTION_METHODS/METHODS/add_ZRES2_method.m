% This script is used to add on the fly the ZRES2 method. This method is
% very similar to Chen's ZRES, but is calculated using Gloor's ZPX2
% algorithm.
% MI must be calculated first.

corr_cov_zscore_ZRES2 = zeros(end_cycle,3);
ZRES2_ALL = zeros(npos,npos,end_cycle);
fcov_ZRES2 = zeros(end_cycle,3);
cov_zscore_ZRES2 = zeros(end_cycle,ncov);
covar_vec_zscore_ZRES2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('ZRES2 cycle %d \n', n);

[~,ZRES2] = MI_to_ZPX(MI_ALL(:,:,n));

 % Here we calculate some statistics.
 ZRES2_ALL(:,:,n) = ZRES2;
 
 [fcov_ZRES2(n,:),cov_zscore_ZRES2(n,:),...
     covar_vec_zscore_ZRES2(:,:,n)] = ... 
     coev_stats_2(ZRES2,ncov,covar_vec); 
  
end

 stat_fcov_ZRES2 = [mean(fcov_ZRES2);std(fcov_ZRES2)];
 mean_cov_zscore_ZRES2 = nanmean(cov_zscore_ZRES2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ZRES2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES2(i,:)','type','Pearson');
 corr_cov_zscore_ZRES2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES2(i,:)','type','Spearman');
 corr_cov_zscore_ZRES2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES2(i,:)','type','Kendall');
end

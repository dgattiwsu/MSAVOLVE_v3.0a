% This script is used to add on the fly the ZPX method. 
% MIP must be calculated first.

corr_cov_zscore_ZPX = zeros(end_cycle,3);
ZPX_ALL = zeros(npos,npos,end_cycle);
fcov_ZPX = zeros(end_cycle,3);
cov_zscore_ZPX = zeros(end_cycle,ncov);
covar_vec_zscore_ZPX = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('ZPX cycle %d \n', n);

[ZPX,~] = MI_to_ZPX(MIP_ALL(:,:,n));

 % Here we calculate some statistics.
 ZPX_ALL(:,:,n) = ZPX;
 
 [fcov_ZPX(n,:),cov_zscore_ZPX(n,:),...
     covar_vec_zscore_ZPX(:,:,n)] = ... 
     coev_stats_2(ZPX,ncov,covar_vec); 
  
end

 stat_fcov_ZPX = [mean(fcov_ZPX);std(fcov_ZPX)];
 mean_cov_zscore_ZPX = nanmean(cov_zscore_ZPX,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ZPX(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX(i,:)','type','Pearson');
 corr_cov_zscore_ZPX(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX(i,:)','type','Spearman');
 corr_cov_zscore_ZPX(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX(i,:)','type','Kendall');
end

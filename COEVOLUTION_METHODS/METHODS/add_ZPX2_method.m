% This script is used to add on the fly the ZPX2 method. 
% MIP must be calculated first.

corr_cov_zscore_ZPX2 = zeros(end_cycle,3);
ZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_ZPX2 = zeros(end_cycle,3);
cov_zscore_ZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_ZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('ZPX2 cycle %d \n', n);

[~,ZPX2] = MI_to_ZPX(MIP_ALL(:,:,n));

 % Here we calculate some statistics.
 ZPX2_ALL(:,:,n) = ZPX2;
 
 [fcov_ZPX2(n,:),cov_zscore_ZPX2(n,:),...
     covar_vec_zscore_ZPX2(:,:,n)] = ... 
     coev_stats_2(ZPX2,ncov,covar_vec); 
  
end

 stat_fcov_ZPX2 = [mean(fcov_ZPX2);std(fcov_ZPX2)];
 mean_cov_zscore_ZPX2 = nanmean(cov_zscore_ZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX2(i,:)','type','Pearson');
 corr_cov_zscore_ZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX2(i,:)','type','Spearman');
 corr_cov_zscore_ZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ZPX2(i,:)','type','Kendall');
end

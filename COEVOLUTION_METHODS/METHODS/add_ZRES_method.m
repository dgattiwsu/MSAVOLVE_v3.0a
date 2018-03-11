% This script is used to add on the fly Chen's ZRES method.
% MI must be calculated first.

corr_cov_zscore_ZRES = zeros(end_cycle,3);
ZRES_ALL = zeros(npos,npos,end_cycle);
fcov_ZRES = zeros(end_cycle,3);
cov_zscore_ZRES = zeros(end_cycle,ncov);
covar_vec_zscore_ZRES = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('ZRES cycle %d \n', n);

[ZRES] = MI_to_ZRES(MI_ALL(:,:,n));

 % Here we calculate some statistics.
 ZRES_ALL(:,:,n) = ZRES;
 
 [fcov_ZRES(n,:),cov_zscore_ZRES(n,:),...
     covar_vec_zscore_ZRES(:,:,n)] = ... 
     coev_stats_2(ZRES,ncov,covar_vec); 
  
end

 stat_fcov_ZRES = [mean(fcov_ZRES);std(fcov_ZRES)];
 mean_cov_zscore_ZRES = nanmean(cov_zscore_ZRES,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ZRES(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES(i,:)','type','Pearson');
 corr_cov_zscore_ZRES(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES(i,:)','type','Spearman');
 corr_cov_zscore_ZRES(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES(i,:)','type','Kendall');
end

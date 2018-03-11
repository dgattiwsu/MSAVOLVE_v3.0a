% This script is used to add on the fly the ngcMI method (MI without gap
% correction). 

corr_cov_zscore_MI = zeros(end_cycle,3);
MI_ALL = zeros(npos,npos,end_cycle);
fcov_MI = zeros(end_cycle,3);
cov_zscore_MI = zeros(end_cycle,ncov);
covar_vec_zscore_MI = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('MI cycle %d \n', n);

[MI] = NMSA_to_ngcMI(MSA_select_ALL(:,:,n));

 % Here we calculate some statistics.
 MI_ALL(:,:,n) = MI;
 
 [fcov_MI(n,:),cov_zscore_MI(n,:),...
     covar_vec_zscore_MI(:,:,n)] = ... 
     coev_stats_2(MI,ncov,covar_vec); 
  
end

 stat_fcov_MI = [mean(fcov_MI);std(fcov_MI)];
 mean_cov_zscore_MI = mean(cov_zscore_MI,2);
 
for i = 1:end_cycle
 corr_cov_zscore_MI(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_MI(i,:)','type','Pearson');
 corr_cov_zscore_MI(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_MI(i,:)','type','Spearman');
 corr_cov_zscore_MI(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_MI(i,:)','type','Kendall');
end

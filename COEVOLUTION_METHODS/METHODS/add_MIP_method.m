% This script is used to add on the fly the MIP method. 
% MI must be calculated first.

corr_cov_zscore_MIP = zeros(end_cycle,3);
MIP_ALL = zeros(npos,npos,end_cycle);
fcov_MIP = zeros(end_cycle,3);
cov_zscore_MIP = zeros(end_cycle,ncov);
covar_vec_zscore_MIP = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('MIP cycle %d \n', n);

[MIP] = MI_to_MIP(MI_ALL(:,:,n));

 % Here we calculate some statistics.
 MIP_ALL(:,:,n) = MIP;
 
 [fcov_MIP(n,:),cov_zscore_MIP(n,:),...
     covar_vec_zscore_MIP(:,:,n)] = ... 
     coev_stats_2(MIP,ncov,covar_vec); 
  
end

 stat_fcov_MIP = [mean(fcov_MIP);std(fcov_MIP)];
 mean_cov_zscore_MIP = nanmean(cov_zscore_MIP,2);
 
for i = 1:end_cycle
 corr_cov_zscore_MIP(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_MIP(i,:)','type','Pearson');
 corr_cov_zscore_MIP(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_MIP(i,:)','type','Spearman');
 corr_cov_zscore_MIP(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_MIP(i,:)','type','Kendall');
end

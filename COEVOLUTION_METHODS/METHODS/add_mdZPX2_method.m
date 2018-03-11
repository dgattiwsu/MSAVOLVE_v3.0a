% This script is used to add on the fly the mdZPX2 method. 
% MIP must be calculated first.

corr_cov_zscore_mdZPX2 = zeros(end_cycle,3);
mdZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_mdZPX2 = zeros(end_cycle,3);
cov_zscore_mdZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_mdZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('mdZPX2 cycle %d \n', n);

[~,mdZPX2] = NMSA_to_mdZPX2(MSA_select(:,:,n),3);

 % Here we calculate some statistics.
 mdZPX2_ALL(:,:,n) = mdZPX2;
 
 [fcov_mdZPX2(n,:),cov_zscore_mdZPX2(n,:),...
     covar_vec_zscore_mdZPX2(:,:,n)] = ... 
     coev_stats_2(mdZPX2,ncov,covar_vec); 
  
end

 stat_fcov_mdZPX2 = [mean(fcov_mdZPX2);std(fcov_mdZPX2)];
 mean_cov_zscore_mdZPX2 = nanmean(cov_zscore_mdZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_mdZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_mdZPX2(i,:)','type','Pearson');
 corr_cov_zscore_mdZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_mdZPX2(i,:)','type','Spearman');
 corr_cov_zscore_mdZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_mdZPX2(i,:)','type','Kendall');
end

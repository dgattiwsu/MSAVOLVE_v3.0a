% This script is used to add on the fly the vZPX2 method. 
% MIP must be calculated first.

corr_cov_zscore_vZPX2 = zeros(end_cycle,3);
vZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_vZPX2 = zeros(end_cycle,3);
cov_zscore_vZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_vZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('vZPX2 cycle %d \n', n);

[~,vZPX2] = NMSA_to_vZPX2(MSA_select(:,:,n),3);

 % Here we calculate some statistics.
 vZPX2_ALL(:,:,n) = vZPX2;
 
 [fcov_vZPX2(n,:),cov_zscore_vZPX2(n,:),...
     covar_vec_zscore_vZPX2(:,:,n)] = ... 
     coev_stats_2(vZPX2,ncov,covar_vec); 
  
end

 stat_fcov_vZPX2 = [mean(fcov_vZPX2);std(fcov_vZPX2)];
 mean_cov_zscore_vZPX2 = nanmean(cov_zscore_vZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_vZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_vZPX2(i,:)','type','Pearson');
 corr_cov_zscore_vZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_vZPX2(i,:)','type','Spearman');
 corr_cov_zscore_vZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_vZPX2(i,:)','type','Kendall');
end

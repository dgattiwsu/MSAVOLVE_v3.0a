% This script is used to add on the fly the gbZPX2 method. 

corr_cov_zscore_gbZPX2 = zeros(end_cycle,3);
gbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_gbZPX2 = zeros(end_cycle,3);
cov_zscore_gbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_gbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('gbZPX2 cycle %d \n', n);

[gbZPX2] = NMSA_to_gbZPX2(MSA_select_ALL(:,:,n));

% Here we calculate the usual statistics gbZPX2.
 gbZPX2_ALL(:,:,n) = gbZPX2;
 
 [fcov_gbZPX2(n,:),cov_zscore_gbZPX2(n,:),...
     covar_vec_zscore_gbZPX2(:,:,n)] = ... 
     coev_stats_2(gbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_gbZPX2 = [mean(fcov_gbZPX2);std(fcov_gbZPX2)];
 mean_cov_zscore_gbZPX2 = nanmean(cov_zscore_gbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_gbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_gbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_gbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_gbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_gbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_gbZPX2(i,:)','type','Kendall');
end

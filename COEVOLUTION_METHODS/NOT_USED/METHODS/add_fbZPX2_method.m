% This script is used to add on the fly the fbZPX2 method. 

corr_cov_zscore_fbZPX2 = zeros(end_cycle,3);
fbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_fbZPX2 = zeros(end_cycle,3);
cov_zscore_fbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_fbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('fbZPX2 cycle %d \n', n);

[fbZPX2] = NMSA_to_fbZPX2(MSA_select_ALL(:,:,n));

 % Here we calculate the usual statistics fbZPX2.
 fbZPX2_ALL(:,:,n) = fbZPX2;
 
 [fcov_fbZPX2(n,:),cov_zscore_fbZPX2(n,:),...
     covar_vec_zscore_fbZPX2(:,:,n)] = ... 
     coev_stats_2(fbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_fbZPX2 = [mean(fcov_fbZPX2);std(fcov_fbZPX2)];
 mean_cov_zscore_fbZPX2 = nanmean(cov_zscore_fbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_fbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_fbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fbZPX2(i,:)','type','Kendall');
end

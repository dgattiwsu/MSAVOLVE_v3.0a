% This script is used to add on the fly the fdgbZPX2 method. 

corr_cov_zscore_fdgbZPX2 = zeros(end_cycle,3);
fdgbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_fdgbZPX2 = zeros(end_cycle,3);
cov_zscore_fdgbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_fdgbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('fdgbZPX2 cycle %d \n', n);

[fdgbZPX2] = NMSA_to_fdgbZPX2(MSA_select_ALL(:,:,n));

% Here we calculate the usual statistics fdgbZPX2.
 fdgbZPX2_ALL(:,:,n) = fdgbZPX2;
 
 [fcov_fdgbZPX2(n,:),cov_zscore_fdgbZPX2(n,:),...
     covar_vec_zscore_fdgbZPX2(:,:,n)] = ... 
     coev_stats_2(fdgbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_fdgbZPX2 = [mean(fcov_fdgbZPX2);std(fcov_fdgbZPX2)];
 mean_cov_zscore_fdgbZPX2 = nanmean(cov_zscore_fdgbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fdgbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fdgbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_fdgbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fdgbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_fdgbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fdgbZPX2(i,:)','type','Kendall');
end

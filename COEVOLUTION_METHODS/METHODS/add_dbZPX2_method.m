% This script is used to add on the fly the dbZPX2 method.

corr_cov_zscore_dbZPX2 = zeros(end_cycle,3);
dbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_dbZPX2 = zeros(end_cycle,3);
cov_zscore_dbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_dbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('dbZPX2 cycle %d \n', n);

[dbZPX2] = NMSA_to_dbZPX2(MSA_select_ALL(:,:,n));

% Here we calculate the usual statistics dbZPX2.
 dbZPX2_ALL(:,:,n) = dbZPX2;
 
 [fcov_dbZPX2(n,:),cov_zscore_dbZPX2(n,:),...
     covar_vec_zscore_dbZPX2(:,:,n)] = ... 
     coev_stats_2(dbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_dbZPX2 = [mean(fcov_dbZPX2);std(fcov_dbZPX2)];
 mean_cov_zscore_dbZPX2 = nanmean(cov_zscore_dbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_dbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_dbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_dbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_dbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_dbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_dbZPX2(i,:)','type','Kendall');
end

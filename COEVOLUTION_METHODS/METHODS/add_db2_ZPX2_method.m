% This script is used to add on the fly the db2_ZPX2 method.

corr_cov_zscore_db2_ZPX2 = zeros(end_cycle,3);
db2_ZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_db2_ZPX2 = zeros(end_cycle,3);
cov_zscore_db2_ZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_db2_ZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('db2_ZPX2 cycle %d \n', n);

[db2_ZPX2] = NMSA_to_db2_ZPX2(MSA_select_ALL(:,:,n),'NOGAPS','DISTANCE','EXCLUDE');

% Here we calculate the usual statistics db2_ZPX2.
 db2_ZPX2_ALL(:,:,n) = db2_ZPX2;
 
 [fcov_db2_ZPX2(n,:),cov_zscore_db2_ZPX2(n,:),...
     covar_vec_zscore_db2_ZPX2(:,:,n)] = ... 
     coev_stats_2(db2_ZPX2,ncov,covar_vec); 
  
end

 stat_fcov_db2_ZPX2 = [mean(fcov_db2_ZPX2);std(fcov_db2_ZPX2)];
 mean_cov_zscore_db2_ZPX2 = nanmean(cov_zscore_db2_ZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_db2_ZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_db2_ZPX2(i,:)','type','Pearson');
 corr_cov_zscore_db2_ZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_db2_ZPX2(i,:)','type','Spearman');
 corr_cov_zscore_db2_ZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_db2_ZPX2(i,:)','type','Kendall');
end

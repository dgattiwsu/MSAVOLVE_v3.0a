% This script is used to add on the fly the fbhZPX2 method. 

corr_cov_zscore_fbhZPX2 = zeros(end_cycle,3);
fbhZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_fbhZPX2 = zeros(end_cycle,3);
cov_zscore_fbhZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_fbhZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('fbhZPX2 cycle %d \n', n);

[fbhZPX2] = NMSA_to_fbhZPX2(MSA_select_ALL(:,:,n));

 % Here we calculate the usual statistics fbhZPX2.
 fbhZPX2_ALL(:,:,n) = fbhZPX2;
 
 [fcov_fbhZPX2(n,:),cov_zscore_fbhZPX2(n,:),...
     covar_vec_zscore_fbhZPX2(:,:,n)] = ... 
     coev_stats_2(fbhZPX2,ncov,covar_vec); 
  
end

 stat_fcov_fbhZPX2 = [mean(fcov_fbhZPX2);std(fcov_fbhZPX2)];
 mean_cov_zscore_fbhZPX2 = nanmean(cov_zscore_fbhZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fbhZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fbhZPX2(i,:)','type','Pearson');
 corr_cov_zscore_fbhZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fbhZPX2(i,:)','type','Spearman');
 corr_cov_zscore_fbhZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fbhZPX2(i,:)','type','Kendall');
end

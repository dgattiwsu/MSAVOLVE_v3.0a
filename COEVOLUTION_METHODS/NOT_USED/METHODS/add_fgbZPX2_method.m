% This script is used to add on the fly the fgbZPX2 method. 

corr_cov_zscore_fgbZPX2 = zeros(end_cycle,3);
fgbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_fgbZPX2 = zeros(end_cycle,3);
cov_zscore_fgbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_fgbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('fgbZPX2 cycle %d \n', n);

[fgbZPX2] = NMSA_to_fgbZPX2(MSA_select_ALL(:,:,n));

% Here we calculate the usual statistics fgbZPX2.
 fgbZPX2_ALL(:,:,n) = fgbZPX2;
 
 [fcov_fgbZPX2(n,:),cov_zscore_fgbZPX2(n,:),...
     covar_vec_zscore_fgbZPX2(:,:,n)] = ... 
     coev_stats_2(fgbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_fgbZPX2 = [mean(fcov_fgbZPX2);std(fcov_fgbZPX2)];
 mean_cov_zscore_fgbZPX2 = nanmean(cov_zscore_fgbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fgbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fgbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_fgbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fgbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_fgbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fgbZPX2(i,:)','type','Kendall');
end

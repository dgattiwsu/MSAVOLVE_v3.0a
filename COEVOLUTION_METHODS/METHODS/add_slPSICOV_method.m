% This script is used to add on the fly the slPSICOV method. 

corr_cov_zscore_slPSICOV = zeros(end_cycle,3);
slPSICOV_ALL = zeros(npos,npos,end_cycle);
fcov_slPSICOV = zeros(end_cycle,3);
cov_zscore_slPSICOV = zeros(end_cycle,ncov);
covar_vec_zscore_slPSICOV = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('slPSICOV cycle %d \n', n);

[slPSICOV] = NMSA_to_slPSICOV(MSA_select_ALL(:,:,n),'NOGAPS',0.9,1.0,...
    20,'SUM','NONE','fro',0.0001,0.35,0.015,100,0.0,0.0001,'BOTH',0,1);

 % Here we calculate the usual statistics for slPSICOV.
 slPSICOV_ALL(:,:,n) = slPSICOV;
 
 [fcov_slPSICOV(n,:),cov_zscore_slPSICOV(n,:),...
     covar_vec_zscore_slPSICOV(:,:,n)] = ... 
     coev_stats_2(slPSICOV,ncov,covar_vec); 
  
end

 stat_fcov_slPSICOV = [mean(fcov_slPSICOV);std(fcov_slPSICOV)];
 mean_cov_zscore_slPSICOV = nanmean(cov_zscore_slPSICOV,2);
 
for i = 1:end_cycle
 corr_cov_zscore_slPSICOV(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_slPSICOV(i,:)','type','Pearson');
 corr_cov_zscore_slPSICOV(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_slPSICOV(i,:)','type','Spearman');
 corr_cov_zscore_slPSICOV(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_slPSICOV(i,:)','type','Kendall');
end

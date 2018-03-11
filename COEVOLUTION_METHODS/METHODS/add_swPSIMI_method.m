% This script is used to add on the fly the swPSIMI method. Please, note
% that the default PSIMI calculation does not include any longer a
% correction for gaps.

corr_cov_zscore_PSIMI = zeros(end_cycle,3);
PSIMI_ALL = zeros(npos,npos,end_cycle);
fcov_PSIMI = zeros(end_cycle,3);
cov_zscore_PSIMI = zeros(end_cycle,ncov);
covar_vec_zscore_PSIMI = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('PSIMI cycle %d \n', n);

[PSIMI] = NMSA_to_swPSIMI(MSA_select_ALL(:,:,n),'NOGAPS',1.0,0.5,21,0.025,0.0,0.0001,'RHO');

 % Here we calculate the usual statistics for PSIMI.
 PSIMI_ALL(:,:,n) = PSIMI;
 
 [fcov_PSIMI(n,:),cov_zscore_PSIMI(n,:),...
     covar_vec_zscore_PSIMI(:,:,n)] = ... 
     coev_stats_2(PSIMI,ncov,covar_vec); 
  
end

 stat_fcov_PSIMI = [mean(fcov_PSIMI);std(fcov_PSIMI)];
 mean_cov_zscore_PSIMI = nanmean(cov_zscore_PSIMI,2);
 
for i = 1:end_cycle
 corr_cov_zscore_PSIMI(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_PSIMI(i,:)','type','Pearson');
 corr_cov_zscore_PSIMI(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_PSIMI(i,:)','type','Spearman');
 corr_cov_zscore_PSIMI(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_PSIMI(i,:)','type','Kendall');
end

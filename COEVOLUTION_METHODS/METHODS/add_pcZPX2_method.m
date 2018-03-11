% This script is used to add on the fly the pcZPX2 method. 
% MIP must be calculated first.

corr_cov_zscore_pcZPX2 = zeros(end_cycle,3);
pcZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_pcZPX2 = zeros(end_cycle,3);
cov_zscore_pcZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_pcZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('pcZPX2 cycle %d \n', n);

[pcZPX2] = NMSA_to_pcZPX2(MSA_select_ALL(:,:,n),0.08);

 % Here we calculate some statistics.
 pcZPX2_ALL(:,:,n) = pcZPX2;
 
 [fcov_pcZPX2(n,:),cov_zscore_pcZPX2(n,:),...
     covar_vec_zscore_pcZPX2(:,:,n)] = ... 
     coev_stats_2(pcZPX2,ncov,covar_vec); 
  
end

 stat_fcov_pcZPX2 = [mean(fcov_pcZPX2);std(fcov_pcZPX2)];
 mean_cov_zscore_pcZPX2 = nanmean(cov_zscore_pcZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_pcZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_pcZPX2(i,:)','type','Pearson');
 corr_cov_zscore_pcZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_pcZPX2(i,:)','type','Spearman');
 corr_cov_zscore_pcZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_pcZPX2(i,:)','type','Kendall');
end

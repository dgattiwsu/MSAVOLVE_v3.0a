% This script is used to add on the fly the nbZPX2 method. Please, note
% that the default nbZPX2 calculation does not include any longer a
% correction for gaps.

corr_cov_zscore_nbZPX2 = zeros(end_cycle,3);
nbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_nbZPX2 = zeros(end_cycle,3);
cov_zscore_nbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_nbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('nbZPX2 cycle %d \n', n);

[nbZPX2] = NMSA_to_nbZPX2(MSA_select_ALL(:,:,n));

 % Here we calculate the usual statistics nbZPX2.
 nbZPX2_ALL(:,:,n) = nbZPX2;
 
 [fcov_nbZPX2(n,:),cov_zscore_nbZPX2(n,:),...
     covar_vec_zscore_nbZPX2(:,:,n)] = ... 
     coev_stats_2(nbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_nbZPX2 = [mean(fcov_nbZPX2);std(fcov_nbZPX2)];
 mean_cov_zscore_nbZPX2 = nanmean(cov_zscore_nbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_nbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_nbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_nbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_nbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_nbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_nbZPX2(i,:)','type','Kendall');
end

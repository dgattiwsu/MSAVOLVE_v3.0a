% This script is used to add on the fly the dgbZPX2 method. 

corr_cov_zscore_dgbZPX2 = zeros(end_cycle,3);
dgbZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_dgbZPX2 = zeros(end_cycle,3);
cov_zscore_dgbZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_dgbZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('dgbZPX2 cycle %d \n', n);

[dgbZPX2] = NMSA_to_dgbZPX2(MSA_select_ALL(:,:,n));

% Here we calculate the usual statistics dgbZPX2.
 dgbZPX2_ALL(:,:,n) = dgbZPX2;
 
 [fcov_dgbZPX2(n,:),cov_zscore_dgbZPX2(n,:),...
     covar_vec_zscore_dgbZPX2(:,:,n)] = ... 
     coev_stats_2(dgbZPX2,ncov,covar_vec); 
  
end

 stat_fcov_dgbZPX2 = [mean(fcov_dgbZPX2);std(fcov_dgbZPX2)];
 mean_cov_zscore_dgbZPX2 = nanmean(cov_zscore_dgbZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_dgbZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_dgbZPX2(i,:)','type','Pearson');
 corr_cov_zscore_dgbZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_dgbZPX2(i,:)','type','Spearman');
 corr_cov_zscore_dgbZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_dgbZPX2(i,:)','type','Kendall');
end

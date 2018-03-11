% This script is used to add on the fly the nb2_ZPX2 method. Please, note
% that the default nb2_ZPX2 calculation does not include any longer a
% correction for gaps.

corr_cov_zscore_nb2_ZPX2 = zeros(end_cycle,3);
nb2_ZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_nb2_ZPX2 = zeros(end_cycle,3);
cov_zscore_nb2_ZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_nb2_ZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('nb2_ZPX2 cycle %d \n', n);

[nb2_ZPX2] = NMSA_to_nb2_ZPX2(MSA_select_ALL(:,:,n),'GAPS','DISTANCE','INCLUDE');

% Here we calculate the usual statistics nb2_ZPX2.
 nb2_ZPX2_ALL(:,:,n) = nb2_ZPX2;
 
 [fcov_nb2_ZPX2(n,:),cov_zscore_nb2_ZPX2(n,:),...
     covar_vec_zscore_nb2_ZPX2(:,:,n)] = ... 
     coev_stats_2(nb2_ZPX2,ncov,covar_vec); 
  
end

 stat_fcov_nb2_ZPX2 = [mean(fcov_nb2_ZPX2);std(fcov_nb2_ZPX2)];
 mean_cov_zscore_nb2_ZPX2 = nanmean(cov_zscore_nb2_ZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_nb2_ZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_nb2_ZPX2(i,:)','type','Pearson');
 corr_cov_zscore_nb2_ZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_nb2_ZPX2(i,:)','type','Spearman');
 corr_cov_zscore_nb2_ZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_nb2_ZPX2(i,:)','type','Kendall');
end

% This script can be used as a template to add on the fly a new coevolution
% method to the analysis. In this example we are adding ZRES calculated
% according to the ZPX method, which we call here the ZRES2 method. If your
% method is 'NEW', replace 'ZRES2' with 'NEW'.
corr_cov_zscore_ZRES2 = zeros(end_cycle,1);
corr2_cov_zscore_ZRES2 = zeros(end_cycle,1);
corr3_cov_zscore_ZRES2 = zeros(end_cycle,1);

ZRES2_ALL = zeros(npos,npos,ncycles);
fcov_ZRES2 = zeros(ncycles,3);
cbcc_COV_ZRES2 = zeros(ncycles,3);
cbcc_recomb_COV_ZRES2 = zeros(ncycles,3);
cbcc_mut_COV_ZRES2 = zeros(ncycles,3);
cbcc_cov_COV_ZRES2 = zeros(ncycles,3);
cbcc_glob_COV_ZRES2 = zeros(ncycles,3);
cov_zscore_ZRES2 = zeros(ncycles,ncov);
covar_vec_zscore_ZRES2 = zeros(ncov,3,ncycles);

for n = 1:end_cycle
 MI = NMSA_to_MI(MSA_select_ALL(:,:,n));
 [~,ZRES2] = MI_to_ZPX(MI);
 ZRES2_ALL(:,:,n) = ZRES2;
 
 [fcov_ZRES2(n,:),cov_zscore_ZRES2(n,:),...
     covar_vec_zscore_ZRES2(:,:,n)] = ... 
     coev_stats_2(ZRES2,ncov,covar_vec); 

 cbcc_COV_ZRES2(n,:) = mean_cbc_corr_all(COV_ALL(:,:,n), ZRES2);
 cbcc_recomb_COV_ZRES2(n,:) = mean_cbc_corr_all(recomb_COV_ALL(:,:,n), ZRES2);
 cbcc_mut_COV_ZRES2(n,:) = mean_cbc_corr_all(mut_COV_ALL(:,:,n), ZRES2);
 cbcc_cov_COV_ZRES2(n,:) = mean_cbc_corr_all(cov_COV_ALL(:,:,n), ZRES2);
 cbcc_glob_COV_ZRES2(n,:) = mean_cbc_corr_all(glob_COV_ALL(:,:,n), ZRES2);
end

 stat_fcov_ZRES2 = [mean(fcov_ZRES2);std(fcov_ZRES2)];
 mean_cov_zscore_ZRES2 = mean(cov_zscore_ZRES2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_ZRES2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES2(i,:)','type','Pearson');
 corr2_cov_zscore_ZRES2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES2(i,:)','type','Spearman');
 corr3_cov_zscore_ZRES2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ZRES2(i,:)','type','Kendall');
end

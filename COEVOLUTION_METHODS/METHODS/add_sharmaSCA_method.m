% This script is used to add on the fly the SSEM_SCA and RSEM_SCA methods as
% described in Rohit SHARMA PhD Thesis. The RSEM method is a little more
% accurate than the SSEM and provides results that are indistinguishable
% from the original SCA method, which is only available from Rama
% Ranganathan as a Matlab Toolbox.

SSEM_SCA_ALL = zeros(npos,npos,end_cycle);
RSEM_SCA_ALL = zeros(npos,npos,end_cycle);
fcov_SSEM_SCA = zeros(end_cycle,3);
fcov_RSEM_SCA = zeros(end_cycle,3);
cov_zscore_SSEM_SCA = zeros(end_cycle,ncov);
cov_zscore_RSEM_SCA = zeros(end_cycle,ncov);
covar_vec_zscore_SSEM_SCA = zeros(ncov,3,end_cycle);
covar_vec_zscore_RSEM_SCA = zeros(ncov,3,end_cycle);
corr_cov_zscore_SSEM_SCA = zeros(end_cycle,3);
corr_cov_zscore_RSEM_SCA = zeros(end_cycle,3);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('SCA cycle %d \n', n);

 c_MSA_select = int2aa(MSA_select_ALL(:,:,n));
 [SSEM_SCA_DDG_mat,~] = SMSA_to_SSEM_SCA(c_MSA_select);
 [~,~,SSEM_SCA] = get_sca(SSEM_SCA_DDG_mat);
  
 [RSEM_SCA_DDG_mat,~] = SMSA_to_RSEM_SCA(c_MSA_select,1000,0.5);
 [~,~,RSEM_SCA] = get_sca(RSEM_SCA_DDG_mat);

 SSEM_SCA_ALL(:,:,n) = SSEM_SCA;
 RSEM_SCA_ALL(:,:,n) = RSEM_SCA;
 
 [fcov_SSEM_SCA(n,:),cov_zscore_SSEM_SCA(n,:),...
     covar_vec_zscore_SSEM_SCA(:,:,n)] = ... 
     coev_stats_2(SSEM_SCA,ncov,covar_vec);
 [fcov_RSEM_SCA(n,:),cov_zscore_RSEM_SCA(n,:),...
     covar_vec_zscore_RSEM_SCA(:,:,n)] = ... 
     coev_stats_2(RSEM_SCA,ncov,covar_vec);
 
 clear SSEM_SCA_DDG_mat RSEM_SCA_DDG_mat
 
end

 stat_fcov_SSEM_SCA = [mean(fcov_SSEM_SCA);std(fcov_SSEM_SCA)];
 mean_cov_zscore_SSEM_SCA = nanmean(cov_zscore_SSEM_SCA,2);
 stat_fcov_RSEM_SCA = [mean(fcov_RSEM_SCA);std(fcov_RSEM_SCA)];
 mean_cov_zscore_RSEM_SCA = nanmean(cov_zscore_RSEM_SCA,2);

for i = 1:end_cycle
 corr_cov_zscore_SSEM_SCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_SSEM_SCA(i,:)','type','Pearson');
 corr_cov_zscore_SSEM_SCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_SSEM_SCA(i,:)','type','Spearman');
 corr_cov_zscore_SSEM_SCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_SSEM_SCA(i,:)','type','Kendall');
 corr_cov_zscore_RSEM_SCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_RSEM_SCA(i,:)','type','Pearson');
 corr_cov_zscore_RSEM_SCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_RSEM_SCA(i,:)','type','Spearman');
 corr_cov_zscore_RSEM_SCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_RSEM_SCA(i,:)','type','Kendall');
end
 
% This script can be used to add on the fly hpPCA. 

% First we remove all external files.
!rm msa.faln dca_mat.matrix

corr_cov_zscore_hpPCA = zeros(end_cycle,3);
hpPCA_ALL = zeros(npos,npos,end_cycle);
fcov_hpPCA = zeros(end_cycle,3);
cov_zscore_hpPCA = zeros(end_cycle,ncov);
covar_vec_zscore_hpPCA = zeros(ncov,3,end_cycle);


for n = 1:end_cycle
fprintf('hpPCA cycle %d \n', n);

[hpPCA] = get_nmsa_covar_vec(MSA_select_ALL(:,:,n),30,'hpPCA');

for i = 1:npos
    hpPCA(i,i) = NaN;
end
 
 % Here we remove the external files generated at each cycle.
 !rm msa.faln 
 
% hpPCA already has a MIP correction, so we only add a ZPX2 correction.
 hpPCA_ZPX2 = MIP_to_ZPX2(hpPCA);

% Gap correction
hpPCA_ZPX2_1 = hpPCA_ZPX2;
hpPCA_ZPX2_2 = hpPCA_ZPX2_1 - min(hpPCA_ZPX2_1(:));
hpPCA_ZPX2 = hpPCA_ZPX2_2.*gW(:,:,n);
  
% Here we calculate the usual statistics for hpPCA.
 hpPCA_ALL(:,:,n) = hpPCA_ZPX2;
 
 % Here we calculate the usual statistics for hpPCA. 
 [fcov_hpPCA(n,:),cov_zscore_hpPCA(n,:),...
     covar_vec_zscore_hpPCA(:,:,n)] = ... 
     coev_stats_2(hpPCA_ZPX2,ncov,covar_vec); 

 % Clean up before the next cycle.
 clear Lambda Vtilde N q p lambda_diag nlambda hpPCA
  
end

 stat_fcov_hpPCA = [nanmean(fcov_hpPCA);std(fcov_hpPCA)];
 mean_cov_zscore_hpPCA = mean(cov_zscore_hpPCA,2);

for i = 1:end_cycle
 corr_cov_zscore_hpPCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_hpPCA(i,:)','type','Pearson');
 corr_cov_zscore_hpPCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_hpPCA(i,:)','type','Spearman');
 corr_cov_zscore_hpPCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_hpPCA(i,:)','type','Kendall');
end
 
% This script can be used to add on the fly gplmDCA. 

% First we remove all external files.
!rm msa.faln dca_mat.matrix

% corr_cov_zscore_gplmDCA = zeros(end_cycle,3);
% gplmDCA_ALL = zeros(npos,npos,end_cycle);
% fcov_gplmDCA = zeros(end_cycle,3);
% cov_zscore_gplmDCA = zeros(end_cycle,ncov);
% covar_vec_zscore_gplmDCA = zeros(ncov,3,end_cycle);

% for n = 1:end_cycle
for n = 59:93
fprintf('gplmDCA cycle %d \n', n);

 [gplmDCA] = get_nmsa_covar_vec(MSA_select_ALL(:,:,n),30,'gplmDCA_asym');

 % Here we convert to ZPX2
 gplmDCA_ZPX2 = MIP_to_ZPX2(gplmDCA);
 
 % Here we calculate a matrix of weights to correct for the presence of
 % gaps.

gapW0 = ones(npos,npos);
gapW1 = correct_coevmat_forgaps(MSA_select_ALL(:,:,n));
gapW2 = gapW1.^2;
gapW3 = gapW1.^3;

% Here we choose how to corrects the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW1,gapW2,gapW3) can be applied.
gW = gapW3;

gplmDCA_ZPX2_1 = gplmDCA_ZPX2;
gplmDCA_ZPX2_2 = gplmDCA_ZPX2_1 - min(gplmDCA_ZPX2_1(:));
gplmDCA_ZPX2 = gplmDCA_ZPX2_2.*gW;

 % Here we remove the external files generated at each cycle.
 clear dca_mat
 !rm msa.faln dca_mat.matrix
 
 % Here we calculate the usual statistics for gplmDCA.
 gplmDCA_ALL(:,:,n) = gplmDCA_ZPX2;
 
 [fcov_gplmDCA(n,:),cov_zscore_gplmDCA(n,:),...
     covar_vec_zscore_gplmDCA(:,:,n)] = ... 
     coev_stats_2(gplmDCA_ZPX2,ncov,covar_vec); 

end

 stat_fcov_gplmDCA = [nanmean(fcov_gplmDCA);std(fcov_gplmDCA)];
 mean_cov_zscore_gplmDCA = mean(cov_zscore_gplmDCA,2);

for i = 1:end_cycle
 corr_cov_zscore_gplmDCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_gplmDCA(i,:)','type','Pearson');
 corr_cov_zscore_gplmDCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_gplmDCA(i,:)','type','Spearman');
 corr_cov_zscore_gplmDCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_gplmDCA(i,:)','type','Kendall');
end
 

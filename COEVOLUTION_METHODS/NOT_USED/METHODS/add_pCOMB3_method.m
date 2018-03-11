% This script is used to add on the fly the pCOMB_min method. If your
% method is 'NEW', replace 'pCOMB_min' with 'NEW'. This script is based on the 
% method developed in 'positional_MI.m'.

corr_cov_zscore_pZPX2_min = zeros(end_cycle,3);
pZPX2_min_ALL = zeros(npos,npos,ncycles);
% pZPX2_min_bestsum_ALL = zeros(end_cycle,2);
fcov_pZPX2_min = zeros(ncycles,3);
cov_zscore_pZPX2_min = zeros(ncycles,ncov);
covar_vec_zscore_pZPX2_min = zeros(ncov,3,ncycles);

corr_cov_zscore_pCOMB_min = zeros(end_cycle,3);
pCOMB_min_ALL = zeros(npos,npos,ncycles);
% pCOMB_min_bestsum_ALL = zeros(end_cycle,2);
fcov_pCOMB_min = zeros(ncycles,3);
cov_zscore_pCOMB_min = zeros(ncycles,ncov);
covar_vec_zscore_pCOMB_min = zeros(ncov,3,ncycles);

% cbcc_COV_pCOMB_min = zeros(ncycles,3);
% cbcc_recomb_COV_pCOMB_min = zeros(ncycles,3);
% cbcc_mut_COV_pCOMB_min = zeros(ncycles,3);
% cbcc_cov_COV_pCOMB_min = zeros(ncycles,3);
% cbcc_glob_COV_pCOMB_min = zeros(ncycles,3);

% sig_pCORR = mean(mean_cov_zscore_pCORR);
% sig_pZPX2 = mean(mean_cov_zscore_pZPX2);
% scale_pCORR_pZPX2 = sig_pZPX2/sig_pCORR;

[nrows,ncols] = size(MSA_select);
% ref_row = 1;

for n = 1:end_cycle
fprintf('pCOMB_min cycle %d \n', n);

% [~,pCOMB_min,~,~,~,~,~,~] = ...
%     NMSA_to_pCOMB_min(MSA_select_ALL(:,:,n),0.7);
[pZPX2_min,~,pCOMB_min,~,~,~] = ...
    NMSA_to_pCOMB2(MSA_select_ALL(:,:,n),0.7,1,1);

% figure;imagesc(pCOMB_min);set(gca,'Ydir','normal');
% sorted_pCOMB_min = sort_matrix_descend(pCOMB_min);
% pCOMB_min_covar_vec = sorted_pCOMB_min(1:ncov,2:3); 
% pCOMB_min_hits = numel(intersect(covar_vec,pCOMB_min_covar_vec,'rows'))/2

 % Here we calculate the usual statistics pCOMB_min.
 pCOMB_min_ALL(:,:,n) = pCOMB_min;
 pZPX2_min_ALL(:,:,n) = pZPX2_min;
%  pCOMB_min_bestsum_ALL(n,1) = mean(pCOMB_min_min_sum(:,1));
%  pCOMB_min_bestsum_ALL(n,2) = mean(pCOMB_min_min_sum(:,2));
 
 [fcov_pCOMB_min(n,:),cov_zscore_pCOMB_min(n,:),...
     covar_vec_zscore_pCOMB_min(:,:,n)] = ... 
     coev_stats_2(pCOMB_min,ncov,covar_vec); 

 [fcov_pZPX2_min(n,:),cov_zscore_pZPX2_min(n,:),...
     covar_vec_zscore_pZPX2_min(:,:,n)] = ... 
     coev_stats_2(pZPX2_min,ncov,covar_vec); 
 
%  cbcc_COV_pCOMB_min(n,:) = mean_cbc_corr_all(COV_ALL(:,:,n), pCOMB_min);
%  cbcc_recomb_COV_pCOMB_min(n,:) = mean_cbc_corr_all(recomb_COV_ALL(:,:,n), pCOMB_min);
%  cbcc_mut_COV_pCOMB_min(n,:) = mean_cbc_corr_all(mut_COV_ALL(:,:,n), pCOMB_min);
%  cbcc_cov_COV_pCOMB_min(n,:) = mean_cbc_corr_all(cov_COV_ALL(:,:,n), pCOMB_min);
%  cbcc_glob_COV_pCOMB_min(n,:) = mean_cbc_corr_all(glob_COV_ALL(:,:,n), pCOMB_min);

end

 stat_fcov_pCOMB_min = [mean(fcov_pCOMB_min);std(fcov_pCOMB_min)];
 mean_cov_zscore_pCOMB_min = mean(cov_zscore_pCOMB_min,2);

 stat_fcov_pZPX2_min = [mean(fcov_pZPX2_min);std(fcov_pZPX2_min)];
 mean_cov_zscore_pZPX2_min = mean(cov_zscore_pZPX2_min,2);
 
 
for i = 1:end_cycle
 corr_cov_zscore_pCOMB_min(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_pCOMB_min(i,:)','type','Pearson');
 corr_cov_zscore_pCOMB_min(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_pCOMB_min(i,:)','type','Spearman');
 corr_cov_zscore_pCOMB_min(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_pCOMB_min(i,:)','type','Kendall');
 corr_cov_zscore_pZPX2_min(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_pZPX2_min(i,:)','type','Pearson');
 corr_cov_zscore_pZPX2_min(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_pZPX2_min(i,:)','type','Spearman');
 corr_cov_zscore_pZPX2_min(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_pZPX2_min(i,:)','type','Kendall');
end

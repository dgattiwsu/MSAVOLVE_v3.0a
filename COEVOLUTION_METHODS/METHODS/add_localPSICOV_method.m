% This script is used to add on the fly the localPSICOV method. Please, note
% that the default localPSICOV calculation does not include any longer a
% correction for gaps.

corr_cov_zscore_localPSICOV = zeros(end_cycle,3);
localPSICOV_ALL = zeros(npos,npos,end_cycle);
fcov_localPSICOV = zeros(end_cycle,3);
cov_zscore_localPSICOV = zeros(end_cycle,ncov);
covar_vec_zscore_localPSICOV = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('localPSICOV cycle %d \n', n);

[~,~,localPSICOV] = ...
NMSA_to_PSICOV(MSA_select_ALL(:,:,n),'NOGAPS',1.0,0.5,'ZERO','FRO',...
                    0.005,0.0,0.0001,'QUIC','RHO',20);

 % Here we calculate the usual statistics for localPSICOV.
 localPSICOV_ALL(:,:,n) = localPSICOV;
 
 [fcov_localPSICOV(n,:),cov_zscore_localPSICOV(n,:),...
     covar_vec_zscore_localPSICOV(:,:,n)] = ... 
     coev_stats_2(localPSICOV,ncov,covar_vec); 
  
end

 stat_fcov_localPSICOV = [mean(fcov_localPSICOV);std(fcov_localPSICOV)];
 mean_cov_zscore_localPSICOV = nanmean(cov_zscore_localPSICOV,2);
 
for i = 1:end_cycle
 corr_cov_zscore_localPSICOV(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_localPSICOV(i,:)','type','Pearson');
 corr_cov_zscore_localPSICOV(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_localPSICOV(i,:)','type','Spearman');
 corr_cov_zscore_localPSICOV(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_localPSICOV(i,:)','type','Kendall');
end

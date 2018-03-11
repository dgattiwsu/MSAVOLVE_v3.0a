% This script is used to add on the fly the fpcZPX2 method. 

corr_cov_zscore_fpcZPX2 = zeros(end_cycle,3);
fpcZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_fpcZPX2 = zeros(end_cycle,3);
cov_zscore_fpcZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_fpcZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('fpcZPX2 cycle %d \n', n);

[fpcZPX2] = NMSA_to_fpcZPX2(MSA_select_ALL(:,:,n),'NOGAPS',1.0,'SDP',1,'AUTO','QUIC',0.005,...
      0.0,0.001,'INVERSE',21,'NONE',0.05);
  
 % Here we calculate some statistics.
 fpcZPX2_ALL(:,:,n) = fpcZPX2;
 
 [fcov_fpcZPX2(n,:),cov_zscore_fpcZPX2(n,:),...
     covar_vec_zscore_fpcZPX2(:,:,n)] = ... 
     coev_stats_2(fpcZPX2,ncov,covar_vec); 
  
end

 stat_fcov_fpcZPX2 = [mean(fcov_fpcZPX2);std(fcov_fpcZPX2)];
 mean_cov_zscore_fpcZPX2 = nanmean(cov_zscore_fpcZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_fpcZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fpcZPX2(i,:)','type','Pearson');
 corr_cov_zscore_fpcZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fpcZPX2(i,:)','type','Spearman');
 corr_cov_zscore_fpcZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fpcZPX2(i,:)','type','Kendall');
end

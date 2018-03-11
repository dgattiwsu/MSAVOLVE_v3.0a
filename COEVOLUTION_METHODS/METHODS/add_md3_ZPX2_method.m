% This script is used to add on the fly the md3_ZPX2 method. 

corr_cov_zscore_md3_ZPX2 = zeros(end_cycle,3);
md3_ZPX2_ALL = zeros(npos,npos,end_cycle);
fcov_md3_ZPX2 = zeros(end_cycle,3);
cov_zscore_md3_ZPX2 = zeros(end_cycle,ncov);
covar_vec_zscore_md3_ZPX2 = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL);

% The following are the default values for mdMI parameters ----------------
if exist('sim_cutoff','var')
else
    sim_cutoff = 0.6;
end
if exist('res_cutoff','var')
else
    res_cutoff = 22;
end
if exist('gapcorr_1','var')
else
    gapcorr_1 = 1.0;
end
if exist('gapcorr_2','var')
else
    gapcorr_2 = 1.0;
end
if exist('gapcorr_3','var')
else
    gapcorr_3 = 3.0;
end
if exist('parproc','var')
else
    par_proc = 12;
end
%--------------------------------------------------------------------------

for n = 1:end_cycle
fprintf('md3_ZPX2 cycle %d \n', n);

[~,md3_ZPX2] = ...
    NMSA_to_mdMI(MSA_select_ALL(:,:,n),'GAPS','3D','FULL',sim_cutoff,1,...
    res_cutoff,gapcorr_1,gapcorr_2,gapcorr_3,par_proc);

 % Here we calculate some statistics.
 md3_ZPX2_ALL(:,:,n) = md3_ZPX2;
 
 [fcov_md3_ZPX2(n,:),cov_zscore_md3_ZPX2(n,:),...
     covar_vec_zscore_md3_ZPX2(:,:,n)] = ... 
     coev_stats_2(md3_ZPX2,ncov,covar_vec); 
  
end

 stat_fcov_md3_ZPX2 = [mean(fcov_md3_ZPX2);std(fcov_md3_ZPX2)];
 mean_cov_zscore_md3_ZPX2 = nanmean(cov_zscore_md3_ZPX2,2);
 
for i = 1:end_cycle
 corr_cov_zscore_md3_ZPX2(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_md3_ZPX2(i,:)','type','Pearson');
 corr_cov_zscore_md3_ZPX2(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_md3_ZPX2(i,:)','type','Spearman');
 corr_cov_zscore_md3_ZPX2(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_md3_ZPX2(i,:)','type','Kendall');
end

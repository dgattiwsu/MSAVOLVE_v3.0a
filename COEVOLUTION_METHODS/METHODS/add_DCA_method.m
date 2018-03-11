% This script can be used to add on the fly wMI and DCA. 

% First we remove all external files.
!rm msa.faln dca_mat.matrix

corr_cov_zscore_DCA = zeros(end_cycle,3);
DCA_ALL = zeros(npos,npos,end_cycle);
fcov_DCA = zeros(end_cycle,3);
cov_zscore_DCA = zeros(end_cycle,ncov);
covar_vec_zscore_DCA = zeros(ncov,3,end_cycle);

corr_cov_zscore_wMI = zeros(end_cycle,3);
wMI_ALL = zeros(npos,npos,end_cycle);
fcov_wMI = zeros(end_cycle,3);
cov_zscore_wMI = zeros(end_cycle,ncov);
covar_vec_zscore_wMI = zeros(ncov,3,end_cycle);

for n = 1:end_cycle
fprintf('DCA cycle %d \n', n);

 % Here we write out the msa in fasta format.   
 nmsa_to_faln(MSA_select_ALL(:,:,n),'msa.faln');
 
 % Here we perform DCA using the original code from Andrea Pagnani.
 dca('msa.faln','dca_mat.matrix');
 
 % Here we read in the dca matrix
 import_dca_matrix('dca_mat.matrix');
 
 % Here we convert the dca_mat.matrix into our internal matrices for both DCA
 % and MI weighted in the same way as DCA
 [ wMI,DCA,~,~ ] = dca_to_wMI_DCA( dca_mat,fcov );
 
 % Here we remove the external files generated at each cycle.
 clear dca_mat
 !rm msa.faln dca_mat.matrix
 
 % Here we calculate the usual statistics for both DCA and wMI.
 DCA_ALL(:,:,n) = DCA;
 
 [fcov_DCA(n,:),cov_zscore_DCA(n,:),...
     covar_vec_zscore_DCA(:,:,n)] = ... 
     coev_stats_2(DCA,ncov,covar_vec); 

 wMI_ALL(:,:,n) = wMI;
 
 [fcov_wMI(n,:),cov_zscore_wMI(n,:),...
     covar_vec_zscore_wMI(:,:,n)] = ... 
     coev_stats_2(wMI,ncov,covar_vec); 
  
end

 stat_fcov_DCA = [nanmean(fcov_DCA);std(fcov_DCA)];
 mean_cov_zscore_DCA = mean(cov_zscore_DCA,2);
 stat_fcov_wMI = [nanmean(fcov_wMI);std(fcov_wMI)];
 mean_cov_zscore_wMI = mean(cov_zscore_wMI,2);

for i = 1:end_cycle
 corr_cov_zscore_DCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_DCA(i,:)','type','Pearson');
 corr_cov_zscore_DCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_DCA(i,:)','type','Spearman');
 corr_cov_zscore_DCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_DCA(i,:)','type','Kendall');
 corr_cov_zscore_wMI(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_wMI(i,:)','type','Pearson');
 corr_cov_zscore_wMI(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_wMI(i,:)','type','Spearman');
 corr_cov_zscore_wMI(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_wMI(i,:)','type','Kendall');
end
 
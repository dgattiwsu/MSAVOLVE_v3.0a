% This script can be used to add on the fly plmDCA. 

% First we remove all external files.
!rm msa.faln dca_mat.matrix

corr_cov_zscore_plmDCA = zeros(end_cycle,3);
plmDCA_ALL = zeros(npos,npos,end_cycle);
fcov_plmDCA = zeros(end_cycle,3);
cov_zscore_plmDCA = zeros(end_cycle,ncov);
covar_vec_zscore_plmDCA = zeros(ncov,3,end_cycle);

for n = 1:end_cycle
fprintf('DCA cycle %d \n', n);

 % Here we write out the msa in fasta format.   
 nmsa_to_faln(MSA_select_ALL(:,:,n),'msa.faln');
 
 % Here we perform plmDCA using the original code.
 plmDCA_symmetric('msa.faln','dca_mat.matrix',0.01,0.01,0.1,12)        

 % Here we read in the dca matrix
 import_dca_matrix('dca_mat.matrix');
 
 % Here we convert the dca_mat.matrix into our internal matrix.
 [ plmDCA,~] = plm_dca_to_wMI_DCA( dca_mat,fcov );

 % Here we remove the external files generated at each cycle.
 clear dca_mat
 !rm msa.faln dca_mat.matrix
 
 % Here we calculate the usual statistics for both DCA and wMI.
 plmDCA_ALL(:,:,n) = plmDCA;
 
 [fcov_plmDCA(n,:),cov_zscore_plmDCA(n,:),...
     covar_vec_zscore_plmDCA(:,:,n)] = ... 
     coev_stats_2(plmDCA,ncov,covar_vec); 

end

 stat_fcov_plmDCA = [nanmean(fcov_plmDCA);std(fcov_plmDCA)];
 mean_cov_zscore_plmDCA = mean(cov_zscore_plmDCA,2);

for i = 1:end_cycle
 corr_cov_zscore_plmDCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_plmDCA(i,:)','type','Pearson');
 corr_cov_zscore_plmDCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_plmDCA(i,:)','type','Spearman');
 corr_cov_zscore_plmDCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_plmDCA(i,:)','type','Kendall');
end
 
% This script adds on the fly the logR coevolution method (logR). All the
% required C and perl programs must be present in the directory from which
% this script is launched.

% First we remove all external files.
!rm msa.aln msa.faln msa.fasta msa.post postmi_mat.matrix
!rm msa.mapping msa_ungapped.aln msa_ungapped.post

corr_cov_zscore_logR = zeros(end_cycle,3);
logR_ALL = zeros(npos,npos,end_cycle);
fcov_logR = zeros(end_cycle,3);
cov_zscore_logR = zeros(end_cycle,ncov);
covar_vec_zscore_logR = zeros(ncov,3,end_cycle);

for n = 1:end_cycle
fprintf('logR cycle %d \n', n);

 % Here we write out the msa in fasta format.
 nmsa_to_faln(MSA_select_ALL(:,:,n),'msa.faln');
 
 % Here we perform logR using the original code from Lukas Burger:
 % lukas.burger@fmi.ch
 !t_coffee -other_pg seq_reformat -in msa.faln -output fasta_aln > msa.fasta
 !perl runContactPredictions.pl msa.fasta
 !mv msa.post postmi_mat.matrix
 
 % Here we read in the logR matrix
 import_postmi_matrix('postmi_mat.matrix');
 
 % Here we convert the postmi_mat.matrix into our internal matrices for both logR
 % and MI weighted in the same way as logR
 [ logR,~] = postmi_to_logR( postmi_mat,npos,fcov );
 ones_ind = logR == 1;
 logR(ones_ind) = NaN;
 
 % Here we remove the internal and external files generated at each cycle.
 clear postmi_mat
 !rm msa.aln msa.faln msa.fasta postmi_mat.matrix
 !rm msa.mapping msa_ungapped.aln msa_ungapped.post
 
 % Here we calculate the usual statistics.
 logR_ALL(:,:,n) = logR;
 
 [fcov_logR(n,:),cov_zscore_logR(n,:),...
     covar_vec_zscore_logR(:,:,n)] = ... 
     coev_stats_2(logR,ncov,covar_vec); 
 
end

 stat_fcov_logR = [mean(fcov_logR);std(fcov_logR)];
 mean_cov_zscore_logR = nanmean(cov_zscore_logR,2);
 
for i = 1:end_cycle
 corr_cov_zscore_logR(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_logR(i,:)','type','Pearson');
 corr_cov_zscore_logR(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_logR(i,:)','type','Spearman');
 corr_cov_zscore_logR(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_logR(i,:)','type','Kendall');
end


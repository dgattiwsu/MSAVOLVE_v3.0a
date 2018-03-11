% This script calculates correlation coefficients between coevolution 
% matrices by first converting them into a linear vector that does not
% include the diagonal elements.

COV = COV_ALL(:,:,n);
cov_COV = cov_COV_ALL(:,:,n);
mut_COV = mut_COV_ALL(:,:,n);
recomb_COV = recomb_COV_ALL(:,:,n);
MI = MI_ALL(:,:,n);
NMI = NMI_ALL(:,:,n);
ZMI = ZMI_ALL(:,:,n);
MIP = MIP_ALL(:,:,n);
ZPX = ZPX_ALL(:,:,n);
ZPX2 = ZPX2_ALL(:,:,n);
ZRES2 = ZRES2_ALL(:,:,n);
ELSC = ELSC_ALL(:,:,n);
ELSC = ELSC_ALL(:,:,n);
McBASC = McBASC_ALL(:,:,n);
omesSCA = omesSCA_ALL(:,:,n);
RAMA_SCA = RAMA_SCA_ALL(:,:,n);
RSEM_SCA = RSEM_SCA_ALL(:,:,n);
SSEM_SCA = SSEM_SCA_ALL(:,:,n);

ut_ind = logical(triu(ones(280),1));

lin_COV = COV(ut_ind);
lin_covCOV = cov_COV(ut_ind);
lin_mutCOV = mut_COV(ut_ind);
lin_recombCOV = recomb_COV(ut_ind);
lin_MI = MI(ut_ind);
lin_NMI = NMI(ut_ind);
lin_MIP = MIP(ut_ind);
lin_ZMI = ZMI(ut_ind);
lin_ZPX = ZPX(ut_ind);
lin_ZPX2 = ZPX2(ut_ind);
lin_ZRES2 = ZRES2(ut_ind);
lin_OMES = OMES(ut_ind);
lin_ELSC = ELSC(ut_ind);
lin_McBASC = McBASC(ut_ind);
lin_omesSCA = omesSCA(ut_ind);
lin_RAMA_SCA = RAMA_SCA(ut_ind);
lin_SSEM_SCA = SSEM_SCA(ut_ind);
lin_RSEM_SCA = RSEM_SCA(ut_ind);

coev_corr_table = zeros(11,6);
coev_corr_table(1,1) = corr(linCOV,lin_MI,'type','Pearson');
coev_corr_table(1,2) = corr(linCOV,lin_MI,'type','Spearman');
coev_corr_table(1,3) = corr(linCOV,lin_MI,'type','Kendall');
coev_corr_table(1,4) = corr(lin_covCOV,lin_MI,'type','Pearson');
coev_corr_table(1,5) = corr(lin_covCOV,lin_MI,'type','Spearman');
coev_corr_table(1,6) = corr(lin_covCOV,lin_MI,'type','Kendall');
coev_corr_table(2,1) = corr(linCOV,lin_MIP,'type','Pearson');
coev_corr_table(2,2) = corr(linCOV,lin_MIP,'type','Spearman');
coev_corr_table(2,3) = corr(linCOV,lin_MIP,'type','Kendall');
coev_corr_table(2,4) = corr(lin_covCOV,lin_MIP,'type','Pearson');
coev_corr_table(2,5) = corr(lin_covCOV,lin_MIP,'type','Spearman');
coev_corr_table(2,6) = corr(lin_covCOV,lin_MIP,'type','Kendall');
coev_corr_table(3,1) = corr(linCOV,lin_ZPX,'type','Pearson');
coev_corr_table(3,2) = corr(linCOV,lin_ZPX,'type','Spearman');
coev_corr_table(3,3) = corr(linCOV,lin_ZPX,'type','Kendall');
coev_corr_table(3,4) = corr(lin_covCOV,lin_ZPX,'type','Pearson');
coev_corr_table(3,5) = corr(lin_covCOV,lin_ZPX,'type','Spearman');
coev_corr_table(3,6) = corr(lin_covCOV,lin_ZPX,'type','Kendall');
coev_corr_table(4,1) = corr(linCOV,lin_ZPX2,'type','Pearson');
coev_corr_table(4,2) = corr(linCOV,lin_ZPX2,'type','Spearman');
coev_corr_table(4,3) = corr(linCOV,lin_ZPX2,'type','Kendall');
coev_corr_table(4,4) = corr(lin_covCOV,lin_ZPX2,'type','Pearson');
coev_corr_table(4,5) = corr(lin_covCOV,lin_ZPX2,'type','Spearman');
coev_corr_table(4,6) = corr(lin_covCOV,lin_ZPX2,'type','Kendall');
coev_corr_table(5,1) = corr(linCOV,lin_ZRES2,'type','Pearson');
coev_corr_table(5,2) = corr(linCOV,lin_ZRES2,'type','Spearman');
coev_corr_table(5,3) = corr(linCOV,lin_ZRES2,'type','Kendall');
coev_corr_table(5,4) = corr(lin_covCOV,lin_ZRES2,'type','Pearson');
coev_corr_table(5,5) = corr(lin_covCOV,lin_ZRES2,'type','Spearman');
coev_corr_table(5,6) = corr(lin_covCOV,lin_ZRES2,'type','Kendall');
coev_corr_table(6,1) = corr(linCOV,lin_ZMI,'type','Pearson');
coev_corr_table(6,2) = corr(linCOV,lin_ZMI,'type','Spearman');
coev_corr_table(6,3) = corr(linCOV,lin_ZMI,'type','Kendall');
coev_corr_table(6,4) = corr(lin_covCOV,lin_ZMI,'type','Pearson');
coev_corr_table(6,5) = corr(lin_covCOV,lin_ZMI,'type','Spearman');
coev_corr_table(6,6) = corr(lin_covCOV,lin_ZMI,'type','Kendall');
coev_corr_table(7,1) = corr(linCOV,lin_OMES,'type','Pearson');
coev_corr_table(7,2) = corr(linCOV,lin_OMES,'type','Spearman');
coev_corr_table(7,3) = corr(linCOV,lin_OMES,'type','Kendall');
coev_corr_table(7,4) = corr(lin_covCOV,lin_OMES,'type','Pearson');
coev_corr_table(7,5) = corr(lin_covCOV,lin_OMES,'type','Spearman');
coev_corr_table(7,6) = corr(lin_covCOV,lin_OMES,'type','Kendall');
coev_corr_table(8,1) = corr(linCOV,lin_ELSC,'type','Pearson');
coev_corr_table(8,2) = corr(linCOV,lin_ELSC,'type','Spearman');
coev_corr_table(8,3) = corr(linCOV,lin_ELSC,'type','Kendall');
coev_corr_table(8,4) = corr(lin_covCOV,lin_ELSC,'type','Pearson');
coev_corr_table(8,5) = corr(lin_covCOV,lin_ELSC,'type','Spearman');
coev_corr_table(8,6) = corr(lin_covCOV,lin_ELSC,'type','Kendall');
coev_corr_table(9,1) = corr(linCOV,lin_McBASC,'type','Pearson');
coev_corr_table(9,2) = corr(linCOV,lin_McBASC,'type','Spearman');
coev_corr_table(9,3) = corr(linCOV,lin_McBASC,'type','Kendall');
coev_corr_table(9,4) = corr(lin_covCOV,lin_McBASC,'type','Pearson');
coev_corr_table(9,5) = corr(lin_covCOV,lin_McBASC,'type','Spearman');
coev_corr_table(9,6) = corr(lin_covCOV,lin_McBASC,'type','Kendall');
coev_corr_table(10,1) = corr(linCOV,lin_omesSCA,'type','Pearson');
coev_corr_table(10,2) = corr(linCOV,lin_omesSCA,'type','Spearman');
coev_corr_table(10,3) = corr(linCOV,lin_omesSCA,'type','Kendall');
coev_corr_table(10,4) = corr(lin_covCOV,lin_omesSCA,'type','Pearson');
coev_corr_table(10,5) = corr(lin_covCOV,lin_omesSCA,'type','Spearman');
coev_corr_table(10,6) = corr(lin_covCOV,lin_omesSCA,'type','Kendall');
coev_corr_table(11,1) = corr(linCOV,lin_RAMA_SCA,'type','Pearson');
coev_corr_table(11,2) = corr(linCOV,lin_RAMA_SCA,'type','Spearman');
coev_corr_table(11,3) = corr(linCOV,lin_RAMA_SCA,'type','Kendall');
coev_corr_table(11,4) = corr(lin_covCOV,lin_RAMA_SCA,'type','Pearson');
coev_corr_table(11,5) = corr(lin_covCOV,lin_RAMA_SCA,'type','Spearman');
coev_corr_table(11,6) = corr(lin_covCOV,lin_RAMA_SCA,'type','Kendall');









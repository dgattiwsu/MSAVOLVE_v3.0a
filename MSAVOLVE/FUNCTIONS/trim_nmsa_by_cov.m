%-- 09/23/2011 03:53:23 PM --%
clear trim_ind out_ind
REF_cov_mat_upper = triu(REF_cov_mat,1);
[trim_ind(:,1),trim_ind(:,2)] = find(REF_cov_mat_upper' > 0.0378);
keep = setdiff([1:nseq],unique(trim_ind(:,1)));
REF_nmsa_t1 = REF_nmsa(keep,:);
REF_msa_t1 = REF_msa(keep);
REF_nmsa_bk = REF_nmsa;
REF_msa_bk = REF_msa;
REF_nmsa = REF_nmsa_t1;
REF_msa = REF_msa_t1;

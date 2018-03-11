function [DDG_mat,DG_mat]=SMSA_to_RSEM_SCA(cmsa,ntrial,fkeep)
% This function calculates the DDG perturbation matrix using the Random 
% Sequence Elimination Method. 'fkeep' is the fraction (e.g. 0.5) of
% sequences that are retained at each trial
%
alphabet=('ACDEFGHIKLMNPQRSTVWY');
scale=100; % used to normalize to 100 sequences
rr_factx = rr_fact(scale);
bg_prob = [0.072658 0.024692 0.050007 0.061087 0.041774 0.071589 0.023392...
0.052691 0.063923 0.089093 0.023150 0.042931 0.052228 0.039871...
0.052012 0.073087 0.055606 0.063321 0.012720 0.032955];
[nseq,npos] = size(cmsa);
naa = length(alphabet);
profile = zeros(naa, npos);
sub_profile = zeros(naa, npos);
DG_mat = zeros(naa,npos);
sub_DG_mat = zeros(naa,npos);
DDG_mat = zeros(naa,npos,ntrial);

% determine profile (normalized to 100 seqs) in the msa and calculate DG_mat.
    for aa = 1:naa
        profile(aa,1:npos) = sum(cmsa == alphabet(aa)).*scale/nseq;
        DG_mat(aa,:) = rr_factx - gammaln(profile(aa,:)+1)...
        - gammaln(100-profile(aa,:)+1) + profile(aa,:)*log(bg_prob(aa))...
        + (scale-profile(aa,:))*log(1-bg_prob(aa));
    end;

    for k = 1:ntrial
    
        n = randi(nseq,round(nseq*fkeep),1);
        sub_nseq = length(n);
        sub_cmsa = cmsa(n,:);
    
% determine profile (normalized to 100 seqs) in sub_msa and calculate 
% sub_DG_mat.
        for aa = 1:naa
        sub_profile(aa,1:npos) = sum(sub_cmsa == alphabet(aa)).*scale/sub_nseq;
        sub_DG_mat(aa,:) = rr_factx - gammaln(sub_profile(aa,:)+1)...
        - gammaln(100-sub_profile(aa,:)+1) + sub_profile(aa,:)*log(bg_prob(aa))...
        + (scale-sub_profile(aa,:))*log(1-bg_prob(aa));
        end;        

        DDG_mat(:,:,k) = DG_mat-sub_DG_mat;

    end
    
end


function [ msa_covar_vec_low,msa_covar_vec_mid,msa_covar_vec_high ] = ...
    pick_covarions_5( msa,msa_entropy,npos,fr_cov,cov_gaps )

% Here we pick the covarying columns. For example, we want 10%(fr_cov) of 
% all the 300 positions of the sequence to be covarying with another 10%, 
% which means a total of 60 columns. We will pick two vectors of 30 columns. 
% In this function we pick three sets of covarying columns with, 
% respectively, the largest or lowest entropies, and entropies around the  
% center of the normal distribution of entropies in the reference msa.

if cov_gaps == 0
% First we find the indices of all the positions for which there are no 
% gaps in the MSA.
msa_nogapsposmat = msa ~= 25;
msa_nogapspos = logical(prod(+msa_nogapsposmat));
msa_nogaps_entropy = msa_entropy;
msa_nogaps_entropy(~msa_nogapspos) = 0;
msa_entropy = msa_nogaps_entropy;
end

ncov = round(npos*fr_cov/100);
msa_covar_vec_high = zeros(ncov,2);
% msa_entropy = Entropy(msa)';
msa_entropy = msa_entropy';
[~,o_msa_entropy_ind_high] = sort(msa_entropy,'descend');
msa_covar_vec_high(:,1) = o_msa_entropy_ind_high(1:2:(ncov*2));
msa_covar_vec_high(:,2) = o_msa_entropy_ind_high(2:2:(ncov*2+1));
entropy_threshold = max(msa_entropy)/100;

% We take the mid-values from the center of the normal distribution of
% entropies above a minimum threshold.
msa_nonzeroentr_ind = msa_entropy > entropy_threshold;
msa_nonzeroentr_pd = fitdist(msa_entropy(msa_nonzeroentr_ind),'normal');
msa_entr_diff = abs(msa_entropy - msa_nonzeroentr_pd.Params(1));
[~,o_msa_entropy_ind_mid_1] = sort(msa_entr_diff,'ascend');
o_msa_entropy_ind_mid_all = o_msa_entropy_ind_mid_1(1:ncov*2);
[~,o_msa_entropy_ind_mid_2] = sort(msa_entropy,'ascend');
o_msa_entropy_ind_mid_3 = ismember(o_msa_entropy_ind_mid_2,...
    o_msa_entropy_ind_mid_all);
o_msa_entropy_ind_mid = o_msa_entropy_ind_mid_2(o_msa_entropy_ind_mid_3);
msa_covar_vec_mid(:,1) = o_msa_entropy_ind_mid(1:2:(ncov*2));
msa_covar_vec_mid(:,2) = o_msa_entropy_ind_mid(2:2:(ncov*2+1));

msa_entr_diff = msa_entropy - entropy_threshold;
[~,o_msa_entropy_ind_low_1] = sort(msa_entr_diff,'ascend');
o_msa_entropy_ind_low_2 = find(msa_entropy > entropy_threshold);
o_msa_entropy_ind_low_3 = ismember(o_msa_entropy_ind_low_1,o_msa_entropy_ind_low_2);
o_msa_entropy_ind_low = o_msa_entropy_ind_low_1(o_msa_entropy_ind_low_3);
msa_covar_vec_low(:,1) = o_msa_entropy_ind_low(1:2:(ncov*2));
msa_covar_vec_low(:,2) = o_msa_entropy_ind_low(2:2:(ncov*2+1));

end


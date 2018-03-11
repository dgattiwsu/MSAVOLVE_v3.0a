function [ nmsa_rel_entropy ] = get_nmsa_rel_entropy( nmsa )
% This function calculates a vector representing the relative entropy of
% each position of a nmsa.

nmsa_hmm_model = nmsa_to_hmm_model(nmsa);
bg_frequencies = nmsa_hmm_model.NullEmission';
% Next we calculate the profile of the MSA without including gaps.
cmsa = int2aa(nmsa);
[nmsa_profile,~] = seqprofile(cmsa,'Gaps','none');
% Here we determine the level of conservation using the traditional 
% definitions of entropy and relative entropy.
nmsa_entropy = Entropy(nmsa);
[~,nmsa_rel_entropy] = rel_entropy_nmsa(nmsa,...
    bg_frequencies,nmsa_profile);

end


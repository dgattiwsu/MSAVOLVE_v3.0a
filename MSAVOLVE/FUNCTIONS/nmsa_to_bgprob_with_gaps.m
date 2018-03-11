function [ msa_bgprob ] = nmsa_to_bgprob_with_gaps( nmsa )
% This function calculates the background probabilities for each aa in a
% msa in matlab numeric format from the msa profile.

[~,npos] = size(nmsa);
msa = int2aa(nmsa);
msa_profile = seqprofile(msa,'gaps','all');
msa_bgprob = sum(msa_profile')/npos;
end


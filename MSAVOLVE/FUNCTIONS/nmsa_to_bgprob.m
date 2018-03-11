function [ msa_bgprob ] = nmsa_to_bgprob( nmsa )
% This function calculates the background probabilities for each aa in a
% msa in matlab numeric format from the msa profile.

[~,npos] = size(nmsa);
msa = int2aa(nmsa);
msa_profile = seqprofile(msa);
msa_bgprob = sum(msa_profile')/npos;
end


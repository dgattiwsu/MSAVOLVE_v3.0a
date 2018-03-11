function [ hmm_model ] = nmsa_to_hmm_model( nmsa )
% This function derives an HMM model from an MSA in Matlab numeric format.

[~,npos] = size(nmsa);
msa_bgprob = nmsa_to_bgprob( nmsa );
hmm_structure = hmmprofstruct(npos,'Alphabet','AA','NullEmission',msa_bgprob);
hmm_model = hmmprofestimate(hmm_structure,int2aa(nmsa));
end

function [ MSA ] = make_ancestor( nseq,npos,PD )
% Here we generate an ancestor protein of npos aa's and we repeat the
% sequence nseq times. Background sequences are imported from a hmm model.
% bg_frequencies = ALL_KDO8PS_hmm_model.NullEmission';

% symbols = [1:1:20]';
% PD = fitdist(symbols,'kernel','width',0.01,'frequency',frequencies);
SEQ = zeros(1,npos);
for i=1:npos
SEQ(i) = round(random(PD));
end
MSA = repmat(SEQ,nseq,1);

end


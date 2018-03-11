function [ MSA ] = make_ancestor_2( nseq,npos,prob_PD )
% Here we generate an ancestor protein of npos aa's and we repeat the
% sequence nseq times. 

SEQ = zeros(1,npos);
for i=1:npos
SEQ(i) = round(random(prob_PD{i}));
end
MSA = repmat(SEQ,nseq,1);

end


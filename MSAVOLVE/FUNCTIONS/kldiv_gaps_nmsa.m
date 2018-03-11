function [ rH_ai,rH_i,D_ai,D_i,Pm_ai ] = ...
    kldiv_gaps_nmsa( nmsa,bg_prob,seq_profile )
% This function returns a traditional expression for relative entropy (rH_i) 
% and the corresponding expression (D_i) used in SCA. 
% PM[f(a)i] is the probability of observing f(a)i in an alignment of M 
% sequences given a background probability q(a).
% The value of D(a)i indicates how unlikely the observed frequency f(a)i of 
% amino acid a at position i would be if a occurred randomly with 
% probability q(a).
[nseq,npos] = size(nmsa);
rel_entropy_1 = zeros(21,npos);
rel_entropy_2 = zeros(21,npos);
rel_entropy_3 = zeros(21,npos);
Pm_ai = zeros(21,npos);
D_ai = zeros(21,npos);
rH_ai = zeros(21,npos);

for i = 1:npos
    rel_entropy_1(:,i) = ...
        nantozero(...
        seq_profile(:,i) .* ...
        log(seq_profile(:,i)./bg_prob)/log(21)...
        );
    rel_entropy_2(:,i) = ...
        nantozero(...
        (1-seq_profile(:,i)) .* ...
        log( (1-seq_profile(:,i))./(1-bg_prob) )/log(21)...
        );
    rel_entropy_3(:,i) = rel_entropy_1(:,i) + rel_entropy_2(:,i);
end

rH_ai = rel_entropy_1;
rH_i = sum(rel_entropy_1);
D_ai = rel_entropy_3;
D_i = sum(rel_entropy_3);
Pm_ai(:) = exp(-npos*rel_entropy_3);


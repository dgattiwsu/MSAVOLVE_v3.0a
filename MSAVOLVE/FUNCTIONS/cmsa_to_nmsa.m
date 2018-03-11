function nmsa = cmsa_to_nmsa(cmsa)
% Converts an msa from single letter to matlab numeric format. Identical to
% matlab 'aa2int' function.

alphabet='ARNDCQEGHILKMFPSTWYVBX_*-';

[nseq,npos] = size(cmsa);
nmsa=zeros(nseq,npos);

for i = 1:25
    symbol = alphabet(i);
    ind = cmsa == symbol;
    nmsa(ind) = i;
end


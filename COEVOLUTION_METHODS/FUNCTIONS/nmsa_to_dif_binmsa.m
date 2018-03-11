function [ dif_binmsa,cum_dif_binmsa,dif_binmsa_2,cum_dif_binmsa_2,...
    dif_binmsa_3,cum_dif_binmsa_3] = ...
    nmsa_to_dif_binmsa( nmsa )
% This function converts an msa from matlab numeric format (nmsa) to various 
% types of binary formats. If only one output is used the differential 
% binary format used in dbZPX2 is produced.

[nrows,ncols] = size(nmsa);
dif_binmsa = false(nrows,ncols);
    for i = 2:nrows
        dif_binmsa(i,:) = nmsa(i,:)~=nmsa(i-1,:);
    end
cum_dif_binmsa = cumsum(dif_binmsa);

dif_binmsa_2 = false(nrows,ncols);
    for i = 2:nrows
        dif_binmsa_2(i,:) = nmsa(i,:)~=nmsa(1,:);
    end
cum_dif_binmsa_2 = cumsum(dif_binmsa_2);

dif_binmsa_3 = false(nrows,ncols);
    for i = 2:nrows
        dif_binmsa_3(i,:) = nmsa(i,:)==nmsa(i-1,:);
    end
cum_dif_binmsa_3 = cumsum(dif_binmsa_3);

end


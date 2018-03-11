function [ nmsa ] = ...
    make_ancestor_3B( nmsa,branch_start,branch_end,npos,prob_PD,freq_PD,ncycles )

[nseq,~] = size(nmsa);
nmsa_1 = nmsa;
nmsa_s = zeros(nseq,npos);
nmsa1_s = zeros(nseq,npos);

for cycle = 1:ncycles
    for i = branch_start:branch_end
        for j=1:npos
        nmsa_1(i,j) = round(random(prob_PD{j}));
        nmsa_n = nmsa_1(i,j);        
        nmsa1_s(i,j) = freq_PD(nmsa_n,j);    
        end
    end
evolved = nmsa1_s >= nmsa_s;
nmsa(evolved) = nmsa_1(evolved);
nmsa_s = nmsa1_s;
end

end



function [ordered] = reorder_nmsa_1(nmsa,nrows,ncols)
%--------------------------------------------------------------------------    
    for k = 1:nrows-1
    nmsa_dif = zeros(nrows,1);
        for i = k:nrows
        nmsa_dif(i) = sum(nmsa(k,:)~=nmsa(i,:));
        end
    [nmsa_dif, sort_msa_ind] = sort(nmsa_dif,'ascend');
    nmsa = nmsa(sort_msa_ind,:);
    end
    ordered = nmsa;
end


function [ gapW ] = correct_coevmat_forgaps_3D( nmsa )
% This function finds the positions in the msa where there are no gaps in both
% columns considered in a MI calculation and calculates a matrix of weights
% to scale the coevolution matrix.

[nseq,npos] = size(nmsa);
gapW = zeros(npos,npos,npos);

for i = 1:npos
    for j = i:npos
        for k = j:npos
        gap = find((nmsa(:,i)~=25) & (nmsa(:,j)~=25) & (nmsa(:,k)~=25));
        gapW(i,j,k) = numel(gap);
        gapW(j,i,k) = gapW(i,j,k);
        gapW(i,k,j) = gapW(i,j,k);
        gapW(k,i,j) = gapW(i,j,k);
        gapW(j,k,i) = gapW(i,j,k);
        gapW(k,j,i) = gapW(i,j,k);
        end
    end       
end

gapW=gapW/nseq;

end


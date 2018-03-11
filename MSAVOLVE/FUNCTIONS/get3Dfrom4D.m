function [COV3D] = get3Dfrom4D(COV4D,cycle,seq_ind)
% This function returns a 2D COV matrix from a 4D COV matrix. 'cycle' is
% the cycle we want to extract. 'seq_ind' is the vector of indices of the
% sequences we want to extract from a specific cycle. For example:
% [COV10_23to45 = get2Dfrom4D(COV3D_ALL,15,[23:45] will extract the COV
% matrix of sequences 23 to 45 from cycle 15.

COV3D = COV4D(:,:,seq_ind,cycle);
end

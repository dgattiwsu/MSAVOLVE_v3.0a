function [ smsa ] = nmsa_to_smsa( nmsa )
% This function converts a nmsa (matlab numeric format) into a structure of
% sequences that can be written out to an ascii file as a fasta or faln
% format.

[rows,~] = size(nmsa);
smsa(1,rows) = struct();
for i = 1:rows
    smsa(1,i).Header = ['Sequence_' int2str(i)];
    smsa(1,i).Sequence = int2aa(nmsa(i,:));
end
end


function [MI_test,MI_test_sum] = NMSA_to_MIsum(nmsa)
% This fucntion calculates the MI matrix of a nmsa in the normal way and
% after converting the nmsa to a different numeric format in which every aa
% in a certain column is replaced by the count of how many times that aa 
% appears in that column. Useful only to check that indeed the two MI
% calculations give almost identical results.

[nrows,ncols] = size(nmsa);
% nmsa_pos_sum = zeros(nrows,ncols);
nmsa_sum = zeros(nrows,ncols);

% rows = 1:nrows;
% for j = 1:nrows
% ref_row = j;
% test_rows = setdiff(rows,ref_row);
% nmsa_pos = false(nrows,ncols);
%     for i = test_rows
%       nmsa_pos(i,:) = nmsa(i,:)~=nmsa(ref_row,:);
%     end
% nmsa_pos_sum = nmsa_pos + nmsa_pos_sum;
% end

for i = 1:20
    nmsa_i = nmsa == i;
    nmsa_i_sum = sum(nmsa_i);
    for j = 1:nrows
        nmsa_ij = nmsa_i(j,:);
        nmsa_sum(j,nmsa_ij) = nmsa_i_sum(nmsa_ij);
    end
end

MI_test = NMSA_to_MI(nmsa);
% MI_test_pos = NMSA_to_MI(nmsa_pos_sum);
MI_test_sum = NMSA_to_MI(nmsa_sum);
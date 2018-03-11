function [ZPX2] = NMSA_to_nbZPX2(nmsa)

% nbZPX2 flowchart
% 1.	Calculate a distance matrix for the sequences of the msa.
% 2.	Reorder the msa by placing as the 1st and 2nd row the two sequences 
%       most similar to each other.
% 3.	Place in the 3rd row the sequence most similar to the sequence in 
%       the 2nd row. 
% 4.	Loop by placing in each consecutive row the sequence most similar 
%       to the sequence in the previous row, until the entire msa is reordered.
% 5.	Convert the msa to binary differential: the 1st row is assigned all 
%       Os. In each consecutive row a 0 is placed at every position at which 
%       the aa is the same as in the previous row, and a 1 at every position 
%       in which the aa is different.
% 6.	Go back to the reordered msa in standard format and change to 0 
%       every position that is a 0 in the binary differential msa: the 
%       result is the normal binary msa or nb_msa.
% 7.	Calculate the MI matrix for the nb_msa.
% 8.	Calculate a ZPX2 matrix from the MI matrix.
% 

% 1.	Calculate a distance matrix for the sequences of the msa. Here more 
%       distant sequences are represented by smaller numbers.

% Replace unusual symbols if any are present
ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

[nrows,ncols] = size(nmsa);

% The following loop is the traditional way of calculating the distance
% matrix:
% dist = zeros(nrows,nrows);
% for i = 1:nrows
%     for j = i:nrows
%         dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
%         dist(j,i) = dist(i,j);
%     end
% end

% However, the following method is much faster than the previous loop. It's
% based on the trick of first converting the nmsa to a long binary format
% as described by Halabi et al. Cell 138: 774-786, 2009.

bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';

% 2.	Here we reorder the msa by placing as the 1st and 2nd row the two  
%       sequences most similar to each other. In practice we only 
%       need to find out the index of the sequence that will be at the top.

triu_dist = triu(dist,1);
[~,max_ind] = max(triu_dist(:));
[firstrow,~] = ind2sub([nrows,nrows],max_ind);

% 3.	The index of the 1st sequence is the one we just found.

        ordered_ind = zeros(nrows,1);
        ordered_ind(1) = firstrow;        
        ind = firstrow;

% 4.	Then we loop by finding for each consecutive row the index of the  
%       sequence most similar to the sequence in the previous row, until 
%       all the indices of the msa are reordered.

        for i = 2:nrows
            % zero the entire test column
            dist(:,ind) = 0;
            % find the index of the closest sequence
            [~,newind] = max(dist(ind,:));
            ordered_ind(i) = newind;
            % the next index will be the index of the sequence just selected
            ind = newind;
        end        

%       Finally we reorder the msa based on the reordered indices.

        ordered = nmsa(ordered_ind,:);
            
% 5.	Here we convert the msa to binary differential: the 1st row is assigned 
%       all Os. In each consecutive row a 0 is placed at every position at 
%       which the aa is the same as in the previous row, and a 1 at every 
%       position in which the aa is different.

    bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
    
    binsum = sum(bin_ordered(:));
    fprintf('Binary Sum = %f \n', binsum);    
    
% 6.	Here we go back to the reordered msa in standard format and change  
%       to 0 every position that is a 0 in the binary differential msa: the 
%       result is the normal binary msa or nb_msa.

    z_ordered = ordered;
    z_ordered(~bin_ordered) = 0;

% 7.	Calculate the MI matrix for the nb_msa.

    [MI] = NMSA_to_MI(z_ordered);
    
% 8.	Calculate a ZPX2 matrix from the MI matrix.

    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
end



function [binmsa] = nmsa_to_binmsa_21q(nmsa)
% Returns each sequence of length L as a vector of size 21L with 0 and 1. 
% Number 21 represents gaps (which would be # 25 in the original in the 
% original Matlab numeric representation of an MSA).

[nseq,npos]=size(nmsa);
ind25 = nmsa == 25;
nmsa(ind25) = 21;
binmsa=zeros(nseq,21*npos);
for i=1:npos 
    for aa=1:21 
        binmsa(:,21*(i-1)+aa)=(nmsa(:,i)==aa); 
    end; 
end;
end


function [pMIP] = MI_to_MIP(pMI,ncols)
%--------------------------------------------------------------------------
% MIP calculation
mean_mat = nanmean(pMI(:));
mean_row = zeros(ncols,1);
MCA_mat = zeros(ncols,ncols);

% Here  we calculate the MCA matrix.

for m = 1:ncols
    mean_row(m) = nanmean(pMI(m,:));
end

for m = 1:ncols
    for n = m:ncols
    MCA_mat(m,n)=(mean_row(m)*mean_row(n))/mean_mat;
    MCA_mat(n,m) = MCA_mat(m,n);    
    end
MCA_mat(m,m) = NaN;
end

% Finally we subtract the MCA matrix from the MI matrix.
pMIP = pMI-MCA_mat;

end


function [pZPX2] = MIP_to_ZPX2(pMIP,ncols)
%--------------------------------------------------------------------------
% ZPX2 calculation
mean_row=zeros(ncols,1);
std_row=zeros(ncols,1);
pZPX2=zeros(ncols,ncols);

for m=1:ncols
    mean_row(m)=nanmean(pMIP(m,:));
    std_row(m)=nanstd(pMIP(m,:));   
end
for m=1:ncols
    for n=m:ncols

    ZPX2_i=(pMIP(m,n)-mean_row(m))/std_row(m);
    ZPX2_j=(pMIP(m,n)-mean_row(n))/std_row(n);
        
    pZPX2(m,n)=ZPX2_i*ZPX2_j;

% Here we correct for the product of two negative ZPX2_i and ZPX2_j, which
% would give the wrong MI. Change of sign is not in the original ZPX2 algorithm by
% Gloor, but is included in the ZRES algorithm by Chen.
    
    if (ZPX2_i<0&&ZPX2_j<0)
        pZPX2(m,n)=-pZPX2(m,n);
    end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end



function [ZPX2,MIP,MI] = ...
    NMSA_to_nb2_ZPX2(nmsa,dist_method,sort_method,gaps_method)
% NMSA_to_nb2_ZPX2 expands NMSA_nbZPX3 by offering the
% option of reordering the msa using the distance matrix as such 
% (sort_method = 'DISTANCE'), the eigen decomposition of the distance 
% matrix (sort_method = 'EIGEN'), or the sums of the columns of the distance
% matrix (sort_method = 'WEIGHTS'). With 'EIGEN' or 'WEIGTHS' the final result 
% may be better for the detection of covarying pairs separated by at
% least 5 positions in sequence. The distance matrix is calculated by first
% converting the MSA to its binary 'long' representation in which every
% column is extended into 20 (dist_method = 'NOGAPS) or 21 (dist_method =
% 'GAPS') different columns, each one representing a different aa. The flag
% 'gaps_method' can be set to 'INCLUDE' or 'EXCLUDE' to include or exclude
% gaps in the MI calculation. It is difficult to recommend a best syntax for
% every sequence: simulations can be very helpful for this purpose. Examples  
% of syntax are: 
% [nb2_ZPX2] = NMSA_to_nb2_ZPX2(nmsa,'GAPS','DISTANCE','INCLUDE');
% [nb2_ZPX2] = NMSA_to_nb2_ZPX2(nmsa,'GAPS','DISTANCE','EXCLUDE'); 
% [nb2_ZPX2] = NMSA_to_nb2_ZPX2(nmsa,'GAPS','WEIGHTS','EXCLUDE'); 
% [nb2_ZPX2] = NMSA_to_nb2_ZPX2(nmsa,'NOGAPS','EIGEN','EXCLUDE'); 
% A result identical to the simpler NMSA_to_nbZPX2 is obtained with: 
% [nb2_ZPX2] = NMSA_to_nb2_ZPX2(nmsa,'GAPS','DISTANCE','INCLUDE');

% Set the default run (which is not necessarily the best option):
if nargin == 1
    dist_method = 'NOGAPS';
    sort_method = 'DISTANCE';
    gaps_method = 'EXCLUDE';
end
    
% Replace unusual symbols
ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

[nrows,ncols] = size(nmsa);

% Here we calculate the distance matrix.

switch dist_method
    case 'NOGAPS'
        if nargin == 2
            sort_method = 'DISTANCE';
            gaps_method = 'EXCLUDE';
        end
        
        bin_ordered = nmsa_to_binmsa_20q(nmsa);        
        dist = bin_ordered*bin_ordered';

        sdist = zeros(nrows,nrows);      
        for i = 1:nrows
            sdist(i,:) = dist(i,:)/dist(i,i);
        end

        udist = triu(sdist,1);
        ldist = tril(sdist,-1)';
        lind = udist < ldist;
        udist(lind) = 0;
        ldist(~lind) = 0;
        mdist = udist + ldist;
        mdist = mdist + mdist' + eye(nrows);

    % Here we find which sequence will be at the top in the reordered msa.
        triu_dist = triu(mdist,1);
        [~,max_ind] = max(triu_dist(:));
        [firstrow,firstcol] = ind2sub([nrows,nrows],max_ind);

    % Choose as 1st the sequence with the smallest number of gaps, if the
    % two sequences are of different length.
    
        if dist(firstcol,firstcol) > dist(firstrow,firstrow)
            firstrow = firstcol;
        end
        
    % Now that we have selected the 1st sequence, the merged distance
    % matrix becomes the final distance matrix.
    
        dist = mdist;
    
    case 'GAPS'
        if nargin == 2
            sort_method = 'DISTANCE';            
            gaps_method = 'INCLUDE';
        end
        
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = bin_ordered*bin_ordered';

    % Here we find which sequence will be at the top in the reordered msa.
        triu_dist = triu(dist,1);
        [~,max_ind] = max(triu_dist(:));
        [firstrow,~] = ind2sub([nrows,nrows],max_ind);

end

%--------------------------------------------------------------------------
% Here we reorder the indices of the msa and of the weights vector. We have
% three options: we can use the distance matrix as such, its eigen
% decomposition, or the columns sums.

switch sort_method
%--------------------------------------------------------------------------
    case 'DISTANCE'

        % The index of the 1st sequence is the one we found before.
        ordered_ind = zeros(nrows,1);
        ordered_ind(1) = firstrow;
        
        ind = firstrow;

        for i = 2:nrows
            % zero the entire test column
            dist(:,ind) = 0;
            % find the index of the closest sequence
            [~,newind] = max(dist(ind,:));
            ordered_ind(i) = newind;
            % the next index will be the index of the sequence just selected
            ind = newind;
        end        
%--------------------------------------------------------------------------
    case 'EIGEN'
    % Here we reorder based on the eigenvector corresponding to the  
    % largest absolute eigenvalue of the distance matrix.
    % [V,~] = eigs(dist,1,'lm');
    [V,D]=eig(dist);
    [sorted_D,ind]=sort(abs(diag(D)),'descend');
    top_V = V(:,ind(1:nrows));
    top_D = sorted_D(1:nrows);

    top_D_mat = zeros(nrows);

    for i = 1:nrows
    top_D_mat(i,i) = top_D(i,1);
    end

    for i=1:nrows
    switch_sign = sign(mean(top_V(:,i)));
    top_V(:,i) = switch_sign*top_V(:,i);
    end
    % V = sum(V,2);
    % switch_sign = sign(mean(V));
    % V = switch_sign*V;
    [~,ordered_ind] = sort(abs(top_V(:,1)),'descend');    
%--------------------------------------------------------------------------
    case 'WEIGHTS'
    % Here we reorder based on the sum of the columns of the distance
    % matrix.        
    dist_sum = sum(dist);
    [~,ordered_ind] = sort(dist_sum,'descend');
    ordered_ind = ordered_ind';
%--------------------------------------------------------------------------    
end

% Here we reorder the msa based on the reordered indices.

    ordered = nmsa(ordered_ind,:);
    
% Here we convert the msa from matlab numeric format to binary
% differential.

bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
    
binsum = sum(bin_ordered(:));
fprintf('Binary Sum = %f \n', binsum);

%Here we set to zero all the indices of the reordered nmsa that are 0 in 
% the binary msa.

z_ordered = ordered;
z_ordered(~bin_ordered) = 0;

% Here we convert to 0 all the positions that have a value of 21.
% This is equivalent to assuming the each time there is a gap in the real
% msa, the position does not change in the reordered msa.

switch gaps_method
    case 'EXCLUDE'
        ind_21 = z_ordered == 21;
        z_ordered(ind_21) = 0;
        post_binsum = numel(find(z_ordered));
        fprintf('Gap Corrected Binary Sum = %f \n', post_binsum);
    case 'INCLUDE'
end

    % MI matrix
    [MI] = NMSA_to_MI(z_ordered);
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
                
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
% would give the wrong MI. Comment: I am not sure the following three lines
% make a difference. Change of sign is not in the original ZPX2 algorithm by
% Gloor, but is included in the ZRES algorithm by Chen. At any rate the
% change of sign seems to affect only i,j positions with very small counts,
% and therefore the final effect of the change is marginal.
    
    if (ZPX2_i<0&&ZPX2_j<0)
        pZPX2(m,n)=-pZPX2(m,n);
    end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end



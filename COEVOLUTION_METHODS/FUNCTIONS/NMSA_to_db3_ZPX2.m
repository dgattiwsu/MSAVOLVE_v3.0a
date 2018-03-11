function [ZPX2] = ...
    NMSA_to_db3_ZPX2(nmsa,nshuffle,dist_method,sort_method,gaps_method)
% This function differs from 'db2' by the fact that it subtract a db
% matrix based on the 'null' evolutionary hypothesis (no correlation
% between columns) from the db map derived from the msa. The null
% hypothesis is tested by shufling 'nshuffle' times every column of the
% msa, and calculating the db matrix for each shuffled msa. Then, all the
% shuffled db matrices are averaged. The original 'db2 function 
% expands the 'db' method by allowing a selection of the
% algorithm used to obtain the distance matrix, by including (dist_method =
% 'GAPS') or not including (dist_method = 'NOGAPS') gaps. After the distance 
% matrix is determined, there are three options to reorder the msa: we can 
% use the distance matrix as such (sort_method = 'DISTANCE'), its eigen 
% decomposition (sort_method = 'EIGEN'), or the sum of its columns 
% (sort_method = 'WEIGHTS'). Then, it is also
% possible to include (gap_method = 'INCLUDE') or exclude (gap_method =
% 'EXCLUDE') the gaps from the differential binary matrix. The default
% synthax, which will give a result essentially identical to NMSA_dbZPX2 is: 
% [dbZPX2] = NMSA_to_db2_ZPX2(nmsa,'GAPS','DISTANCE','INCLUDE');

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

[MI_ref,COV_ref] = NMSA_to_MI_COV(nmsa,dist_method,sort_method,...
            gaps_method,nrows,ncols);
        
MI_all = zeros(ncols,ncols,nshuffle);
COV_all = zeros(ncols,ncols,nshuffle);
shuffled_nmsa = zeros(nrows,ncols);

for i = 1:nshuffle
    for j = 1:ncols
        shuffle_ind = randperm(nrows);
        shuffled_nmsa(:,j) = nmsa(shuffle_ind,j);
    end
    fprintf('Shuffled MSA no. %d \n', i);    
    [MI_all(:,:,i),COV_all(:,:,i)] = NMSA_to_MI_COV(shuffled_nmsa,...
        dist_method,sort_method,...
        gaps_method,nrows,ncols);
end

    MI_mean = nanmean(MI_all,3);
    COV_mean = nanmean(COV_all,3);
           
    MI = (MI_ref - MI_mean);
    COV = (COV_ref - COV_mean);

% Here we merge MI and correlation matrix by linear regression. 

% Best scale for MI/COV

    scale = scale_matrices(COV,MI);

    % fprintf('Best Scale MI/COV = %f \n', scale);

% We combine the two matrices.

    MI = MI + scale*COV;

% Here we calculate the ZPX2 matrix.
    
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
    % fprintf('Posterior Mean Binary MSA Entropy = %f \n', pH);
            
end


function [MI,COV] = ...
            NMSA_to_MI_COV(nmsa,dist_method,sort_method,...
            gaps_method,nrows,ncols)

% Here we calculate the distance matrix 

switch dist_method
    case 'NOGAPS'
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
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;

    % Here we find which sequence will be at the top in the reordered msa.
        triu_dist = triu(dist,1);
        [~,max_ind] = max(triu_dist(:));
        [firstrow,~] = ind2sub([nrows,nrows],max_ind);
end

%--------------------------------------------------------------------------
% Here we reorder the indices of the msa. We have three options: we can use 
% the distance matrix as such, its eigen decomposition, or the sum of its 
% columns.

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

%--------------------------------------------------------------------------
% Here we reorder the msa based on the reordered indices.

    ordered = nmsa(ordered_ind,:);
    
%--------------------------------------------------------------------------

% Here we convert the msa from matlab numeric format to binary
% differential.

    bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
    
% Here we convert to 0 all the positions that have a value of 21 in the 
% standard msa and a value of 1 in the binary msa.
% This is equivalent to assuming the each time there is a gap in the real
% msa, the position does not change in the reordered msa.

switch gaps_method
    case 'EXCLUDE'
        ind_21 = ordered == 21;
        z_ordered = logical(bin_ordered .* ind_21);
        bin_ordered = + bin_ordered;
        % fprintf('Binary Sum = %f \n', sum(sum(bin_ordered)));
        bin_ordered(z_ordered) = 0;
        % fprintf('Gap Corrected Binary Sum = %f \n', sum(sum(bin_ordered)));        
    case 'INCLUDE'
end

% MI matrix and mean column entropy
    
    [MI,~] = BNMSA_to_MI(bin_ordered,nrows,ncols);
    COV = cov(bin_ordered);
    

end


function [scale] = scale_matrices(mat1,mat2)
% This function returns the best scale such that scale*mat1 = mat2.
[rows,cols] = size(mat1);
template = ones(rows,cols);
upper = triu(template,1);
triu_ind = template == upper;
mat1 = mat1(triu_ind);
mat2 = mat2(triu_ind);
scale = mat1\mat2;
% scale1 = mat1\mat2;
% scale1 = (mat1'*mat2)/(mat1'*mat1);
% scale2 = mat2\mat1;
% scale2 = (mat1'*mat2)/(mat2'*mat2);
end


function [pMI,mH] = BNMSA_to_MI(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% MI calculation; the following are the rules used in generating the
% arrays.
% pI(1) == 1 ; pI(2) == 0 
% pJ(1) == 1 ; pJ(2) == 0
% pIJ(1) == 1,1
% pIJ(2) == 1,0
% pIJ(3) == 0,1
% pIJ(4) == 0,0
bnmsa = + bnmsa;
const = log2(nrows);
H = zeros(1,ncols);
jH = zeros(ncols,ncols);
pMI = zeros(ncols,ncols);
pIJ = zeros(4,1);
% epIJ = zeros(4,1);

pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        % H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
        H(I) = const - (sum(pI(ind,I).*(log2(pI(ind,I)))))/nrows;        
end
mH = mean(H);

for I = 1:ncols
    for J = I:ncols
        pIJ(1) = bnmsa(:,I)'*bnmsa(:,J);
        pIJ(2) = pI(1,I) - pIJ(1);        
        pIJ(3) = pI(1,J) - pIJ(1);
        pIJ(4) = nrows - pIJ(1) - pIJ(2) - pIJ(3);
%         epIJ(1) = pI(1,I)*pI(1,J);
%         epIJ(2) = pI(1,I)*pI(2,J);
%         epIJ(3) = pI(2,I)*pI(1,J);
%         epIJ(4) = pI(2,I)*pI(2,J);
        ind = find(pIJ);
%         jH(I,J) = -(sum(pIJ(ind) .* (log2(pIJ(ind)) - const)))/nrows;
        jH(I,J) = const - (sum(pIJ(ind) .* (log2(pIJ(ind)))))/nrows;
        pMI(I,J) = H(I) + H(J) - jH(I,J);
%         pMI(I,J) = (sum(pIJ(ind) .* (log2(pIJ(ind) ./ epIJ(ind)) + const)))/nrows;
    end
end
pMI = pMI+pMI';
for m = 1:ncols
pMI(m,m)=NaN;
end

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



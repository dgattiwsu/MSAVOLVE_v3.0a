function [ZPX2] = NMSA_to_dgbZPX2(nmsa)
% In this function we reorder the msa such that there is the minimum 
% distance between each sequence and the following one.
% First we determine the distance matrix: here more distant sequences are
% represented by smaller numbers.

[nrows,ncols] = size(nmsa);

% dist = zeros(nrows,nrows);
% for i = 1:nrows
%     for j = i:nrows
%         dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
%         dist(j,i) = dist(i,j);
%     end
%     % We need to set the diagonal to 0 in order to find the largest
%     % non-diagonal terms.
%     dist(i,i) = 0;
% end

% Here we calculate the distance matrix. The following method is much
% faster than the previous for loop.
bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';

triu_dist = triu(dist,1);
[~,max_ind] = max(triu_dist(:));
[firstrow,~] = ind2sub([nrows,nrows],max_ind);

nmsa1 = nmsa(firstrow:nrows,:);
nmsa2 = nmsa(1:firstrow-1,:);
nmsa = [nmsa1;nmsa2];

%--------------------------------------------------------------------------
% Start of 1st reordering algorithm
%--------------------------------------------------------------------------

ordered = reorder_nmsa_1(nmsa,nrows,ncols);

%--------------------------------------------------------------------------
% Before we proceed further we first calculate the bMI matrix.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);

[bMI,pH] = BNMSA_to_MI(bin_ordered,bin_ordered_ind,nrows,ncols);

% Here we first convert the msa from matlab numeric to long format binary.

ordered = nmsa_to_binmsa_21q(ordered);

% Here we convert the long binary to binary differential.

[nrows,ncols] = size(ordered);
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);

% Here we calculate the cov matrix only for the indices included in the
% vector 'bin_ordered_ind'.

    gCOV = zeros(ncols,ncols);
    [ogCOV] = cov(bin_ordered(:,bin_ordered_ind));
    gCOV(bin_ordered_ind,bin_ordered_ind) = ogCOV;
    
    % Here we recover the original values of nrows and ncols

    [~,ncols] = size(nmsa);
    COV = zeros(ncols,ncols);
    for oi = 1:ncols
        i = 21*(oi-1) + 1;
        for oj = oi:ncols
            j = 21*(oj-1)+1;
            COV(oi,oj) = norm(gCOV(i:i+20,j:j+20),'fro');
            COV(oj,oi) = COV(oi,oj);
        end
    end

    bMI_1 = bMI;
    pH_1 = pH;
    COV_1 = COV;
%--------------------------------------------------------------------------
% End the 1st reordering algorithm
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Start of 2nd reordering algorithm
%--------------------------------------------------------------------------

ordered = reorder_nmsa_3(nmsa,nrows,ncols);

%--------------------------------------------------------------------------
% Before we proceed further we first calculate the bMI matrix.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);

[bMI,pH] = BNMSA_to_MI(bin_ordered,bin_ordered_ind,nrows,ncols);

% Here we first convert the msa from matlab numeric to long format binary.

ordered = nmsa_to_binmsa_21q(ordered);

% Here we convert the long binary to binary differential.

[nrows,ncols] = size(ordered);
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);

% Here we calculate the cov matrix only for the indices included in the
% vector 'bin_ordered_ind'.

    gCOV = zeros(ncols,ncols);
    [ogCOV] = cov(bin_ordered(:,bin_ordered_ind));
    gCOV(bin_ordered_ind,bin_ordered_ind) = ogCOV;
    
    % Here we recover the original values of nrows and ncols

    [~,ncols] = size(nmsa);
    COV = zeros(ncols,ncols);
    for oi = 1:ncols
        i = 21*(oi-1) + 1;
        for oj = oi:ncols
            j = 21*(oj-1)+1;
            COV(oi,oj) = norm(gCOV(i:i+20,j:j+20),'fro');
            COV(oj,oi) = COV(oi,oj);
        end
    end

    bMI_2 = bMI;
    pH_2 = pH;
    COV_2 = COV;
%--------------------------------------------------------------------------
% End the 1st reordering algorithm
%--------------------------------------------------------------------------

    % Here we first merge everything from the two runs.

    pH = (pH_1 + pH_2)/2;
    MI = bMI_1 + bMI_2;
    COV = COV_1 + COV_2;
    
    % Here we merge the MI and COV matrices by linear regression.
    % Best scale for COV/MI.
    
    scale = scale_matrices(MI,COV);

    fprintf('best scale COV/MI = %f \n', scale);

    % Finally we combine the two matrices.

    MI = COV + scale*MI;
    
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);

    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);

    fprintf('a. binary MSA sum = %f \n', sum(sum(bin_ordered)));
    fprintf('b. posterior mean binary MSA entropy = %f \n', pH);

end


function [ordered] = reorder_nmsa_1(nmsa,nrows,ncols)
%--------------------------------------------------------------------------    
    for k = 1:nrows-1
    nmsa_dif = zeros(nrows,1);
        for i = k:nrows
        nmsa_dif(i) = sum(nmsa(k,:)~=nmsa(i,:));
        end
    [~, sort_msa_ind] = sort(nmsa_dif,'ascend');
    nmsa = nmsa(sort_msa_ind,:);
    end
    ordered = nmsa;
    
end


function [ordered] = reorder_nmsa_2(nmsa,nrows,ncols)
%--------------------------------------------------------------------------
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
    dist(i,i) = 0;
end

ordered = zeros(nrows,ncols);

% The 1st sequence is set by default
ordered(1,:) = nmsa(1,:);
ind = 1;

for i = 2:nrows
    % zero the entire test column
    dist(:,ind) = 0;
    % find the index of the closest sequence
    [~,newind] = max(dist(ind,:));
    ordered(i,:) = nmsa(newind,:);
    % the next index will be the index of the sequence just selected
    ind = newind;
end

end


function [ordered] = reorder_nmsa_3(nmsa,nrows,ncols)
%--------------------------------------------------------------------------
% for i = 1:nrows
%     for j = i:nrows
%         dist(i,j) = sum(nmsa(i,:)~=nmsa(j,:));
%         dist(j,i) = dist(i,j);
%     end
% end

% Here we calculate the distance matrix. The following method is much
% faster than the previous for loop.
bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';
dist = ncols-dist;

ordered = zeros(nrows,ncols);

% The 1st sequence is set by default
ordered(1,:) = nmsa(1,:);
ind = 1;

for i = 2:nrows
    % zero the entire test column
    dist(:,ind) = NaN;
    % find the index of the closest sequence
    [~,newind] = min(dist(ind,:));
    ordered(i,:) = nmsa(newind,:);
    % the next index will be the index of the sequence just selected
    ind = newind;
end

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


function [binmsa] = nmsa_to_binmsa_21q(nmsa)
% Returns each sequence of length L as a vector of size 21L with 0 and 1. 
% Number 21 represents gaps in the Matlab numeric representation of an MSA.

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


function [pMI,mH] = BNMSA_to_MI(bnmsa,bin_ordered_ind,nrows,ncols)
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
pMI = zeros(ncols,ncols);
pIJ = zeros(4,1);
epIJ = zeros(4,1);

pI = zeros(2,ncols);
for I = bin_ordered_ind
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
end
mH = mean(H);

count = 0;
for I = bin_ordered_ind
    count = count +1;
    for J = bin_ordered_ind(count:end)
        pIJ(1) = bnmsa(:,I)'*bnmsa(:,J);
        pIJ(2) = pI(1,I) - pIJ(1);        
        pIJ(3) = pI(1,J) - pIJ(1);
        pIJ(4) = nrows - pIJ(1) - pIJ(2) - pIJ(3);
        epIJ(1) = pI(1,I)*pI(1,J);
        epIJ(2) = pI(1,I)*pI(2,J);
        epIJ(3) = pI(2,I)*pI(1,J);
        epIJ(4) = pI(2,I)*pI(2,J);
        ind = find(pIJ);
        pMI(I,J) = (sum(pIJ(ind) .* (log2(pIJ(ind) ./ epIJ(ind)) + const)))/nrows;
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



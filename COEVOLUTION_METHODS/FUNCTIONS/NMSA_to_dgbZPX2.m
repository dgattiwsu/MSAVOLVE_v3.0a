function [ZPX2] = NMSA_to_dgbZPX2(nmsa)

% dgbZPX2 flowchart
% 1.	Calculate a distance matrix for the sequences of the msa.
% 2.	Reorder the msa by placing as the 1st and 2nd row the two sequences 
%       most similar to each other.
% 3.	Place in the 3rd row the sequence most similar to the sequence in 
%       the 2nd row. 
% 4.	Loop by placing in each consecutive row the sequence most similar 
%       to the sequence in the previous row, until the entire msa is reordered.
% 5.	Convert the msa to binary differential: the 1st row is assigned all 
%       Os. In each consecutive row a 0 is placed at every position at 
%       which the aa is the same as in the previous row, and a 1 at every 
%       position in which the aa is different.
% 6.	Calculate the MI matrix for the binary differential msa.
% 7.	Convert each row of the reordered msa to long binary format in 
%       which every symbol (including gaps) is represented by a vector of 
%       twenty 0s and a single 1.
% 8.	Convert the long binary format msa to binary differential as in step 5.
% 9.	Calculate the positional covariance matrix (largeCOV) of the msa in 
%       differential long binary format.
% 10.	Calculate the Frobenius norm of every submatrix of largeCOV that 
%       corresponds to the covariance between two columns of the reordered 
%       msa. This generates a positional covariance matrix (COV) of the 
%       reordered msa.
% 11.	Scale MI and COV matrices by linear regression.
% 12.	Sum scaled MI and COV matrices.
% 13.	Calculate a ZPX2 matrix from the summed matrix.
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
    
%       Here we find out which columns are not all 0s.

    bin_ordered_sum = sum(bin_ordered);
    bin_ordered_ind = find(bin_ordered_sum);
    
% 6.	Here we calculate the MI matrix for the binary differential msa. We
%       use a special function tailored for binary arrays and we only need
%       to do the calculation for the columns which are not all 0s.

    [MI,pH] = BNMSA_to_MI(bin_ordered,bin_ordered_ind,nrows,ncols);
    
% 7.	Here we convert each row of the reordered msa to long binary format  
%       in which every symbol (including gaps) is represented by a vector  
%       of twenty 0s and a single 1.

    ordered = nmsa_to_binmsa_21q(ordered);

% 8.    Here we convert the long binary to binary differential.

    [nrows,ncols] = size(ordered);
    bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
    bin_ordered_sum = sum(bin_ordered);
    bin_ordered_ind = find(bin_ordered_sum);

% 9.	Here we calculate the positional covariance matrix (largeCOV) of  
%       the msa in differential long binary format. In practice we calculate  
%       the cov matrix only for the indices included in the vector 
%       'bin_ordered_ind'.

    gCOV = zeros(ncols,ncols);
    [ogCOV] = cov(bin_ordered(:,bin_ordered_ind));
    gCOV(bin_ordered_ind,bin_ordered_ind) = ogCOV;
    
% 10.	Here we calculate the Frobenius norm of every submatrix of largeCOV  
%       that corresponds to the covariance between two columns of the reordered 
%       msa. This generates a positional covariance matrix (COV) of the 
%       reordered msa.

%   First we zero the diagonal of the large covariance matrix

        for i = 1:ncols
            gCOV(i,i) = 0;
        end

%   Here we recover the original values of nrows and ncols

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

% 11.	Here we scale MI and COV matrices by linear regression.
    
    scale = scale_matrices(COV,MI);

    fprintf('best scale MI/COV = %f \n', scale);

% 12.	Here we sum the scaled MI and COV matrices.

    MI = MI + scale*COV;
        
% 13.	Here we calculate a ZPX2 matrix from the summed matrix.
    
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
    fprintf('a. long binary MSA sum = %f \n', sum(sum(bin_ordered)));
    fprintf('b. posterior mean binary MSA entropy = %f \n', pH);

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
        % H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
        H(I) = const - (sum(pI(ind,I).*(log2(pI(ind,I)))))/nrows;
        
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
        % pMI(I,J) = (sum(pIJ(ind) .* (log2(pIJ(ind) ./ epIJ(ind)) + const)))/nrows;
        pMI(I,J) = const + (sum(pIJ(ind) .* (log2(pIJ(ind) ./ epIJ(ind)))))/nrows;
    end
end
pMI = pMI+pMI';
for m = 1:ncols
pMI(m,m)=NaN;
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
% would give the wrong MI. Change of sign is not in the original ZPX2 
% algorithm by Gloor, but is included in the ZRES algorithm by Chen. 
    
    if (ZPX2_i<0&&ZPX2_j<0)
        pZPX2(m,n)=-pZPX2(m,n);
    end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end



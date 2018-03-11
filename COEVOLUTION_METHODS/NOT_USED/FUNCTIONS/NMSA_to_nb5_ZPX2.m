function [ZPX2,MIP,MI] = NMSA_to_nb5_ZPX2(nmsa,Lmat,lambda,delta,method)
% This function differs from 'nbZPX2' for ignoring gaps in the calculation 
% of the distance matrix. Furthermore it converts the final Mi matrix to a
% matrix of partial correlation coefficients.

[nrows,ncols] = size(nmsa);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        not_gap = nmsa(i,:) ~= 25 & nmsa(j,:) ~= 25;
        % dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
        dist(i,j) = sum(nmsa(i,not_gap)==nmsa(j,not_gap));
        dist(j,i) = dist(i,j);
    end
end

triu_dist = triu(dist,1);
[~,max_ind] = max(triu_dist(:));
[firstrow,firstcol] = ind2sub([nrows,nrows],max_ind);

% Here we select as the top sequence the sequence with the smallest
% entropy. Not clear whether the following 4 lines make a difference.
% rowsum = Entropy(nmsa(firstrow,:)');
% colsum = Entropy(nmsa(firstcol,:)');
% if rowsum >= colsum
%     firstrow = firstcol;
% end

nmsa1 = nmsa(firstrow:nrows,:);
nmsa2 = nmsa(1:firstrow-1,:);
nmsa = [nmsa1;nmsa2];

%--------------------------------------------------------------------------
% Here we reorder based on the first algorithm.
ordered = reorder_nmsa_1(nmsa,nrows,ncols);
%--------------------------------------------------------------------------

% Here we convert the msa from matlab numeric format to binary
% differential.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end

%Here we set to zero all the indices of the reordered nmsa that are 0 in 
% the binary msa. 
z_ordered = ordered;
z_ordered(~bin_ordered) = 0;

    % MI matrix
    [MI] = NMSA_to_MI(z_ordered);
    
%--------------------------------------------------------------------------    
% Here we determine the partial correlation matrix using the 'QUIC' method.
% Recommended input:
% Lmat = 0.02
% lambda = 0.0
% delta = 0.0001
% method = 'QUIC'

[RHO] = COV_to_RHO(MI,ncols,Lmat,lambda,delta,method);

%--------------------------------------------------------------------------
    
    
    % MIP matrix
    MIP = MI_to_MIP(RHO,ncols);
    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
%    fprintf('a. binary MSA sum = %f \n', sum(sum(bin_ordered)));
%    fprintf('b. posterior mean binary MSA entropy = %f \n', pH);
            
end


function [ordered] = reorder_nmsa_1(nmsa,nrows,ncols)
%--------------------------------------------------------------------------    
    for k = 1:nrows-1
    nmsa_dif = zeros(nrows,1);
        for i = k:nrows
        not_gap = nmsa(k,:) ~= 25 & nmsa(i,:) ~= 25;
        % nmsa_dif(i) = sum(nmsa(k,:)~=nmsa(i,:));
        nmsa_dif(i) = sum(nmsa(k,not_gap)~=nmsa(i,not_gap));
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
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)~=nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
end

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
% This function returns the best scale for mat1/mat2.
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


function [ MI_mat ] = NMSA_to_MI( nmsa )
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (.nmsa). The function uses  the "MutualInformation" function from 
% Dwinnel package. It differs from 'NMSA_to_gcMI' because it does not apply 
% a correction for gaps. It is the same function as NMSA_to_ngcMI, where
% ngc = no gap correction.

[~,ncols] = size(nmsa);
MI_mat = zeros(ncols,ncols);

for i = 1:ncols
    for j = i:ncols
    MI_mat(i,j) = MutualInformation(nmsa(:,i),nmsa(:,j));
    MI_mat(j,i) = MI_mat(i,j);    
    end
    % MI_mat(i,i)=NaN;
end

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
pMI = zeros(ncols,ncols);
pIJ = zeros(4,1);
epIJ = zeros(4,1);

pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
end
mH = mean(H);

for I = 1:ncols
    for J = I:ncols
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


function [pIJ_mat] = BNMSA_to_pIJ(bnmsa,bin_ordered_ind,nrows,ncols)
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
pIJ_mat = zeros(ncols,ncols);
% pIJ = zeros(4,1);

pI = zeros(1,ncols);
for I = bin_ordered_ind
        pI(I) = nnz(bnmsa(:,I));
end

count = 0;
for I = bin_ordered_ind
    count = count +1;
    for J = bin_ordered_ind(count:end)
        pIJ = bnmsa(:,I)'*bnmsa(:,J);
        if pIJ == 0
        else
            epIJ = pI(I)*pI(J);
            pIJ_mat(I,J) = (pIJ * (log2(pIJ / epIJ) + const))/nrows;
        end
        
    end
end
pIJ_mat = pIJ_mat+pIJ_mat';
for m = 1:ncols
pIJ_mat(m,m)=NaN;
end

end


function [RHO,opt,iter] = COV_to_RHO(COV,ncols,Lmat,lambda,delta,method)
% This function calculates the partial correlation matrix for the columns
% of an MSA. The standard 'inv', QUIC, or GLASSO methods can be used to
% invert the covariance matrix. Lmat is a value between 1 and 0 controlling
% the 'sparseness' of the concentration matrix (smaller values produce a
% less sparse inverse, and take more time. Larger values increase the
% sparseness, with the inverse ultimately being populated only in the
% diagonal). The default value of Lmat is 0.04. Lmat can also be a
% regularization matrix of the same dimensions as the covariance matrix to
% be inverted.
    
    imatrix = eye(ncols);
    fmatrix = mean(diag(COV))*imatrix;
    
    % Test for positive definitiveness
    sCOV = COV;
    while lambda <= 1.0
        try
            sCOV = lambda*fmatrix + (1-lambda)*COV;    
            test_chol = chol(sCOV);
        end
        if exist('test_chol','var')
            break
        end
        lambda = lambda + delta;
        % sCOV = lambda*fmatrix + (1-lambda)*COV;    
    end

    % sCOV = lambda*fmatrix + (1-lambda)*COV;    
    % sCOV = -lambda*fmatrix + lambda*COV;    
    % sCOV = COV;
    % Lmat = 0.06 % default value
    
    if exist('Lmat','var')
    else
        Lmat = 0.04;
    end
    
        switch method
        
            case 'INVERSE'
            cov_inv = inv(sCOV);
            cov_rev = sCOV;
                    
            case 'QUIC'                
            [cov_inv cov_rev opt cputime iter dGap] = ...
                QUIC('default', sCOV, Lmat, 1e-6, 2, 200);

            case 'GLASSO'
            % function [w, theta, iter, avgTol, hasError] = ...
            % glasso(numVars, s, computePath, lambda, approximate, ...
            % warmInit, verbose, penalDiag, tolThreshold, ...
            % maxIter, w, theta)
            [cov_rev, cov_inv, iter, avgTol, hasError] = ...
                glasso(ncols, sCOV, 0, Lmat.*ones(ncols), 0, ...
                0, 0, 1, 1e-4, 1e4, zeros(ncols), zeros(ncols));
            
        end
           
    RHO = zeros(ncols,ncols);
    for i = 1:ncols
        for j = i:ncols
            RHO(i,j) = -cov_inv(i,j)/sqrt(cov_inv(i,i)*cov_inv(j,j));
            RHO(j,i) = RHO(i,j);
        end
    end
    RHO(i,i) = NaN;    
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



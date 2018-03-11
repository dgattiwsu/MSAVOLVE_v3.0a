function [pcZPX2,pcMIP,cov_inv] = ...
NMSA_to_pc2_ZPX2(nmsa,Lmat,lambda,delta,inverse_method,pc_method,nsymbols)
% Differs from NMSA_to_pc1a_ZPX2 and pc1b_ZPX2 by the fact that it
% calculates the sparse inverse on covariance matrix of the binary 'long'
% format, rather than on the recovered shrunk covariance matrix.
%
% This function calculates the partial correlation matrix for the columns
% of an MSA. The MSA is first converted to its binary 'long' representation
% in which every column is extended into 21 different columns, each one
% representing a different aa. Then, the covariance matrix is calculated,
% and an inverse method is used to invert it. Three option are available
% for this purpose: standard 'INVERSE', 'GLASSO', or 'QUIC'. Lmat is a
% value between 1 and 0 controlling the 'sparseness' of the concentration
% matrix (smaller values produce a less sparse inverse, and take more time.
% Larger values increase the sparseness, with the inverse ultimately being
% populated only in the diagonal). The default value of Lmat is 0.02. Lmat
% can also be a regularization matrix of the same dimensions as the
% covariance matrix to be inverted. The value of nsymbols determines
% whether we include or not gaps in the calculation of the 'shrunk' inverse
% covariance matrix (nsymbols = 21 includes gaps; nsymbols = 20 excludes
% gaps). 'pc_method' determines whether we calculate the MIP matrix from
% the inverse covariance or from the RHO matrix. It affects only the sign
% of the MIP matrix, but not that of the ZPX2 matrix. Recommended usage:
% [pcZPX2] = NMSA_to_pcZPX2(nmsa,0.02,0.0,0.0001,'QUIC','RHO',20);

bin_ordered = nmsa_to_binmsa_21q(nmsa);

[~,ncols] = size(bin_ordered);
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);

% Here we calculate the COV matrix only for the indices included in the
% vector 'bin_ordered_ind'.

    COV = zeros(ncols,ncols);
    ogCOV = cov(bin_ordered(:,bin_ordered_ind));
    COV(bin_ordered_ind,bin_ordered_ind) = ogCOV;

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
    end

    % Lmat = 0.02 % default value
    
    if exist('Lmat','var')
        % Lmat = Lmat*NMSA_to_MI(nmsa) + eye(ncols);
    else
        Lmat = 0.02;
    end
    
        switch inverse_method
        
            case 'INVERSE'
            gcov_inv = inv(sCOV);
            gcov_rev = sCOV;
                    
            case 'QUIC'                
            [gcov_inv gcov_rev opt cputime iter dGap] = ...
                QUIC('default', sCOV, Lmat, 1e-6, 2, 200);

            case 'GLASSO'
            % function [w, theta, iter, avgTol, hasError] = ...
            % glasso(numVars, s, computePath, lambda, approximate, ...
            % warmInit, verbose, penalDiag, tolThreshold, ...
            % maxIter, w, theta)
            [gcov_rev, gcov_inv, iter, avgTol, hasError] = ...
                glasso(ncols, sCOV, 0, Lmat.*ones(ncols), 0, ...
                0, 0, 1, 1e-4, 1e4, zeros(ncols), zeros(ncols));
            
        end
    
        
    % Here we recover the original values of nrows and ncols
    shift = nsymbols - 1;        
    [~,ncols] = size(nmsa);
    cov_inv = zeros(ncols,ncols);
    for oi = 1:ncols
        i = 21*(oi-1) + 1;
        for oj = oi:ncols
            j = 21*(oj-1)+1;
            % Frobenius norm
            % cov_inv(oi,oj) = norm(gcov_inv(i:i+shift,j:j+shift),'fro');
            % cov_inv(oj,oi) = cov_inv(oi,oj);       
            % L1 norm
            cov_inv(oi,oj) = norm(gcov_inv(i:i+shift,j:j+shift),1);
            cov_inv(oj,oi) = cov_inv(oi,oj);       
           
        end
    end
   
    rho = zeros(ncols,ncols);
    for i = 1:ncols
        for j = i:ncols
            rho(i,j) = -cov_inv(i,j)/sqrt(cov_inv(i,i)*cov_inv(j,j));
            rho(j,i) = rho(i,j);
        end
    end
           
    % MIP matrix
    switch pc_method
        case 'RHO'
        pcMIP = MI_to_MIP(rho,ncols);
        case 'INVERSE'
        pcMIP = MI_to_MIP(cov_inv,ncols);
    end

    % ZPX2 matrix
    pcZPX2 = MIP_to_ZPX2(pcMIP,ncols);
    
    pcZPX2 = real(pcZPX2);
    
end


function [pcorr] = pcorr(nmsa,ncols)
% This function calculate a matrix of partial correlation coefficients for
% the columns of the binary msa.
nmsa = +nmsa;
sum_col_vec = sum(nmsa,2);
nC = sum_col_vec/norm(sum_col_vec);
% cC = sum_col_vec - mean(sum_col_vec);
% nC = cC/norm(cC);
% pcorr = NaN(ncols,ncols);
pcorr = zeros(ncols,ncols);
for i = 1:ncols
    A = nmsa(:,i);
    pA = A - (A'*nC)*nC;
    cpA = pA - mean(pA);
    ncpA = cpA/norm(cpA);
    for j = i:ncols
        B = nmsa(:,j);
        pB = B - (B'*nC)*nC;
        cpB = pB - mean(pB);
        ncpB = cpB/norm(cpB);
        pcorr(i,j) = ncpA'*ncpB;
        pcorr(j,i) = pcorr(i,j);
    end
end 
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
    
%     if (ZPX2_i<0&&ZPX2_j<0)
%         pZPX2(m,n)=-pZPX2(m,n);
%     end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end



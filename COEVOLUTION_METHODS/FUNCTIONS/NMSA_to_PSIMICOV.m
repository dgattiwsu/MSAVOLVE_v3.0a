function [pcZPX2,pcMIP,ppvZPX2,ppvMIP,gCOV,cov_inv,rho,COV,sCOV,lambda] = ...
    NMSA_to_PSIMICOV(nmsa,long_method,cov_diag,recov_method,...
    Lmat,lambda,delta,inverse_method,pc_method,nsymbols)
% Matlab experimental implementations of the PSICOV algorithm. It
% calculates a partial correlation matrix for the columns of an MSA. The
% MSA is first converted to its binary 'long' representation in which every
% column is extended into 20 (dist_method = 'NOGAPS) or 21 (dist_method =
% 'GAPS) different columns, each one representing a different aa. A
% threshold is applied for the similarity between sequences. A pseudocount
% (psc_lambda) is also applied. Then, the covariance matrix is calculated
% and is converted to the equivalent covariance matrix for the 'short'
% standard representation. Before this conversion it is possible to keep
% (cov_diag = 'KEEP') or zero (cov_diag = 'ZERO') [recommended] the
% diagonal of the large covariance matrix. The conversion to the 'short'
% format can be carried out by taking the Frobenius norm (recov_method =
% 'FRO') or the L1 norm (recov_method = 'L1') of each submatrix. The value
% of nsymbols determines whether we include or not gaps in the calculation
% of the 'short' covariance matrix (nsymbols = 21 includes gaps; nsymbols =
% 20 excludes gaps). Finally, an inverse method is used to invert the
% covariance matrix. Three option are available for this purpose: standard
% 'INVERSE', 'GLASSO', or 'QUIC'. Lmat is a value between 1 and 0
% controlling the 'sparseness' of the concentration matrix (smaller values
% produce a less sparse inverse, and take more time. Larger values increase
% the sparseness, with the inverse ultimately being populated only in the
% diagonal). The default value of Lmat is 0.02. Suitable ranges are between
% 0.025 and 0.005. Lmat can also be a regularization matrix of the same
% dimensions as the covariance matrix to be inverted. 'pc_method'
% determines whether we calculate the MIP matrix from the inverse
% covariance or from the RHO matrix. It appears to affects only the sign of
% the MIP matrix, but not that of the ZPX2 matrix. The ppvZPX2 and ppvMIP
% matrices include the correction produced by the application of a logistic
% fit to the MIP data. Possible usage: [pcZPX2,pcMIP,ppvZPX2,ppvMIP] =
% NMSA_to_PSICOV(nmsa,'NOGAPS',0.62,1,'ZERO','FRO',0.005,0.0,0.0001,'QUIC','RHO',20);
% To reproduce PSICOV use: [pcZPX2,pcMIP,ppvZPX2,ppvMIP] =
% NMSA_to_PSICOV(nmsa,'NOGAPS',0.62,1,'ZERO','FRO',0.015,0.0,0.0001,'GLASSO','INVERSE',20);

switch long_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(nmsa);
        jump = 20;
        nsymbols = 20; % value is changed to 20 even if 21 is requested.
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        jump = 21;
end

[nrows,ncols] = size(bin_ordered);

% Here we find the indices of the non-zero columns.
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);
n_bin_ordered_ind = length(bin_ordered_ind);
s_bin_ordered = bin_ordered(:,bin_ordered_ind);
gMI = zeros(ncols,ncols);
gCOV = zeros(ncols,ncols);

% MI matrix
ogMI = BNMSA_to_MI(s_bin_ordered,nrows,n_bin_ordered_ind);
gMI(bin_ordered_ind,bin_ordered_ind) = ogMI;

% Covariance matrix
ogCOV = cov(s_bin_ordered);
gCOV(bin_ordered_ind,bin_ordered_ind) = ogCOV;

% Here we zero the diagonal of the large covariance matrix

switch cov_diag
    case 'ZERO'
        for i = 1:ncols
            gMI(i,i) = 0;
            gCOV(i,i) = 0;
        end
    case 'KEEP'
end
        
    % Here we recover the original values of nrows and ncols
    shift = nsymbols - 1;    
    [~,ncols] = size(nmsa);
    MI = zeros(ncols,ncols);
    COV = zeros(ncols,ncols);
    for oi = 1:ncols
        i = jump*(oi-1) + 1;
        for oj = oi:ncols
            j = jump*(oj-1)+1;
            
            % Here we define the submatrix 
            submi = gMI(i:i+shift,j:j+shift);
            subcov = gCOV(i:i+shift,j:j+shift);
            
            % and zero its diagonal
            % for k = 1:nsymbols
            %     subcov(k,k) = 0;
            % end
            
            switch recov_method
                case 'FRO'
                % Frobenius norm
                MI(oi,oj) = norm(submi,'fro');
                COV(oi,oj) = norm(subcov,'fro');
                case 'SPECTRAL'
                % Spectral norm
                [~,MI(oi,oj),~] =svds(submi,1);
                [~,COV(oi,oj),~] =svds(subcov,1);
                case 'L1'
                % L1 norm
                MI(oi,oj) = norm(submi,1);
                COV(oi,oj) = norm(subcov,1);
                case '2'
                % 2 norm
                MI(oi,oj) = norm(submi,2);
                COV(oi,oj) = norm(subcov,2);
                case 'SUM'
                % Simple sum
                MI(oi,oj) = sum(sum(submi));
                COV(oi,oj) = sum(sum(subcov));
            end                       
            MI(oj,oi) = MI(oi,oj);
            COV(oj,oi) = COV(oi,oj);
        end
    end
    
    
% Here we scale MI and COV matrices by linear regression.
    
    scale = scale_matrices(COV,MI);

    fprintf('best scale MI/COV = %f \n', scale);

% Here we sum the scaled MI and COV matrices.

    COV = MI + scale*COV;
       
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
        % Here Lmat can be a 'regularization' Rij matrix that gives more
        % weight to pair known to be important.
        % Lmat = Lmat*NMSA_to_MI(nmsa) + imatrix;
        % Lmat = Lmat*NMSA_to_MI(nmsa);        
        % Lmat = Lmat*COV + imatrix;
    else
        Lmat = 0.02;
    end
    
        switch inverse_method
        
            case 'INVERSE'
            cov_inv = inv(sCOV);
            cov_rev = sCOV;
                    
            case 'QUIC'                
            [cov_inv cov_rev opt cputime iter dGap] = ...
                QUIC('default', sCOV, Lmat, 1e-6, 2, 200);

        end
    
    nonzeros = (nnz(cov_inv)-ncols)/(ncols*(ncols-1));
    fprintf('Sparsity = %f \n', nonzeros);

        
    rho = zeros(ncols,ncols);
    for i = 1:ncols
        for j = i:ncols
            rho(i,j) = -cov_inv(i,j)/sqrt(cov_inv(i,i)*cov_inv(j,j));
            rho(j,i) = rho(i,j);
        end
    end
    
    % MIP matrix
    n_rho = rho;
    n_cov_inv = cov_inv;
    
    for i = 1:ncols
        n_rho(i,i) = NaN;
        n_cov_inv(i,i) = NaN;
    end
    
    switch pc_method
        case 'RHO'
                pcMIP = MI_to_MIP(n_rho,ncols);
        case 'INVERSE'
                pcMIP = MI_to_MIP(-n_cov_inv,ncols);
    end

    ppvMIP = mat_to_ppv(pcMIP);
    
    % ZPX2 matrix
    pcZPX2 = MIP_to_ZPX2(pcMIP,ncols);    
    ppvZPX2 = mat_to_ppv(pcZPX2);
    
    
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
        pMI(J,I) = pMI(I,J);
    end
end

% for m = 1:ncols
% pMI(m,m)=NaN;
% pMI(m,m)=0;
% end

end


function [scale] = scale_matrices(mat1,mat2)
%--------------------------------------------------------------------------
% This function returns the best scale such that scale*mat1 = mat2.
[rows,cols] = size(mat1);
template = ones(rows,cols);
upper = triu(template,1);
triu_ind = template == upper;
mat1 = mat1(triu_ind);
mat2 = mat2(triu_ind);
% To find the best scale we use the Matlab 'backslash' operator. However,
% notice that exactly the same result can be obtained as a simple ratio of
% dot products, which is the standard way of scaling two vectors.
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


function [ppv_mat] = mat_to_ppv(mat)

% Here we fit a logistic distribution to the data.

[~,cols] = size(mat);
data = mat(:);
ndata = length(data);
all_ind = 1:ndata;

% Logistic fit Matlab style
% [param] = fitdist(data,'logistic');
% data_mean = param.Params(1);
% data_std = param.Params(2);
% z_data = (data - data_mean)/data_std;
% ppv_data = z_data;

% Logistic fit PSICOV style
data_mean = nanmean(data);
data_std = nanstd(data);
z_data = (data - data_mean)/data_std; 
ppv_data = 0.904 ./ (1.0 + 16.61 * exp(-0.8105 * z_data));

ppv_mat = zeros(cols);
[i,j] = ind2sub(size(mat),all_ind);

for n = 1:ndata
    ppv_mat(i(n),j(n)) = ppv_data(n);
end

end

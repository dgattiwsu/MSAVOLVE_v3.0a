function [pcZPX2,pcMIP,ppvZPX2,ppvMIP,ZPX2,MIP,cov_inv,rho,COV,sCOV,lambda] = ...
    NMSA_to_vertPSIMI(nmsa,vert_method,Lmat,lambda,delta,...
    inverse_method,pc_method)
% Matlab experimental implementations of a 'vertical' PSIMI algorithm. 
% It calculates a fast binary MI matrix for the columns of an MSA. The
% MSA is first converted to a binary 'long' representation in which every
% column is extended into a longer column with 20 (vert_method = 'NOGAPS) or 
% 21 (vert_method = 'GAPS) times the original number of rows. Finally, an 
% inverse method is used to invert the MI matrix. Two option are available 
% for this purpose: standard 'INVERSE', or 'QUIC'. Lmat is a value between 
% 1 and 0 controlling the 'sparseness' of the concentration matrix (smaller 
% values produce a less sparse inverse, and take more time. Larger values 
% increase the sparseness, with the inverse ultimately being populated only 
% in the diagonal). Suitable ranges for Lmat are between 0.025 and 0.005. 
% Lmat can also be a regularization matrix of the same dimensions as the 
% MI matrix to be inverted. 'pc_method' determines whether we calculate the 
% MIP matrix from the inverse covariance or from the RHO matrix. The ppvZPX2 
% and ppvMIP matrices include the correction produced by the application of 
% a logistic fit. Possible usage: [pcZPX2,pcMIP,ppvZPX2,ppvMIP,ZPX2,MIP] =
% NMSA_to_PSICOV(nmsa,'NOGAPS',0.01,0.0,0.0001,'QUIC','RHO');

switch vert_method
    case 'NOGAPS'
        v_bin_ordered = nmsa_to_v_binmsa_20q(nmsa);
        
    case 'GAPS'
        v_bin_ordered = nmsa_to_v_binmsa_21q(nmsa);
end

[v_nrows,v_ncols] = size(v_bin_ordered);

% Here we determine the MI, MIP, and ZPX2 matrix before any weighing for 
% similarity between sequences.

MI = BNMSA_to_MI(v_bin_ordered,v_nrows,v_ncols);
MIP = MI_to_MIP(MI,v_ncols);
ZPX2 = MIP_to_ZPX2(MIP,v_ncols);

COV = MI;

    imatrix = eye(v_ncols);
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
    
    nonzeros = (nnz(cov_inv)-v_ncols)/(v_ncols*(v_ncols-1));
    fprintf('Sparsity = %f \n', nonzeros);

        
    rho = zeros(v_ncols,v_ncols);
    for i = 1:v_ncols
        for j = i:v_ncols
            rho(i,j) = -cov_inv(i,j)/sqrt(cov_inv(i,i)*cov_inv(j,j));
            rho(j,i) = rho(i,j);
        end
    end
    
    % MIP matrix
    n_rho = rho;
    n_cov_inv = cov_inv;
    
    for i = 1:v_ncols
        n_rho(i,i) = NaN;
        n_cov_inv(i,i) = NaN;
    end
    
    switch pc_method
        case 'RHO'
                pcMIP = MI_to_MIP(n_rho,v_ncols);
        case 'INVERSE'
                pcMIP = MI_to_MIP(-n_cov_inv,v_ncols);
    end

    ppvMIP = mat_to_ppv(pcMIP);
    
    % ZPX2 matrix
    pcZPX2 = MIP_to_ZPX2(pcMIP,v_ncols);        
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

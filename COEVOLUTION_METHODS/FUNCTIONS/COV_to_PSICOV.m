function [pcZPX2,pcMIP,ppvZPX2,ppvMIP,cov_inv,rho,COV,sCOV,lambda] = ...
    COV_to_PSICOV(COV,Lmat,lambda,delta,inverse_method,pc_method)
% This function uses a sparse inverse method to invert a covariance 
% matrix and ultimately obtain a partial correlation matrix. Three option
% are available for this purpose: standard 'INVERSE', 'GLASSO', or 'QUIC'.
% Lmat is a value between 1 and 0 controlling the 'sparseness' of the
% concentration matrix (smaller values produce a less sparse inverse, and
% take more time. Larger values increase the sparseness, with the inverse
% ultimately being populated only in the diagonal). The default value of
% Lmat is 0.02. Suitable ranges are between 0.025 and 0.005. Lmat can also
% be a regularization matrix of the same dimensions as the covariance
% matrix to be inverted. 'lambda' and 'delta' define the conditioning of
% the covariance matrix prior to the calculation of the sparse inverse. The
% recommended values are lambda = 0.0 and delta = 0.0001. Larger values of
% lambda and delta will make the calculation of the inverse faster, but
% possibly less accurate. 'pc_method' determines whether we calculate the
% MIP matrix from the inverse covariance or from the RHO matrix. It appears
% to affects only the sign of the MIP matrix, but not that of the ZPX2
% matrix. The ppvZPX2 and ppvMIP matrices include the correction produced
% by the application of a logistic fit to the MIP data.
% Example syntax:
% [pcZPX2,~,ppvZPX2,~] = COV_to_PSICOV(COV,0.005,0.0,0.0001,'QUIC','RHO');

% Here we test for positive definitiveness
    [~,ncols] = size(COV);
    imatrix = eye(ncols);
    fmatrix = mean(diag(COV))*imatrix;
    
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
        % weight to pairs known to be important.
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

            case 'GLASSO'
            % function [w, theta, iter, avgTol, hasError] = ...
            % glasso(numVars, s, computePath, lambda, approximate, ...
            % warmInit, verbose, penalDiag, tolThreshold, ...
            % maxIter, w, theta)
            [cov_rev, cov_inv, iter, avgTol, hasError] = ...
                glasso(ncols, sCOV, 0, Lmat.*ones(ncols), 0, ...
                0, 0, 1, 1e-4, 1e4, zeros(ncols), zeros(ncols));
            
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
    pcZPX2 = real(pcZPX2);        
    ppvZPX2 = mat_to_ppv(pcZPX2);
        
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



function [pZP,pZP2] = MIP_to_ZP(pMIP,ncols)
%--------------------------------------------------------------------------
% ZP calculation based on Dickson,...,Gloor (PLoS ONE 5(6):e11082,2012).
% ZP2 is ZP squared so that it can be compared directly to ZPX2.

pZP=zeros(ncols,ncols);
pZP2=zeros(ncols,ncols);

mean_mat = nanmean(pMIP(:));
std_mat = nanstd(pMIP(:));
    
for i=1:ncols
    for j=i:ncols

    pZP(i,j) = (pMIP(i,j)-mean_mat)/std_mat;
    pZP2(i,j) = pZP(i,j)^2;
    pZP(j,i) = pZP(j,i);        
    pZP2(j,i) = pZP2(j,i);        
    
    end
    
    pZP(m,m)=NaN;
    pZP2(m,m)=NaN;

end
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

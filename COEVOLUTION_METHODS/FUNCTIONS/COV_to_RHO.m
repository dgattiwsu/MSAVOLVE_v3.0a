function [rho,cov_inv,COV,sCOV,lambda] = ...
    COV_to_RHO(COV,Lmat,lambda,delta,inverse_method)
% This function uses a sparse inverse method to invert a covariance matrix
% and ultimately obtain a partial correlation matrix or an inverse matrix
% suitable for further calculation of MIP and ZPX2 values. Three option are
% available for this purpose: standard 'INVERSE', 'GLASSO', or 'QUIC'. Lmat
% is a value between 1 and 0 controlling the 'sparseness' of the
% concentration matrix (smaller values produce a less sparse inverse, and
% take more time. Larger values increase the sparseness, with the inverse
% ultimately being populated only in the diagonal). The default value of
% Lmat is 0.02. Suitable ranges are between 0.025 and 0.005. Lmat can also
% be a regularization matrix of the same dimensions as the covariance
% matrix to be inverted. 'lambda' and 'delta' define the conditioning of
% the covariance matrix prior to the calculation of the sparse inverse. The
% recommended values are lambda = 0.0 and delta = 0.0001. Larger values of
% lambda and delta will make the calculation of the inverse faster, but
% possibly less accurate. The sign of all the 'cov_inv' matrix terms is
% changed to make it compatible for further MIP/ZPX2 calculations.
% Example syntax:
% [rho,cov_inv] = COV_to_RHO(COV,0.005,0.0,0.0001,'QUIC');

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
    
    cov_inv = -cov_inv;
    
    for i = 1:ncols
        rho(i,i) = NaN;
        cov_inv(i,i) = NaN;
    end
    
end



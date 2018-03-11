function [ pcZPX2,pcMIP,Fi,fnmsa,COV,sCOV,cov_inv,rho,lambda ] = ...
    NMSA_to_fpc2_ZPX2(nmsa,Lmat,lambda,delta,inverse_method,pc_method,nsymbols)
% This function calculates a partial correlation matrix (RHO) without first 
% converting the MSA to a long binary format. This is achieved by
% converting every symbol in the msa to the frequency of that symbol in
% each column. It differs from fpcZPX2 as no weights are applied to the
% frequencies to correct for the similarities between sequences. Since it
% does not calculate a distance matrix, it is slightly faster for large msa's.

[nrows,ncols] = size(nmsa);

% Replace unusual symbols
ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

Fi = zeros(21,ncols);
nmsa_3 = zeros(nrows,ncols,21);

% Here we create a 3d matrix in which every layer has the dimensions of the
% nmsa and represents one of 20 symbols (if gaps are excluded: nsymbols = 20). 
% Then every row of each layer is scaled by the weight of that sequence.

for i = 1:nsymbols
    nmsa_3(:,:,i) = nmsa == i;
    for j = 1:nrows
        nmsa_3(j,:,i) = nmsa_3(j,:,i);
    end
end

% nmsa_1 = sum(nmsa_3,3);

for i = 1:21
    layer = squeeze(nmsa_3(:,:,i));
    Fi(i,:) = sum(layer)/nrows;
end

% After removing the frequencies of the gaps we scale the remaining
% frequencies.
if nsymbols == 20
sumFi = sum(Fi);
scaleFi = repmat(sumFi,21,1);
Fi = Fi./scaleFi;
% sumFi = sum(Fi);
end

% Here we create a new msa in which every symbol is replaced by the
% frequence of that symbol in that column of the msa.

fnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        fnmsa(i,j) = Fi(row(j),j);
        end
end

% MI = NMSA_to_MI(fnmsa);

COV = cov(fnmsa);

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
        pcMIP = MI_to_MIP(n_rho);
        case 'INVERSE'
        pcMIP = MI_to_MIP(n_cov_inv);
    end
    
    % ZPX2 matrix
    [~,pcZPX2] = MI_to_ZPX(pcMIP);
    
    pcZPX2 = real(pcZPX2);
    
        
end


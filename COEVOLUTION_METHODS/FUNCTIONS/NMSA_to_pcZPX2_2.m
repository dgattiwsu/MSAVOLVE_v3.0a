function [pcZPX2,pcMIP,cov_inv,rho,COV,sCOV,lambda] = ...
    NMSA_to_pcZPX2_2(nmsa,threshold,psc_lambda,...
    Lmat,lambda,delta,inverse_method,pc_method,nsymbols)
% This function is the 1st of two alternative implementations of the PSICOV
% algorithm. It calculates a partial correlation matrix for the columns of
% an MSA. The MSA is first converted to its binary 'long' representation in
% which every column is extended into 21 different columns, each one
% representing a different aa. A threshold is applied for the similarity
% between sequences. A pseudocount (psc_lambda) is also applied. Then, the
% covariance matrix is calculated and is shrunk to the equivalent
% covariance matrix for the 'short' standard representation. The value of
% nsymbols determines whether we include or not gaps in the calculation of
% the 'shrunk' covariance matrix (nsymbols = 21 includes gaps; nsymbols =
% 20 excludes gaps). Finally, an inverse method is used to invert the
% covariance matrix. Three option are available for this purpose: standard
% 'INVERSE', 'GLASSO', or 'QUIC'. Lmat is a value between 1 and 0
% controlling the 'sparseness' of the concentration matrix (smaller values
% produce a less sparse inverse, and take more time. Larger values increase
% the sparseness, with the inverse ultimately being populated only in the
% diagonal). The default value of Lmat is 0.02. Suitable ranges are between
% 0.015 and 0.025. Lmat can also be a regularization matrix of the same
% dimensions as the covariance matrix to be inverted. 'pc_method'
% determines whether we calculate the MIP matrix from the inverse
% covariance or from the RHO matrix. It affects only the sign of the MIP
% matrix, but not that of the ZPX2 matrix. Recommended usage:
% [pcZPX2,pcMIP] =
% NMSA_to_pcZPX2(nmsa,0.8,1,0.015,0.0,0.0001,'QUIC','RHO',20);

bin_ordered = nmsa_to_binmsa_21q(nmsa);
[~,ncols] = size(bin_ordered);

% Here we calculate the distance matrix on the transpose and apply the
% weights for the similarity between sequences.
if threshold < 1.0
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);
n_bin_ordered_ind = length(bin_ordered_ind);

dist = corr(bin_ordered(:,bin_ordered_ind)');
[W,Meff] = nmsa_to_W(nmsa,dist,threshold);
W_mat = repmat(W,1,length(bin_ordered_ind));
loq = psc_lambda/21;
loq2 = psc_lambda/21^2;
l_Meff = psc_lambda + Meff;
w_bin_ordered = W_mat.*bin_ordered(:,bin_ordered_ind);
s_bin_ordered = bin_ordered(:,bin_ordered_ind);
Fi = (loq + sum(w_bin_ordered))/l_Meff;
Fij = zeros(n_bin_ordered_ind);
ogCOV = zeros(n_bin_ordered_ind);
gCOV = zeros(ncols,ncols);

for i = 1:n_bin_ordered_ind
    for j = i:n_bin_ordered_ind
        Fij(i,j) = (loq2 + sum(w_bin_ordered(:,i) .* s_bin_ordered(:,j)))/l_Meff;
        Fij(j,i) = Fij(i,j);
        ogCOV(i,j) = Fij(i,j) - Fi(i)*Fi(j);
        ogCOV(j,i) = ogCOV(i,j);
    end
end

gCOV(bin_ordered_ind,bin_ordered_ind) = ogCOV;

else

bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);

% Here we calculate the COV matrix only for the indices included in the
% vector 'bin_ordered_ind'.

    gCOV = zeros(ncols,ncols);
    ogCOV = cov(bin_ordered(:,bin_ordered_ind));
    gCOV(bin_ordered_ind,bin_ordered_ind) = ogCOV;

end

    % Here we recover the original values of nrows and ncols
    shift = nsymbols - 1;    
    [~,ncols] = size(nmsa);
    COV = zeros(ncols,ncols);
    for oi = 1:ncols
        i = 21*(oi-1) + 1;
        for oj = oi:ncols
            j = 21*(oj-1)+1;
            % Frobenius norm
            COV(oi,oj) = norm(gCOV(i:i+shift,j:j+shift),'fro');           
            % Spectral norm
            % subcov = gCOV(i:i+shift,j:j+shift);
            % [~,COV(oi,oj),~] =svds(subcov,1);
            % L1 norm
            % COV(oi,oj) = norm(gCOV(i:i+shift,j:j+shift),1);
            % 2 norm
            % COV(oi,oj) = norm(gCOV(i:i+shift,j:j+shift),2);
            % Simple sum
            % COV(oi,oj) = sum(sum(gCOV(i:i+shift,j:j+shift)));
            COV(oj,oi) = COV(oi,oj);
        end
    end
    
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
        pcMIP = MI_to_MIP(n_rho,ncols);
        case 'INVERSE'
        pcMIP = MI_to_MIP(-n_cov_inv,ncols);
    end

    % ZPX2 matrix
    pcZPX2 = MIP_to_ZPX2(pcMIP,ncols);
    
    pcZPX2 = real(pcZPX2);
    
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


function [W,Meff] = nmsa_to_W(nmsa,dist,threshold)
%
[nrows,ncols] = size(nmsa);

W = ones(nrows,1);

dist_threshold = dist > threshold; 

W = 1./sum(dist_threshold)';
Meff=sum(W);

fprintf('Meff = %f \n', Meff);

end


function [fnmsa] = nmsa_to_fnmsa(nmsa,W)

% Here we create a 3d matrix in which every layer has the dimensions of the
% nmsa and represents one of 20 symbols (if gaps are excluded). Then every row
% of each layer is scaled by the weight of that sequence.
symbols = unique(nmsa);
nsymbols = length(symbols);

for i = 1:nsymbols
    nmsa_3(:,:,i) = nmsa == i;
    for j = 1:nrows
        nmsa_3(j,:,i) = nmsa_3(j,:,i)*W(j);
    end
end

% nmsa_1 = sum(nmsa_3,3);

for i = 1:21
    layer = squeeze(nmsa_3(:,:,i));
    % Pseudocount.
    Fi(i,:) = (loq + sum(layer))/(loq + Meff);
end

% Here we remove the frequencies of the gaps and we scale the remaining
% frequencies.
if nsymbols == 20
sumFi = sum(Fi);
scaleFi = repmat(sumFi,21,1);
Fi = Fi./scaleFi;
% sumFi = sum(Fi);
end

% Here we create a new msa in which every symbol is replaced by the
% weigthed frequence of that symbol in that column of the msa.
fnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        fnmsa(i,j) = Fi(row(j),j);
        end
end


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
    % In this case leave the diagonal unchanged
    % MI_mat(i,i) = NaN;
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
    
%     if (ZPX2_i<0&&ZPX2_j<0)
%         pZPX2(m,n)=-pZPX2(m,n);
%     end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end



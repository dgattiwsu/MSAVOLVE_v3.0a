function [pcZPX2,pcMIP,ppvZPX2,ppvMIP,cov_inv,rho,COV,sCOV,lambda] = ...
    NMSA_to_slicePSICOV(nmsa,dist_method,threshold,psc_lambda,nsymbols,...
    norm_type,Lmat,lambda,delta,inverse_method,pc_method)
% SlicedPSICOV algorithm. It differs from NMSA_to_PSICOV because the
% covariance matrix is calculated directly from vertical 'slices' of a
% matrix, every 'slice' representing a column of msa in which the second
% dimension is 20 or 21 depending on whether gaps are included or not. The
% choice of the distance method affects the inclusion ('GAPS') or exclusion
% ('NOGAPS') of gaps in the calculation of the similarity between
% sequences. A 'threshold' is applied for the similarity between sequences.
% A pseudocount (psc_lambda) is also applied based on the number of symbols
% (nsymbols) desired (20 or 21). If psc_lambda = 0, no pseudocount is
% applied. Then, we create a 3d matrix in which every layer has the
% dimensions of the nmsa and represents one of 20 symbols (if gaps are
% excluded), and every row of each layer is scaled by the weight of that
% sequence. The covariance matrix is calculated directly by
% calculating the cross-covariance matrix between the different slices of
% the matrix and taking the norm (norm_type = |2|'fro'|) of each
% cross-covariance matrix ('2' norm is the largest singular value of the
% matrix; 'frobenius' norm is the 2-norm of the linearized matrix). Finally,
% an inverse method is used to invert the covariance matrix. Three option
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
% Possible syntax for good results with all pairs :
% [pcZPX2,~,ppvZPX2,~] = NMSA_to_slicePSICOV(nmsa,'NOGAPS',0.90,21,'fro',...
%    0.005,0.0,0.0001,'QUIC','RHO');
% Possible syntax for good results with pairs separated by at least 20 
% residues in sequence:
% [pcZPX2,~,ppvZPX2,~] = NMSA_to_slicePSICOV(nmsa,'NOGAPS',0.90,20,'fro',...
%    0.005,0.0,0.0001,'QUIC','RHO');

[nrows,ncols] = size(nmsa);

% First we replace unusual symbols possibly present in the nmsa

ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

switch dist_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(nmsa);
        dist = bin_ordered*bin_ordered';

        sdist = zeros(nrows,nrows);      
        for i = 1:nrows
            sdist(i,:) = dist(i,:)/dist(i,i);
        end

        udist = triu(sdist,1);
        ldist = tril(sdist,-1)';
        lind = udist < ldist;
        udist(lind) = 0;
        ldist(~lind) = 0;
        mdist = udist + ldist;
        dist = mdist + mdist' + eye(nrows);
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;
end

% Here we apply the weights for the similarity between sequences.

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff=round(sum(W));
end

fprintf('Meff = %d \n', Meff);

% Here we create a 3d matrix in which every layer has the dimensions of the
% nmsa and represents one of 20 symbols (if gaps are excluded). Then every
% row of each layer is scaled by the weight of that sequence.

% Pseudocount correction
loq = psc_lambda/nsymbols;
l_Meff = psc_lambda + Meff;
% If frequency of a symbol in column is 0 then the corrected frequency must
% sum up to loq/l_Meff: since we are not calculating frequency by we are
% multiplying the vectors directly we add loq/L_meff to each element of
% the vector.
psc = loq/l_Meff;

nmsa_3 = zeros(nrows,ncols,nsymbols);
for i = 1:nsymbols    
    nmsa_3(:,:,i) = nmsa == i;    
    
    % The following loop can be replaced by a shorter syntax.
    % for j = 1:nrows
    %     nmsa_3(j,:,i) = nmsa_3(j,:,i)*W(j);
    % end
    
    % Which also implements the pseudocount.
    nmsa_3(:,:,i) = nmsa_3(:,:,i).*W(:,ones(1,ncols)) + psc;    
end

% Here we calculate a cross-covariance matrix slice by slice.
COV = zeros(ncols,ncols);
c_nmsa_3 = zeros(nrows,ncols,nsymbols);

% First we center all the vectors.
for i = 1:ncols
    for j = 1:nsymbols
        c_nmsa_3(:,i,j) = nmsa_3(:,i,j) - mean(nmsa_3(:,i,j));
    end
end

% Then we calculate all the dot products and extract the norm.
for i = 1:ncols
    slice_1 = squeeze(c_nmsa_3(:,i,:));
    for j = i:ncols
        slice_2 = squeeze(c_nmsa_3(:,j,:));
        COV(i,j) = norm((slice_1'*slice_2)/nrows,norm_type);
        COV(j,i) = COV(i,j);
    end
end

% Here we test for positive definitiveness
    
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


function [binmsa] = nmsa_to_binmsa_20q(nmsa)
% Returns each sequence of length L as a vector of size 20L with 0 and 1. 
% Gaps (which would be # 25 in the original Matlab numeric representation 
% of an MSA are ignored. Thus, if at a certain position there is a gap that
% position will be converted into a vector of 20 0s.

[nseq,npos]=size(nmsa);
binmsa=zeros(nseq,20*npos);
for i=1:npos 
    for aa=1:20 
        binmsa(:,20*(i-1)+aa)=(nmsa(:,i)==aa); 
    end; 
end;

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

function [pcZPX2,lgZPX2,sparsity_value,pcMIP,lgMIP,cov_inv,rho,COV,sCOV,lambda] = ...
    NMSA_to_slPSICOV_3(nmsa,dist_method,threshold,psc_lambda,nsymbols,...
    slice_scale,col_scale,norm_type,Lmat,target_sparsity,delta_sparsity,...
    max_iter,lambda,delta,pc_method,gW_e1,gW_e2)
%
% SlicedPSICOV algorithm no.3. It differs from the default, because the
% cross-product between slices is calculated in both orders and then the
% two are summed (e.g.: crosscov = slice1'*slice2 + slice2'*slice1). The
% results are inconsistent and altogether not better than the default
% algorithm.
%
% SlicedPSICOV algorithm. It differs from NMSA_to_PSICOV because the
% covariance matrix is calculated directly from vertical 'slices' of a 3d
% matrix in which every layer has the dimensions of the nmsa and represents
% one of 20 (if gaps are excluded) or 21 symbols (if gaps are included).
% Therefore, every 'slice' represents a column of msa in which the second
% dimension is 20 or 21 depending on whether gaps are included or not. The
% choice of the distance method affects the inclusion ('GAPS') or exclusion
% ('NOGAPS') of gaps in the calculation of the similarity between
% sequences. A 'threshold' is applied for the similarity between sequences.
% Every row of each layer is scaled by the weight of that sequence. A
% pseudocount (psc_lambda) is also applied for the symbols that are not
% present in each column.  If psc_lambda = 0, no pseudocount is applied.
% 'slice_scale' determines whether the slices are divided by the number of
% observations.  This number can be the sum of the values in each slice
% ('SUM'), the total number of rows ('NROWS'), the sum of all the row
% weights ('WEIGHTS'), or 1 ('NONE'). 'col_scale' ('GAPS'|'NONE')
% determines whether at the end each slice is scaled by the percentage of
% non-gaps rows. The covariance matrix is calculated directly by
% calculating the cross-covariance matrix between the different slices of
% the matrix. 'norm_type' determines how the norm of each cross-covariance
% matrix is calculated: possible values are: norm_type = 1|2|'fro'|'L1'.
% '1' norm is the largest value that the sum of the elements of any column
% can take; '2' norm is the largest singular value of the matrix;
% 'frobenius' norm is the 2-norm of the linearized matrix; L1 is the sum of
% the absolute values of the linearized matrix. The 'QUIC' method is used
% to invert the covariance matrix. Lmat is a value between 1 and 0
% controlling the 'sparseness' of the concentration matrix (smaller values
% produce a less sparse inverse, and take more time. Larger values increase
% the sparseness, with the inverse ultimately being populated only in the
% diagonal). The default value of Lmat is 0.02. Suitable ranges are between
% 0.5 and 0.0001. Regardless of the initial value assigne to Lmat it will
% be refined in order to obtain the target sparsity requested in the input.
% Lmat can also be a regularization matrix of the same dimensions as the
% covariance matrix to be inverted. Lmat is followed by two inputs: the
% 'target_sparsity' value, and the 'delta_sparsity' value. Target+/-delta
% define the range of values within which we want the final sparsity of the
% inverse to be. Reasonable values appear to be 0.3 and 0.015. The maximum
% number of QUIC iteration allowed to reach the target sparsity is set with
% 'max_iter'. 'lambda' and 'delta' define the conditioning of the
% covariance matrix prior to the calculation of the sparse inverse. The
% recommended values are lambda = 0.0 and delta = 0.0001. Larger values of
% lambda and delta will make the calculation of the inverse faster, but
% possibly less accurate. 'pc_method' = RHO|INVERSE|BOTH determines whether
% we calculate the MIP matrix from the inverse covariance, from the RHO
% matrix, or from both and then we merge the two matrices. 'gW_e1' is the
% exponent of the gap weights applied before the ZPX2 matrix is calculated:
% 0 gives a weight of 1 to everything; 1 leaves the weight as calculated by
% the 'correct_coevmat_forgaps' function; powers between 0 and 1 increase
% the weight, powers > 1 make the weight progressively smaller. 'gW_e2' is
% the exponent of the weights computed after the ZPX2 matrix.The lgZPX2 and
% lgMIP matrices include the correction produced by the application of a
% logistic fit to the MIP data.
%
% 'NMSA_to_slPSICOV' usually requires fairly high values of the final
% sparsity (0.3-0.4) as compared to the original PSICOV, to achieve the
% optimal result.
%
% Recommended syntax if used with scale = 'SUM'; IMPORTANT: in this case
% Lmat must be small: 
% [slPSICOV,sl_lgPSICOV] = NMSA_to_slPSICOV(nmsa,'NOGAPS',0.90,1,20,...
% 'SUM','NONE','fro',0.0001,0.3,0.015,100,0.0,0.0001,'RHO',0,1);
% Recommended syntax if used with scale = 'NONE'; IMPORTANT: in this case
% Lmat must be large: 
% [slPSICOV,sl_lgPSICOV] = NMSA_to_slPSICOV(nmsa,'NOGAPS',0.90,1,20,...
% 'NONE','NONE','fro',0.4,0.3,0.015,100,0.0,0.0001,'RHO',0,1);

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

% Here we calculate the weights for the similarity between sequences.

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

nmsa_3 = zeros(nrows,ncols,nsymbols);
for i = 1:nsymbols    
    nmsa_3(:,:,i) = nmsa == i;        
end

% Here we add a total extra fraction of 1 aa equal to the pseudocount
% value to the columns of each slice that count 0; the fraction is added
% only to the rows which are not gaps. Then every row of each slice is
% scaled by the weight of that sequence, based on the similarity threshold.
% Finally we determine the mean of each column counting only the rows
% in which there are no gaps, and we subtract this mean from each column.

for i = 1:ncols
    slice = squeeze(nmsa_3(:,i,:));
    slice_zeros = slice == 0;
    row_gap_ind = nmsa(:,i) == 21;
    row_nongap_ind = nmsa(:,i) ~= 21;
    slice_zeros(row_gap_ind,:) = false;
    col_sum = sum(slice,1);
    col_nongap_ind = col_sum ~= 0;
    slice_zeros(:,col_nongap_ind) = false;
    loq = psc_lambda/sum(slice_zeros(:));
    slice(slice_zeros) = loq;
    slice = slice.*W(:,ones(1,nsymbols));
    % The following line replaces dividing by the no. of observations.
    switch slice_scale
        case 'SUM'
            slice = slice/sum(slice(:));
        case 'NROWS'
            slice = slice/nrows;
        case 'WEIGHTS'
            slice = slice/sum(W);
        case 'NONE'
    end    
    slice_mean = mean(slice(row_nongap_ind,:));
    slice_mean_mat = zeros(nrows,nsymbols);
    row_non_gap_ind = find(row_nongap_ind);
    for k = 1:length(row_non_gap_ind)
        k_ind = row_non_gap_ind(k);
        slice_mean_mat(k_ind,:) = slice_mean;
    end
    slice = slice - slice_mean_mat;
    % Here we scale each slice by the number of gaps in the original column
    % of the msa.
    sumW = sum(W);
    if strcmp(col_scale,'GAPS')
        slice = slice*(sum(row_nongap_ind.*W)/sumW);
    end    
    nmsa_3(:,i,:) = slice;
end

% Here we calculate a cross-covariance matrix slice by slice: first we
% calculate all the dot products and then extract the norm.

COV = zeros(ncols,ncols);

        if strcmp(norm_type,'L1')
            for i = 1:ncols
                slice_1 = squeeze(nmsa_3(:,i,:));
                for j = i:ncols
                    slice_2 = squeeze(nmsa_3(:,j,:));
                    cross_cov = slice_1'*slice_2 + slice_2'*slice_1;
                    COV(i,j) = sum(abs(cross_cov(:)));
                    COV(j,i) = COV(i,j);
                end
            end
        else    
            for i = 1:ncols
                slice_1 = squeeze(nmsa_3(:,i,:));
                for j = i:ncols
                    slice_2 = squeeze(nmsa_3(:,j,:));
                    COV(i,j) = norm((slice_1'*slice_2 + slice_2'*slice_1),norm_type);
                    COV(j,i) = COV(i,j);
                end
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
    
    sparsity_range = [(target_sparsity - delta_sparsity) (target_sparsity + delta_sparsity)]; 
    sparsity_value = 0;
    iter_count = 0;
    while (sparsity_value < sparsity_range(1) || sparsity_value > sparsity_range(2)) && ...
            iter_count < max_iter
        
            iter_count = iter_count +1;
            iter_count_fract = 1 - (iter_count/max_iter);
            
            [cov_inv cov_rev opt cputime iter dGap] = ...
                QUIC('default', sCOV, Lmat, 1e-6, 2, 200);
            
    sparsity_value = (nnz(cov_inv)-ncols)/(ncols*(ncols-1));    % number of nonzero

    fprintf('Sparsity = %f \n', sparsity_value);
    fprintf('Lmat = %f \n', Lmat);
                                                                % pixels
    if sparsity_value == 1
        break        
    elseif sparsity_value < sparsity_range(1)
        if sparsity_value <= 1e-6
            % Here we step linearly as the last resort, hoping we make it
            % to the target sparsity before we reach the maximum number of 
            % iterations allowed.
            target_dif = 0.1;
            target_dif = target_dif*iter_count_fract;        
            
            % Lmat = 0.9*Lmat;
            Lmat = Lmat - Lmat*target_dif;
            
        else
        target_dif = 1-(sparsity_value/target_sparsity);
        target_dif = target_dif*iter_count_fract;
        Lmat = Lmat - Lmat*target_dif;
        end
    else
        target_dif = 1-(target_sparsity/sparsity_value);
        target_dif = target_dif*iter_count_fract;        
        Lmat = Lmat + Lmat*target_dif; 
    end
    
    fprintf('Target_difference = %f \n', target_dif);
    
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
        case 'BOTH'
                pcMIP1 = MI_to_MIP(n_rho,ncols);
                pcMIP2 = MI_to_MIP(-n_cov_inv,ncols);               
                scale = scale_matrices(pcMIP1,pcMIP2);
                fprintf('best scale MIP2/MIP1 = %f \n', scale);
                pcMIP = pcMIP2 + scale*pcMIP1;
    end

    lgMIP = mat_to_lg_mat(pcMIP);
    
% 1st Correction for gaps -------------------------------------------------
if gW_e1 == 0
    fprintf('1st Gap weights section skipped. \n');
    else
    fprintf('1st Gap weights section started. \n');
    tic   
    gapW = correct_coevmat_forgaps(nmsa);
    gW = gapW.^gW_e1;

    pcMIP = pcMIP.*gW;
    lgMIP = lgMIP.*gW;

    gapW_time = toc;
    fprintf('1st Gap weights section completed in %8.4f minutes. \n', gapW_time/60);
end        
% ---------------------------------------------------------------------

    % ZPX2 matrix
    pcZPX2 = MIP_to_ZPX2(pcMIP,ncols);    
    pcZPX2 = real(pcZPX2);        
    lgZPX2 = mat_to_lg_mat(pcZPX2);
        
% 2nd correction for gaps -------------------------------------------------
if gW_e2 == 0
    fprintf('2nd Gap weighting skipped. \n');
else
    tic
    if gW_e1 == 0
        gapW = correct_coevmat_forgaps(nmsa);
    end
    
    gW2 = gapW.^gW_e2;
    
    pcZPX2_orig = pcZPX2; 
    pcZPX2_orig_2 = pcZPX2_orig - min(pcZPX2_orig(:));
    pcZPX2 = pcZPX2_orig_2.*gW2;

    lgZPX2_orig = lgZPX2; 
    lgZPX2_orig_2 = lgZPX2_orig - min(lgZPX2_orig(:));
    lgZPX2 = lgZPX2_orig_2.*gW2;
        
    gW2_time = toc;
    fprintf('2nd Gap weighting completed in %8.4f seconds. \n', gW2_time);
end        
% -------------------------------------------------------------------------
    
        
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


function [ gapW ] = correct_coevmat_forgaps( nmsa )
% This function finds the positions in the msa where there are no gaps in both
% columns considered in a MI calculation and calculates a matrix of weights
% to scale the coevolution matrix.

[nseq,npos] = size(nmsa);
gapW = zeros(npos,npos);

for i = 1:npos
    for j = i:npos
    gap = find((nmsa(:,i)~=21) & (nmsa(:,j)~=21));
    gapW(i,j) = numel(gap);
    gapW(j,i) = gapW(i,j);
    end
end

gapW=gapW/nseq;
       
end


function [lg_mat] = mat_to_lg_mat(mat)

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
% lg_data = z_data;

% Logistic fit PSICOV style
data_mean = nanmean(data);
data_std = nanstd(data);
z_data = (data - data_mean)/data_std; 
lg_data = 0.904 ./ (1.0 + 16.61 * exp(-0.8105 * z_data));

lg_mat = zeros(cols);
[i,j] = ind2sub(size(mat),all_ind);

for n = 1:ndata
    lg_mat(i(n),j(n)) = lg_data(n);
end

end

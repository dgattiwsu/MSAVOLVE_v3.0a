function [pcZPX2,pcMIP,ppvZPX2,ppvMIP,cov_inv,...
    rho,COV,sCOV,lambda,MI,W,cov_rev,opt,cputime,iter,dGap] = ...
    NMSA_to_swPSIMI_2(nmsa,dist_method,...
    threshold,psc_lambda,nsymbols,Lmat,lambda,delta,pc_method)
% Matlab experimental implementations of the PSIMI algorithm with weights to 
% account for sequence similarity and inclusion of pseudocounts. First, it
% calculates a MI matrix from the columns of an MSA. The similarity between 
% sequences is derived from a distance matrix, which is calculated including 
% (dist_method = 'GAPS') or excluding (dist_method = 'NOGAPS') any gaps. 
% Similarity is taken into account on the basis of the 'threshold' value
% (for example 0.8 for 80% similarity). The psudocount is controlled by the 
% parameter 'psc_lambda'. If threshold = 1, all sequences receive a weight  
% of 1 (recommended), and if psc_lambda is set to zero no pseudocount is  
% applied. 'nsymbols' (20 or 21) is a divider that scales psc_lambda. The 
% QUIC sparse inverse method is used to invert the MI matrix. Lmat values 
% between 1 and 0 control the 'sparseness' of the concentration matrix: smaller 
% values produce a less sparse inverse, and take more time. Larger values 
% increase the sparseness, with the inverse ultimately being populated only 
% in the diagonal. The default value of Lmat is 0.02. Suitable ranges are 
% between 0.025 and 0.005. Lmat can also be a regularization matrix of the 
% same dimensions as the MI matrix to be inverted. 'lambda' and 'delta' 
% define the conditioning of the MI matrix prior to the calculation 
% of the sparse inverse. The recommended values are lambda = 0.0 and 
% delta = 0.0001. Larger values of lambda and delta will make the calculation 
% of the inverse faster, but possibly less accurate. 'pc_method'
% determines whether we calculate the MIP matrix from the sparse inverse
% MI or from the RHO matrix. The ppvZPX2 and ppvMIP
% matrices include the correction produced by the application of a logistic
% fit. Possible usage: [pcZPX2,pcMIP,ppvZPX2,ppvMIP] =
% NMSA_to_swPSIMI(nmsa,'GAPS',1.0,0.5.21,0.02,0.0,0.0001,'RHO');
%
% It differs from NMSA_to_swPSIMI because it produces a coevolution matrix
% in which the diagonal is not NaN. This is useful for a spectral analysis
% of the coevolution matrix.


% Replace unusual symbols
ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

[nrows,ncols] = size(nmsa);

% Here we calculate the distance matrix and weights for the similarity
% between sequences.
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
        mdist = mdist + mdist' + eye(nrows);
       
    % The merged distance matrix becomes the final distance matrix.
    
        dist = mdist;
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;

end

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff=round(sum(W));
end

fprintf('Meff = %d \n', Meff);

    % MI matrix
    [MI] = NMSA_to_MI(nmsa,W,psc_lambda,nsymbols);
    
    % For consistency with the PSICOV code we rename here the MI matrix.
    COV = MI;
    
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
        % weight to pairs known to be important.
        % Lmat = Lmat*NMSA_to_MI(nmsa) + imatrix;
        % Lmat = Lmat*NMSA_to_MI(nmsa);        
        % Lmat = Lmat*COV + imatrix;
    else
        Lmat = 0.02;
    end
    
    % QUIC sparse inverse method.
    [cov_inv,cov_rev,opt,cputime,iter,dGap] = ...
            QUIC('default', sCOV, Lmat, 1e-6, 2, 200);
    
    nonzeros = (nnz(cov_inv)-ncols)/(ncols*(ncols-1));
    fprintf('Sparsity = %f \n', nonzeros);

    % Partial correlation matrix.    
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
    
%     for i = 1:ncols
%         n_rho(i,i) = NaN;
%         n_cov_inv(i,i) = NaN;
%     end
    
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


function [ MI_mat ] = NMSA_to_MI(nmsa,W,psc_lambda,nsymbols)

[~,ncols] = size(nmsa);
MI_mat = zeros(ncols,ncols);

for i = 1:ncols
    for j = i:ncols
    MI_mat(i,j) = MutualInformation(nmsa(:,i),nmsa(:,j),W,psc_lambda,nsymbols);
    MI_mat(j,i) = MI_mat(i,j);    
    end
    % MI_mat(i,i)=NaN;
end

end


function I = MutualInformation(X,Y,W,psc_lambda,nsymbols)

% MutualInformation: returns mutual information (in bits) of the 'X' and
% 'Y'. Modified from Will Dwinnell's original.

    I = Entropy(X,W,psc_lambda,nsymbols,1) + ...
        Entropy(Y,W,psc_lambda,nsymbols,1) - ...
        JointEntropy([X Y],W,psc_lambda,nsymbols,2);

end


function H = Entropy(X,W,psc_lambda,nsymbols,power)

% Entropy: Returns entropy (in bits) of each column of 'X' weigthed by the
% similarity between each sequence. Modified from Will Dwinnel's original.

% Establish size of data
[~,m] = size(X);

% Housekeeping
H = zeros(1,m);

for Column = 1:m,
    % Assemble observed alphabet
    Alphabet = unique(X(:,Column));
	
    % Housekeeping
    local_symbols = length(Alphabet);

    Count = zeros(local_symbols,1);
    loq = psc_lambda/(nsymbols^power);

%--Weighted frequencies----------------------------------------------------    

%     % Calculate sample frequencies
%     for symbol = 1:local_symbols
%         Count_vector = (X(:,Column) == Alphabet(symbol)).*W;
%         Count(symbol) = loq + sum(Count_vector);
%     end
% 	
%     % Calculate sample class probabilities: hard to see which of the two
%     % lines below is more correct theoretically. Using 'length(X)' seems to
%     % require more uniform values of the QUIC subroutine 'lambda'.
%     % l_Meff = lambda + sum(Count);
%     l_Meff = psc_lambda + length(X);    
%     P = Count / l_Meff;
% 
%     % Calculate entropy in bits
%     % Note: floating point underflow is not an issue since we are
%     % dealing only with the observed alphabet.
%     
%     H(Column) = -sum(P .* log2(P));
    
%--Weighted MI-------------------------------------------------------------

    W_vector = zeros(local_symbols,1);
    	
    % Calculate sample frequencies
    for symbol = 1:local_symbols
        Count_vector = (X(:,Column) == Alphabet(symbol));
        W_vector(symbol) = sum(W(Count_vector))/sum(Count_vector);
        Count(symbol) = loq + sum(Count_vector);
    end
	
    % Calculate sample class probabilities
    l_Meff = psc_lambda + sum(Count);
    % P = Count .* W_vector / l_Meff;
    P = Count / l_Meff;

    % Calculate entropy in bits
    % Note: floating point underflow is not an issue since we are
    % dealing only with the observed alphabet.
    
    % H(Column) = -sum(P .* log2(P));
    H(Column) = -sum(P .* W_vector .* log2(P));

end

end


function H = JointEntropy(X,W,psc_lambda,nsymbols,power)
% JointEntropy: Returns joint entropy (in bits) of each column of 'X'.
% by Will Dwinnell
%
% H = JointEntropy(X)
%
% H = calculated joint entropy (in bits)
% X = data to be analyzed
%
% Last modified: Aug-29-2006

% Sort to get identical records together
[X,ind] = sortrows(X);

% Find elemental differences from predecessors
DeltaRow = (X(2:end,:) ~= X(1:end-1,:));

% Summarize by record
Delta = [1; any(DeltaRow')'];

% Generate vector symbol indices
cumDelta = cumsum(Delta);
VectorX = cumDelta(ind);

% Calculate entropy the usual way on the vector symbols
H = Entropy(VectorX,W,psc_lambda,nsymbols,power);
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

% Here we calculate the MCA and W matrices
for m = 1:ncols
    for n = m:ncols
    MCA_mat(m,n)=(mean_row(m)*mean_row(n))/mean_mat;
    MCA_mat(n,m) = MCA_mat(m,n);    
    end
% MCA_mat(m,m) = NaN;
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
% Gloor, but is included in the ZRES algorithm by Chen.
    
    if (ZPX2_i<0&&ZPX2_j<0)
        pZPX2(m,n)=-pZPX2(m,n);
    end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    % pZPX2(m,m)=NaN;

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

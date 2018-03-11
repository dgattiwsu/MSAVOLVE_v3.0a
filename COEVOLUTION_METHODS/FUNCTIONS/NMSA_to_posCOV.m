function [COV,sCOV] = ...
    NMSA_to_posCOV(nmsa,threshold,psc_lambda,recov_method,lambda,delta,nsymbols)
% This function calculates a covariance matrix for the columns of an MSA.
% The MSA is first converted to its binary 'long' representation in which
% every column is extended into 21 different columns, each one representing
% a different aa. A threshold is applied for the similarity between
% sequences. A pseudocount (psc_lambda) is also applied. Then, the
% covariance matrix is calculated and is shrunk to the equivalent
% covariance matrix for the 'short' standard representation. The value of
% nsymbols determines whether we include or not gaps in the calculation of
% the 'shrunk' covariance matrix (nsymbols = 21 includes gaps; nsymbols =
% 20 excludes gaps). sCOV is the closest invertible covariance matrix
% obtained by 'shrinking' the original covariance matrix. Lambda and delta
% are the shrinking parameters. Usage: [posCOV] =
% NMSA_to_posCOV(nmsa,0.9,0.5,'FRO',0.0,0.001,20);

bin_ordered = nmsa_to_binmsa_21q(nmsa);
[nrows,ncols] = size(bin_ordered);

% Here we calculate the distance matrix on the transpose and apply the
% weights for the similarity between sequences.
bin_ordered_sum = sum(bin_ordered);
bin_ordered_ind = find(bin_ordered_sum);
n_bin_ordered_ind = length(bin_ordered_ind);

dist = bin_ordered*bin_ordered';

if threshold < 1.0
    [W,Meff] = nmsa_to_W(nmsa,dist,threshold);
else
    W = ones(nrows,1);
    Meff = nrows;
end

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

    % Here we recover the original values of nrows and ncols
    shift = nsymbols - 1;    
    [~,ncols] = size(nmsa);
    COV = zeros(ncols,ncols);
    for oi = 1:ncols
        i = 21*(oi-1) + 1;
        for oj = oi:ncols
            j = 21*(oj-1)+1;
            switch recov_method
                case 'FRO'
                % Frobenius norm
                COV(oi,oj) = norm(gCOV(i:i+shift,j:j+shift),'fro');
                case 'SPECTRAL'
                % Spectral norm
                subcov = gCOV(i:i+shift,j:j+shift);
                [~,COV(oi,oj),~] =svds(subcov,1);
                case 'L1'
                % L1 norm
                COV(oi,oj) = norm(gCOV(i:i+shift,j:j+shift),1);
                case '2'
                % 2 norm
                COV(oi,oj) = norm(gCOV(i:i+shift,j:j+shift),2);
                case 'SUM'
                % Simple sum
                COV(oi,oj) = sum(sum(gCOV(i:i+shift,j:j+shift)));
            end
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

dist_threshold = dist > threshold*ncols; 

W = 1./sum(dist_threshold)';
Meff=sum(W);

fprintf('Meff = %f \n', Meff);

end

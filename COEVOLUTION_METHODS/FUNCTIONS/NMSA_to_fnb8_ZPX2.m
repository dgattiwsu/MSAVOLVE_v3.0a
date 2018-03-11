function [ZPX2,lgZPX2,ZPX2_orig,MI,sW,gW,gW1,gW2] = ...
    NMSA_to_fnb8_ZPX2(nmsa,dist_method,threshold,psc_lambda,...
    null_symbol,gW_e1,gW_e2,gW_e3)
% This variation of the nb algorithm implements a 'fastMI' calculation that
% is ~15 times faster than standard MI, weights for the similarity between
% sequences, a pseudocount, and weights for the amount of gaps in each
% column. The similarity between sequences can be determined including gaps
% ('dist_method'= 'GAPS') or excluding gaps ('dist_method'= 'NOGAPS'), on
% the basis of a 'threshold' value. In general, including gaps gives better
% results. A pseudocount is also applied for the symbols that are not
% present in each column by adding a pseudocount value (psc_lambda) to the
% corrected count of equivalent sequences in the alignment. If psc_lambda =
% 0, no pseudocount is applied. There is a first round of gap weighting
% before the calculation of the ZPX2 matrix. 'gW_e1' is the exponent of the
% weights for this round: 0 gives a weight of 1 to everything. 1 leaves the
% weight as calculated by the 'correct_coevmat_forgaps' function; higher
% powers make the weight progressively smaller. A second round of gap
% weighting is carried out after the calculation of the ZPX2 matrix:
% 'gW_e2' is the exponent of the weights for this round. 'null_symbol' is a
% numeric value (between 1 and 21) that can be used to remove a symbol from the count (for
% example, use 21 to remove the gap count). If set to any value higher than
% 21, every symbol is counted.
%
% Recommended syntax:
%
% 1. We do a double gap correction before the MIP matrix and after the 
% ZPX2 matrix. This option gives the best match to the contact map when
% 'all' residues are considered, and excellent match when considering
% residues separated by 20 or more intervening positions.
% [nbZPX2] = NMSA_to_fnb_ZPX2(nmsa,'GAPS',0.9,1,21,2,0,3);
%
% 2. We do a gap correction after the calculation of the MIP matrix and/or ZPX2
% matrix. This option gives excellent match to the contact map when 'all'
% residues are considered, and slightly better match when considering
% residues separated by 20 or more intervening positions. This syntax is
% less prone to give artifacts when there are many gaps in the alignment.
% [nbZPX2] = NMSA_to_fnb_ZPX2(nmsa,'GAPS',0.9,1,22, 0,0,3);
% or:
% [nbZPX2] = NMSA_to_fnb_ZPX2(nmsa,'GAPS',0.9,1,22,0,1,3);

[nrows,ncols] = size(nmsa);

% Here we replace unusual symbols possibly present in the nmsa

ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

% Here we determine the distance matrix: more distant sequences are 
% represented by smaller numbers. This is only to determine which sequence
% will be at the top of the msa.

bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';

triu_dist = triu(dist,1);
[~,max_ind] = max(triu_dist(:));
[firstrow,~] = ind2sub([nrows,nrows],max_ind);

% Here we select the top sequence. 

nmsa1 = nmsa(firstrow:nrows,:);
nmsa2 = nmsa(1:firstrow-1,:);
nmsa = [nmsa1;nmsa2];

%--------------------------------------------------------------------------
% Here we reorder the msa.
ordered = reorder_nmsa_1(nmsa,nrows);
%--------------------------------------------------------------------------

% Here we calculate the weights for the similarity between sequences. We
% need to determine a new distance matrix for the reordered msa.

switch dist_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(ordered);
        dist = bin_ordered*bin_ordered';
        nsymbols = 20;

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
        bin_ordered = nmsa_to_binmsa_21q(ordered);
        dist = (bin_ordered*bin_ordered')/ncols;
        nsymbols = 21;
end

global W loq loq2 l_Meff const symbol
symbol = null_symbol+1;
const = log2(l_Meff);

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff = round(sum(W));
    % A higher pseudocount appears to give better results with the nb 
    % algorithm.
    % psc_lambda = psc_lambda*Meff/nrows;  
    l_Meff = psc_lambda + Meff;
    loq = psc_lambda/(l_Meff*nsymbols);
    % loq2 = psc_lambda/(l_Meff*nsymbols^2);
    loq2 = psc_lambda/(l_Meff*nsymbols*2);
    sW = W;
end

fprintf('Meff = %d \n', Meff);

% Here we convert the msa from matlab numeric format to binary
% differential.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end

binsum = sum(bin_ordered(:));
fprintf('Binary Sum = %f \n', binsum);    
    
% Here we set to zero all the indices of the reordered msa that are 0 in 
% the binary msa. 
z_ordered = ordered;
z_ordered(~bin_ordered) = 0;

%--------------------------------------------------------------------------
% MI matrix calculation with fastMI: 1 is added to every position of the
% msa because fastMI can't take 0's.
    z_ordered = z_ordered+1;                 
    % MI = NMSA_to_fastMI(z_ordered);
    MI = NMSA_to_2D_MI(z_ordered,ncols);        
    
%--------------------------------------------------------------------------


% General gap correction --------------------------------------------------
    gW = correct_coevmat_forgaps(nmsa);

% 1st correction for gaps -------------------------------------------------
if gW_e1 == 0
    fprintf('1st Gap weighting skipped. \n');
else
    tic   
    gW1 = gW.^gW_e1;

    MI = MI.*gW1;
    
    gW1_time = toc;
    fprintf('1st Gap weighting completed in %8.4f seconds. \n', gW1_time);
end        
% -------------------------------------------------------------------------
        
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    
% 2nd correction for gaps -------------------------------------------------
if gW_e2 == 0
    fprintf('2nd Gap weighting skipped. \n');
    else
    % fprintf('2nd Gap weights section started. \n');
    tic   
    gW2 = gW.^gW_e2;

    MIP = MIP.*gW2;

    gW2_time = toc;
    fprintf('1st Gap weighting completed in %8.4f seconds. \n', gW2_time);
end        
        
    % ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    lgZPX2 = mat_to_lg_mat(ZPX2);
    
% 3rd correction for gaps -------------------------------------------------
if gW_e3 == 0
    fprintf('3rd Gap weighting skipped. \n');
else
    tic    
    gW3 = gW.^gW_e3;
    
    ZPX2_orig = ZPX2; 
    ZPX2_orig_2 = ZPX2_orig - min(ZPX2_orig(:));
    ZPX2 = ZPX2_orig_2.*gW3;
    lgZPX2_orig = lgZPX2; 
    lgZPX2_orig_2 = lgZPX2_orig - min(lgZPX2_orig(:));
    lgZPX2 = lgZPX2_orig_2.*gW3;
        
    gW3_time = toc;
    fprintf('3rd Gap weighting completed in %8.4f seconds. \n', gW3_time);
end        
% -------------------------------------------------------------------------
    
end


function [binmsa] = nmsa_to_binmsa_20q(nmsa)
% Returns each sequence of length L as a vector of size 20L with 0 and 1. 
% Gaps (which would be # 25 in the original Matlab numeric representation 
% of an MSA) are ignored.

[nseq,npos]=size(nmsa);
binmsa=zeros(nseq,20*npos);
for i=1:npos 
    for aa=1:20 
        binmsa(:,20*(i-1)+aa)=(nmsa(:,i)==aa); 
    end; 
end;

end

function [binmsa] = nmsa_to_binmsa_21q(nmsa)
%--------------------------------------------------------------------------    
% Returns each sequence of length L as a vector of size 21L with 0 and 1. 
% Number 21 represents gaps (which would be # 25 in the original
% Matlab numeric representation of an MSA.

[nseq,npos]=size(nmsa);
ind25 = nmsa == 25;
nmsa(ind25) = 21;
binmsa=zeros(nseq,21*npos);
for i=1:npos 
    for aa=1:21 
        binmsa(:,21*(i-1)+aa)=(nmsa(:,i)==aa); 
    end 
end

end


function [ordered] = reorder_nmsa_1(nmsa,nrows)
%--------------------------------------------------------------------------    
    for k = 1:nrows-1
    nmsa_dif = zeros(nrows,1);
        for i = k:nrows
        nmsa_dif(i) = sum(nmsa(k,:)~=nmsa(i,:));
        end
    [~, sort_msa_ind] = sort(nmsa_dif,'ascend');
    nmsa = nmsa(sort_msa_ind,:);
    end
    ordered = nmsa;
    
end


function [ MI_mat ] = NMSA_to_fastMI( nmsa )
%--------------------------------------------------------------------------
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (nmsa). The function uses  the "fastMI" algorithm. 

[~,ncols] = size(nmsa);
MI_mat = zeros(ncols,ncols);

for i = 1:ncols
    for j = i:ncols
    MI_mat(i,j) = fastMI(nmsa(:,i),nmsa(:,j));
    MI_mat(j,i) = MI_mat(i,j);    
    end
    MI_mat(i,i)=NaN;
end

end


function MI = fastMI(X,Y)
%--------------------------------------------------------------------------
% Mutual information between two column vectors, X and Y, having integer
% values. W is the weight vector accounting for the similarity between
% sequences. We use the scaled value l_Meff for the number of observations
% instead of the total number N. A pseudocount is added to the count of all
% observed pairs.

global W loq2 l_Meff

N = l_Meff;
const = log2(N);
h = accumarray([X Y], W); 
xy_prod = sum(h,2)*sum(h,1);
xy_ind = h~=0;
h(xy_ind) = h(xy_ind) + loq2;  % pseudocount added
MI = const + (sum(h(xy_ind) .* log2(h(xy_ind)./xy_prod(xy_ind))))/N;

end


function [MI] = NMSA_to_2D_MI(nmsa,ncols)
%--------------------------------------------------------------------------
% MI calculation.

global W loq loq2 l_Meff const symbol

MI = NaN(ncols,ncols);
% MI_3d = NaN(ncols,ncols,ncols);


fprintf('Hi section started. \n');
time1 = tic;
Hi = NaN(ncols,1);
for i = 1:ncols
        i_col = nmsa(:,i);
        i_ind = any(i_col == symbol,2);
        i_col = i_col(~i_ind);
        i_W = W(~i_ind);
        
        h = accumarray(i_col, i_W);
                
        i_marg = nonzeros(h) + loq;     % pseudocount
        
        % Entropy of i
        Hi(i) = const - sum(i_marg.*log2(i_marg))/l_Meff;
end

Hi_time = toc(time1);
fprintf('Hi section completed in %8.4f minutes. \n', Hi_time/60);


fprintf('Hij section started. \n');
time2 = tic;
Hij = NaN(ncols,ncols);
for i = 1:ncols    
    for j = i+1:ncols
        ij_col = nmsa(:,[i j]);
        ij_ind = any(ij_col == symbol,2);
        ij_col = ij_col(~ij_ind,:);
        ij_W = W(~ij_ind);
        
        h = accumarray(ij_col, ij_W);        
                
        ij_marg = nonzeros(h) + loq2;   % pseudocount

        % Joint entropy of i and j
        Hij(i,j) = const - sum(ij_marg.*log2(ij_marg))/l_Meff;         
        Hij(j,i) = Hij(i,j);
        
        % Mutual Information of i and j
        MI(i,j) = Hi(i) + Hi(j) - Hij(i,j);
        MI(j,i) = MI(i,j);
        
    end
end

Hij_time = toc(time2);
fprintf('Hij section completed in %9.4f minutes. \n', Hij_time/60);
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
% would give the wrong MI. Change of sign is not in the original ZPX2
% algorithm by Gloor, but is included in the ZRES algorithm by Chen.
    
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

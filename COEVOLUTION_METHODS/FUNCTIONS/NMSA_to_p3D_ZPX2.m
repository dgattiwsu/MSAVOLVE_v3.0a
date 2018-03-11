function [ZPX2,md3_ZPX2,lgZPX2,md3_lgZPX2,ZPX2_orig,md3_ZPX2_orig,MI,wMI,...
    md3_MI,wmd3_MI,sW,gW1,gW2,gW3] = ...
    NMSA_to_p3D_ZPX2(nmsa,dist_method,threshold,psc_lambda,...
    null_symbol,gW_e1,gW_e2,gW_e3,nprocs)
% This implementation of the multidimensional MI algorithm uses a 'fastMI'
% calculation that is ~15 times faster than standard MI, weights for the
% similarity between sequences, a pseudocount, and weights for the amount
% of gaps in each column. The function calculates the 'standard' ZPX2
% matrix, and the 'direct coupling' md3_ZPX2, in which the the 'interaction
% information' between each i,j,k, column is subtracted from MI(i;j). The
% final ZPX2/md3_ZPX2 maps can be further processed using a logistic fit
% PSICOV style to produce the lgZPX2/md3_lgZPX2 matrices. The similarity
% between sequences can be determined including gaps ('dist_method'=
% 'GAPS') or excluding gaps ('dist_method'= 'NOGAPS'), on the basis of a
% 'threshold' value. In general, including gaps gives better results. A
% pseudocount is also applied for the symbols that are not present in each
% column by adding a pseudocount value (psc_lambda) to the scaled count of
% equivalent sequences in the alignment. 'gW_e1' is the exponent of the
% weights on the 1st round of gap weighting before the MIP matrix is
% calculated: 0 gives a weight of 1 to everything; 1 leaves the weight as
% calculated by the 'correct_coevmat_forgaps' function; powers between 0
% and 1 increase the weight, powers > 1 make the weight progressively
% smaller. Additional 2 rounds of gap weighting are possible after the MIP
% matrix is calculated ('gW_e2'), and after the ZPX2 matrix is calculated
% ('gW_e3'). 'null_symbol' is a numeric value (between 1 and 21) that can
% be used to remove a symbol from the count (for example, use 21 to remove
% the gap count). If set to any value higher than 21, every symbol is
% counted.
%
% Recommended syntax:
%
% 1. We include gaps in the calculation of the distance matrix, but we do
% not count the gaps in the marginal frequencies; finally we do gap
% correction after the MIP matrix and/or the ZPX2 matrix. This syntax is
% more prone to give artifacts when there are large blocks of gaps in the
% alignment.
% [ZPX2,md3_ZPX2,md3_lgZPX2] = NMSA_to_3D_ZPX2(nmsa,'GAPS',0.9,1,21,0,0,3);
%
% 2. We count the gaps both in the calculation of the distance matrix and
% of the marginal frequencies, and we do gap correction after the MIP
% matrix and/or the ZPX2 matrix. This syntax is less prone to give
% artifacts when there are large blocks of gaps in the alignment.
% [ZPX2,md3_ZPX2,md3_lgZPX2] = NMSA_to_3D_ZPX2(nmsa,'GAPS',0.9,1,22,0,1,3);
%
%--------------------------------------------------------------------------

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

% Here we calculate the weights for the similarity between sequences. 

switch dist_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(nmsa);
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
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;
        nsymbols = 21;
end


global W loq loq2 loq3 l_Meff const symbol
symbol = null_symbol;

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff = round(sum(W));
    % The pseudocount is scaled to the number of equivalent sequences.
    psc_lambda = psc_lambda*Meff/nrows;
end

    l_Meff = psc_lambda + Meff;
    loq = psc_lambda/(l_Meff*nsymbols);
    loq2 = psc_lambda/(l_Meff*nsymbols*2);
    loq3 = psc_lambda/(l_Meff*nsymbols*3);
    % loq2 = psc_lambda/(l_Meff*nsymbols^2);
    % loq3 = psc_lambda/(l_Meff*nsymbols^3);    
    sW = W;    

fprintf('Meff = %d \n', Meff);

const = log2(l_Meff);

    if matlabpool('size') > 0
        matlabpool close
    end

    if nprocs > 1
        matlabpool(nprocs)
    end

[MI,md3_MI] = NMSA_to_md3_MI(nmsa,ncols);

    if nprocs > 1
        matlabpool close
    end

% General gap correction --------------------------------------------------
    gapW = correct_coevmat_forgaps(nmsa);

% 1st Correction for gaps -------------------------------------------------
if gW_e1 == 0
    wMI = MI;
    wmd3_MI = md3_MI;
    fprintf('1st Gap weights section skipped. \n');
    else
    % fprintf('1st Gap weights section started. \n');
    tic   
    gW1 = gapW.^gW_e1;

    wMI = MI.*gW1;
    wmd3_MI = md3_MI.*gW1;

    gapW1_time = toc;
    fprintf('1st Gap weights section completed in %8.4f minutes. \n', gapW1_time/60);
end        
% ---------------------------------------------------------------------

% Before any corrections we set the diagonal of the coevolution matrices to 
% all NaN's.

for i = 1:ncols
    wMI(i,i) = NaN;
    wmd3_MI(i,i) = NaN;
end

% MIP
MIP = MI_to_MIP(wMI,ncols);
md3_MIP = MI_to_MIP(wmd3_MI,ncols);

% 2nd Correction for gaps -------------------------------------------------
if gW_e2 == 0
    fprintf('2nd Gap weights section skipped. \n');
    else
    % fprintf('2nd Gap weights section started. \n');
    tic   
    gW2 = gapW.^gW_e2;

    MIP = MIP.*gW2;
    md3_MIP = md3_MIP.*gW2;

    gapW2_time = toc;
    fprintf('1st Gap weights section completed in %8.4f minutes. \n', gapW2_time/60);
end        
% ---------------------------------------------------------------------

% ZPX2
ZPX2 = MIP_to_ZPX2(MIP,ncols);
md3_ZPX2 = MIP_to_ZPX2(md3_MIP,ncols);
% Logistic correction PSICOV style.
lgZPX2 = mat_to_lg_mat(ZPX2);
md3_lgZPX2 = mat_to_lg_mat(md3_ZPX2);

% 3rd correction for gaps -------------------------------------------------
if gW_e3 == 0
    fprintf('3rd Gap weighting skipped. \n');
else
    tic
    gW3 = gapW.^gW_e3;
    
    ZPX2_orig = ZPX2; 
    ZPX2_orig_2 = ZPX2_orig - min(ZPX2_orig(:));
    ZPX2 = ZPX2_orig_2.*gW3;

    lgZPX2_orig = lgZPX2; 
    lgZPX2_orig_2 = lgZPX2_orig - min(lgZPX2_orig(:));
    lgZPX2 = lgZPX2_orig_2.*gW3;
        
    md3_ZPX2_orig = md3_ZPX2; 
    md3_ZPX2_orig_2 = md3_ZPX2_orig - min(md3_ZPX2_orig(:));
    md3_ZPX2 = md3_ZPX2_orig_2.*gW3;

    md3_lgZPX2_orig = md3_lgZPX2; 
    md3_lgZPX2_orig_2 = md3_lgZPX2_orig - min(md3_lgZPX2_orig(:));
    md3_lgZPX2 = md3_lgZPX2_orig_2.*gW3;
        
    gW3_time = toc;
    fprintf('3rd Gap weighting completed in %8.4f seconds. \n', gW3_time);
end        
% -------------------------------------------------------------------------

end



function [binmsa] = nmsa_to_binmsa_20q(nmsa)
%--------------------------------------------------------------------------    
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


function [MI,mdMI] = NMSA_to_md3_MI(nmsa,ncols)
%--------------------------------------------------------------------------
% MI, MI_3d calculation.

global W loq loq2 loq3 l_Meff const symbol

MI = NaN(ncols,ncols);
% MI_3d = NaN(ncols,ncols,ncols);


fprintf('Hi section started. \n');
tic
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

Hi_time = toc;
fprintf('Hi section completed in %8.4f minutes. \n', Hi_time/60);


fprintf('Hij section started. \n');
tic
Hij = NaN(ncols,ncols);
for i = 1:ncols    
    for j = i:ncols
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

Hij_time = toc;
fprintf('Hij section completed in %9.4f minutes. \n', Hij_time/60);


fprintf('Hijk section started. \n');

if ncols < 200
    total_count = 0;
    for i = 1:ncols
        for j = i:ncols        
            for k = j:ncols
             total_count = total_count + 1;         % Exact count.  
            end        
        end               
    end
else
    total_count = (ncols^3 + ncols^2 + ncols)/6;    % Approximate count.
end
    
count = 0;
time = 0;

Hijk = NaN(ncols,ncols,ncols);
psymbol = symbol;
pW = W;
ploq3 = loq3;
pconst = const;
pl_Meff = l_Meff;

for i = 1:ncols
    tic
    
    parfor j = i:ncols
        Hij_col = NaN(ncols,1);
        pnmsa = nmsa;
        
        for k = j:ncols
            count = count +1;
            
            ijk_col = pnmsa(:,[i j k]);
            ijk_ind = any(ijk_col == psymbol,2);
            ijk_W = pW;
            h = accumarray(ijk_col(~ijk_ind,:), ijk_W(~ijk_ind));
                        
            ijk_marg = nonzeros(h) + ploq3;      % pseudocount

            % Joint Entropy of i,j and k
            Hij_col(k) = pconst - sum(ijk_marg.*log2(ijk_marg))/pl_Meff;
        end
        
        Hijk(i,j,:) = Hij_col;
    end
    
    cycle_time = toc;
    time = time + cycle_time;  
    progress = count/total_count;
    count_left = total_count - count;
    fprintf('Percent completion = %5.3f \n', progress);
    time_left = count_left*(time/count)/60;
    fprintf('Time elapsed: %6.1f minutes \n',time/60);            
    fprintf('Est. time to Hijk completion: %6.1f minutes \n',time_left);
    
end

for i = 1:ncols
    for j = i:ncols
        for k = j:ncols
            % We store also all the permutations. 
            Hijk(i,k,j) = Hijk(i,j,k);
            Hijk(j,i,k) = Hijk(i,j,k);
            Hijk(j,k,i) = Hijk(i,j,k);
            Hijk(k,i,j) = Hijk(i,j,k);
            Hijk(k,j,i) = Hijk(i,j,k);            
        end        
    end                
end

tic
mdMI = NaN(ncols,ncols);
all_ind = (1:ncols);
for i = 1:ncols
    for j = i:ncols
        k_ind = setdiff(all_ind,[i j]);
        
            % The following is the MI minus the interaction information.
            % mdMI(i,j,k) = Hij(i,k) - Hijk(i,k,j) - Hi(k) + Hij(k,j);
            % Since we are averaging over the 3rd dimension, only i and j
            % are permuted.
                        
            mdMI(i,j) = nanmean(Hij(j,k_ind)' + Hij(i,k_ind)' ...
                        - Hi(k_ind) - squeeze(Hijk(i,j,k_ind)));
            mdMI(j,i) = mdMI(i,j);
                    
    end
end

mdMI_time = toc;
fprintf('Hijk/mdMI section completed in %9.4f minutes. \n', (time+mdMI_time)/60);

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
% Gloor, but is included in the ZRES algorithm by Chen.
    
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



function [ZPX2,MIP,MI,gW] = NMSA_to_fnb5_ZPX2(nmsa,gapW_exp)
% This is the most effective nb algorithm, implementing an improved
% sequence resorting method and a new fast MI calculation. 

% First we determine the distance matrix: more distant sequences are 
% represented by smaller numbers.

[nrows,ncols] = size(nmsa);

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

% Here we convert the msa from matlab numeric format to binary
% differential.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end

binsum = sum(bin_ordered(:));
fprintf('Binary Sum = %f \n', binsum);    
    
% Here we set to zero all the indices of the reordered nmsa that are 0 in 
% the binary msa. 
z_ordered = ordered;
z_ordered(~bin_ordered) = 0;

z_ordered = z_ordered+1; % fastMI can't take 0's.

    % MI matrix
    MI = NMSA_to_fastMI(z_ordered);
    
% Correct for gaps --------------------------------------------------------
if gapW_exp == 0
    fprintf('Gap weights section skipped. \n');
    else
    fprintf('Gap weights section started. \n');
    tic   
    gapW = correct_coevmat_forgaps(nmsa);
    gW = gapW.^gapW_exp;

    MI = MI.*gW;
    gapW_time = toc;
    fprintf('Gap weights section completed in %8.4f seconds. \n', gapW_time);
end        
% -------------------------------------------------------------------------
        
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
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
% format (.nmsa). The function uses  the "fastMI" function. It does not
% apply a correction for gaps.

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
% Mutual information between two column vectors, X and Y, having
% integer values. 

N = size(X,1);
const = log2(N);
h = accumarray([X Y], 1); 
xy_prod = sum(h,2)*sum(h,1);
xy_ind = h~=0;
MI = const + (sum(h(xy_ind) .* log2(h(xy_ind)./xy_prod(xy_ind))))/N;

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



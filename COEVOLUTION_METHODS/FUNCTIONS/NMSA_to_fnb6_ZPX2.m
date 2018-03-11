function [ZPX2,MIP,MI] = NMSA_to_fnb6_ZPX2(nmsa,MI_method)
% This function combines the nb and md algorithms, implementing an improved
% sequence resorting method and a new fast MI calculation. The MI method
% can be:
% odMI: 'one dimensional MI'. This is the standard fastMI calculation.
% mdMI: 'multi-dimensional MI.
% TI: 'total interaction'

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
    switch MI_method
        case 'odMI'
            MI = NMSA_to_fastMI(z_ordered);
        case 'mdMI'
            [MI] = NMSA_to_mdMI(z_ordered);
        case 'TI'
            [~,MI] = NMSA_to_mdMI(z_ordered);            
    end
    
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


function [vMI,tMI] = NMSA_to_mdMI(nmsa)

% This function calculates the 'standard' MI matrix, the 'direct
% coupling' vMI, and the 'total interaction' tMI. The direct coupling
% vMI_3d is the negative of the 'interaction information' minus the mutual
% information. tMI_3d is the 'total interaction' between the three
% variables.

    [nrows,ncols] = size(nmsa);
    const = log2(nrows);
    
    [~,vMI_3d,tMI_3d] = NMSA_to_vMI_3d(nmsa,nrows,ncols,const);
    
    % Correct for gaps along the 3rd dimension-----------------------------
    % vMI_3d_orig = vMI_3d;
    % tMI_3d_orig = tMI_3d;
    
    % gapW = correct_coevmat_forgaps(nmsa);
    % gW = gapW.^gapW_exp;
    % lgW = mean(gW);
    
    % for k = 1:ncols
    % vMI_3d(:,:,k) = vMI_3d(:,:,k)*lgW(k);
    % tMI_3d(:,:,k) = tMI_3d(:,:,k)*lgW(k);
    % end
    % ---------------------------------------------------------------------
   
    vMI = squeeze(nanmean(vMI_3d,3));
    tMI = squeeze(nanmean(tMI_3d,3));
                       
end


function [MI,vMI_3d,tMI_3d] = NMSA_to_vMI_3d(nmsa,nrows,ncols,const)
%--------------------------------------------------------------------------
% MI, vMI_3d, tMI_3d calculation.

MI = NaN(ncols,ncols);
vMI_3d = NaN(ncols,ncols,ncols);
tMI_3d = NaN(ncols,ncols,ncols);

Hi = NaN(ncols,1);
for i = 1:ncols
        h = accumarray(nmsa(:,i), 1);        
        i_marg = h(h~=0);
        
        % Entropy of i
        Hi(i) = const - sum(i_marg.*log2(i_marg))/nrows;
end

Hij = NaN(ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols

        h = accumarray(nmsa(:,[i j]), 1); 
        ij_marg = h(h~=0);

        % Joint Entropy of i and j
        Hij(i,j) = const - sum(ij_marg.*log2(ij_marg))/nrows;         
        Hij(j,i) = Hij(i,j);
        
        % Mutual Information of i and j
        MI(i,j) = Hi(i) + Hi(j) - Hij(i,j);
        MI(j,i) = MI(i,j);
        
    end
end

Hijk = NaN(ncols,ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols
        for k = j+1:ncols
        
            h = accumarray(nmsa(:,[i j k]), 1); 
            ijk_marg = h(h~=0);

            % Joint Entropy of i,j and k
            Hijk(i,j,k) = const - sum(ijk_marg.*log2(ijk_marg))/nrows;         

            % We store also all the permutations. 
            Hijk(i,k,j) = Hijk(i,j,k);
            Hijk(j,i,k) = Hijk(i,j,k);
            Hijk(j,k,i) = Hijk(i,j,k);
            Hijk(k,i,j) = Hijk(i,j,k);
            Hijk(k,j,i) = Hijk(i,j,k);
            
        end        
    end
end

all_ind = (1:ncols);
for i = 1:ncols
    for j = i+1:ncols
        k_ind = setdiff(all_ind,[i j]);
        for k = k_ind
        
            % The following is the interaction information minus the MI.
            % vMI_3d(i,j,k) = Hi(k) - Hij(j,k) - Hij(i,k) + Hijk(i,j,k);
            % We take directly the negative of it, which is equal to 
            % Sunil Srinivasa expression.
            % vMI_3d(i,j,k) = Hij(i,k) - Hijk(i,k,j) - Hi(k) + Hij(k,j);
            % Since we are going to average over the 3rd dimension, here
            % only i and j can be permuted.
            
            vMI_3d(i,j,k) = Hij(j,k) + Hij(i,k) - Hi(k) - Hijk(i,j,k);
            vMI_3d(j,i,k) = vMI_3d(i,j,k);
            tMI_3d(i,j,k) = Hi(i) + Hi(j) + Hi(k) - Hijk(i,j,k);
            tMI_3d(j,i,k) = tMI_3d(i,j,k);
                    
        end        
    end
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



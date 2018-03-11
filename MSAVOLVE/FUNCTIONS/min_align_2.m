[nrows,ncols] = size(nmsa);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
%    dist(i,i) = 0;
end

dist2 = dist==min(dist(:));
[row,col] = find(triu(dist2,1));

ordered = zeros(nrows,ncols);
ordered(1,:) = nmsa(row,:);
ind = 1;
for i = 2:nrows
    % zero the entire test column
    dist(:,ind) = 0;
    % find the index of the closest sequence
    [~,newind] = max(dist(ind,:));
    ordered(i,:) = nmsa(newind,:);
    % the next index will be the index of the sequence just selected
    ind = newind;
end

ordered_dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        ordered_dist(i,j) = sum(ordered(i,:)==ordered(j,:));
        ordered_dist(j,i) = ordered_dist(i,j);
    end
    ordered_dist(i,i) = ncols;
end

figure;imagesc(ordered_dist);

% Here we convert the msa from matlab numeric format to binary
% differential.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end

% Here we select the best level in the stack based on its entropy and we 
% calculate the best ZPX2 matrix from that level.

    % MI matrix and mean column entropy
    [MI,pH] = BNMSA_to_MI(bin_ordered,nrows,ncols);
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);

    
    fprintf('a. binary MSA sum = %f \n', sum(sum(bin_ordered)));
    fprintf('b. posterior mean binary MSA entropy = %f \n', pH);
    
    
    

%figure;imagesc(dist);
inv_dist = 1-(dist/236);
phylotree = seqlinkage(inv_dist);
view(phylotree);

figure;imagesc(ordered_dist);
inv_ordered_dist = 1-(ordered_dist/236);
phylotree_ordered = seqlinkage(inv_ordered_dist);
view(phylotree_ordered);

figure;imagesc(bin_ordered_dist);
%inv_bin_ord_dist = 1-(bin_ordered_dist/max(bin_ordered_dist(:)));
%figure;imagesc(inv_bin_ord_dist);
phylotree_bin_ord = seqlinkage(bin_ordered_dist);
view(phylotree_bin_ord);
plot(phylotree_bin_ord,'Type','equaldaylight')


function [pMI,mH] = BNMSA_to_MI(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% MI calculation; the following are the rules used in generating the
% arrays.
% pI(1) == 1 ; pI(2) == 0 
% pJ(1) == 1 ; pJ(2) == 0
% pIJ(1) == 1,1
% pIJ(2) == 1,0
% pIJ(3) == 0,1
% pIJ(4) == 0,0
bnmsa = + bnmsa;
const = log2(nrows);
H = zeros(1,ncols);
pMI = zeros(ncols,ncols);
pIJ = zeros(4,1);
epIJ = zeros(4,1);

pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
end
mH = mean(H);

for I = 1:ncols
    for J = I:ncols
        pIJ(1) = bnmsa(:,I)'*bnmsa(:,J);
        pIJ(2) = pI(1,I) - pIJ(1);        
        pIJ(3) = pI(1,J) - pIJ(1);
        pIJ(4) = nrows - pIJ(1) - pIJ(2) - pIJ(3);
        epIJ(1) = pI(1,I)*pI(1,J);
        epIJ(2) = pI(1,I)*pI(2,J);
        epIJ(3) = pI(2,I)*pI(1,J);
        epIJ(4) = pI(2,I)*pI(2,J);
        ind = find(pIJ);
        pMI(I,J) = (sum(pIJ(ind) .* (log2(pIJ(ind) ./ epIJ(ind)) + const)))/nrows;
    end
end
pMI = pMI+pMI';
for m = 1:ncols
pMI(m,m)=NaN;
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



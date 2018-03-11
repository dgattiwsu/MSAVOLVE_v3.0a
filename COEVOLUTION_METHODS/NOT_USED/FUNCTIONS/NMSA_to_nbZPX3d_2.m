function [ZPX3d,MIP3d,MI3d] = NMSA_to_nbZPX3d_2(nmsa)
% In this function we reorder the msa such that there is the minimum 
% distance between each sequence and the following one.
% First we determine the distance matrix: here more distant sequences are
% represented by smaller numbers.

[nrows,ncols] = size(nmsa);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
end

triu_dist = triu(dist,1);
[~,max_ind] = max(triu_dist(:));
[firstrow,firstcol] = ind2sub([nrows,nrows],max_ind);

% Here we select as the top sequence the sequence with the smallest
% entropy. Not clear whether the following 4 lines make a difference.
% rowsum = Entropy(nmsa(firstrow,:)');
% colsum = Entropy(nmsa(firstcol,:)');
% if rowsum >= colsum
%     firstrow = firstcol;
% end

nmsa1 = nmsa(firstrow:nrows,:);
nmsa2 = nmsa(1:firstrow-1,:);
nmsa = [nmsa1;nmsa2];

%--------------------------------------------------------------------------
% Here we reorder based on the first algorithm.
ordered = reorder_nmsa_1(nmsa,nrows,ncols);
%--------------------------------------------------------------------------

% Here we convert the msa from matlab numeric format to binary
% differential.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end

%Here we set to zero all the indices of the reordered nmsa that are 0 in 
% the binary msa. 
z_ordered = ordered;
z_ordered(~bin_ordered) = 0;

% Here we calculate the ZPX2 matrix.

    % MI matrix and mean column entropy
    [MI3d] = NMSA_to_MI3d(z_ordered,nrows,ncols);
    % MIP3d matrix
    MIP3d = MI_to_MIP3d(MI3d,ncols);
    % and finally ZPX3 matrix
    ZPX3d = MIP_to_ZPX3d(MIP3d,ncols);
% 
%     fprintf('a. binary MSA sum = %f \n', sum(sum(bin_ordered)));
%     fprintf('b. posterior mean binary MSA entropy = %f \n', pH);
            
end


function [ordered] = reorder_nmsa_1(nmsa,nrows,ncols)
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


function [ordered] = reorder_nmsa_3(nmsa,nrows,ncols)
%--------------------------------------------------------------------------
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)~=nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
end

ordered = zeros(nrows,ncols);

% The 1st sequence is set by default
ordered(1,:) = nmsa(1,:);
ind = 1;

for i = 2:nrows
    % zero the entire test column
    dist(:,ind) = NaN;
    % find the index of the closest sequence
    [~,newind] = min(dist(ind,:));
    ordered(i,:) = nmsa(newind,:);
    % the next index will be the index of the sequence just selected
    ind = newind;
end

end

function [MI3d] = NMSA_to_MI3d(nmsa,nrows,ncols)
%--------------------------------------------------------------------------
% MI3d calculation.

MI3d = zeros(ncols,ncols,ncols);
MI3d_unique = zeros(ncols,ncols,ncols);

Hi = NaN(ncols);
for i = 1:ncols
        Hi(i) = Entropy(nmsa(:,i));
end

Hij = NaN(ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols
        
        Hij(i,j) = JointEntropy(nmsa(:,[i j]));
        
    end
end

Hijk = NaN(ncols,ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols        
        for k = j+1:ncols
        
            Hijk(i,j,k) = JointEntropy(nmsa(:,[i j k]));
            
            MI3d(i,j,k) = Hi(i) + Hi(j) + Hi(k) - ...
                (Hij(i,j) + Hij(j,k) + Hij(i,k)) + Hijk(i,j,k);
                
        % This is the unique part of the matrix.
        MI3d_unique(i,j,k) = MI3d(i,j,k);
        
        % We store also all the permutations. This will be useful in
        % calculating the ZPX3 matrix.
        
        MI3d(i,k,j) = MI3d(i,j,k);
        MI3d(j,i,k) = MI3d(i,j,k);
        MI3d(j,k,i) = MI3d(i,j,k);
        MI3d(k,i,j) = MI3d(i,j,k);
        MI3d(k,j,i) = MI3d(i,j,k);
        
        end        
    end
end

% zero_ind = MI3d == 0;
% MI3d(zero_ind) = NaN;

end


function [pMI3d] = BNMSA_to_MI3d(bnmsa,nrows,ncols)
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
% const = log2(nrows);
% H = zeros(1,ncols);
% pMI = zeros(ncols,ncols);
const2 = log2(nrows^2);
pMI3d = zeros(ncols,ncols,ncols);
pMI3d_unique = zeros(ncols,ncols,ncols);

% pIJ = zeros(4,1);
% epIJ = zeros(4,1);
pIJ = zeros(ncols,ncols,4);
epIJ = zeros(ncols,ncols,4);
pIJK = zeros(ncols,ncols,ncols,8);
epIJK = zeros(ncols,ncols,ncols,8);


pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
%         ind = find(pI(:,I));
%         H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
end
% mH = mean(H);

for I = 1:ncols
    for J = I+1:ncols
        
        pIJ(I,J,1) = bnmsa(:,I)'*bnmsa(:,J);
        pIJ(I,J,2) = pI(1,I) - pIJ(I,J,1);        
        pIJ(I,J,3) = pI(1,J) - pIJ(I,J,1);
        pIJ(I,J,4) = nrows - pIJ(I,J,1) - pIJ(I,J,2) - pIJ(I,J,3);
        
        epIJ(I,J,1) = pI(1,I)*pI(1,J);
        epIJ(I,J,2) = pI(1,I)*pI(2,J);
        epIJ(I,J,3) = pI(2,I)*pI(1,J);
        epIJ(I,J,4) = pI(2,I)*pI(2,J);
        
%         ind_orig = pIJ(I,J,:);
%         ind = find(ind_orig);
% %        ind = find(pIJ);
%         pMI(I,J) = (sum(pIJ(I,J,ind) .* (log2(pIJ(I,J,ind) ...
%             ./ epIJ(I,J,ind)) + const)))/nrows;
        
% Next we calculate the joint probability of each triplet
% % pIJK(1) == 0,0,0
% % pIJK(2) == 0,0,1
% % pIJK(3) == 0,1,0
% % pIJK(4) == 0,1,1
% % pIJK(5) == 1,0,0
% % pIJK(6) == 1,0,1
% % pIJK(7) == 1,1,0
% % pIJK(8) == 1,1,1

        for K = J+1:ncols
        pIJK(I,J,K,:) = JointProb3(bnmsa(:,[I J K]));
        epIJK(I,J,K,1) = epIJ(I,J,4)*pI(2,K); % 0 0 0
        epIJK(I,J,K,2) = epIJ(I,J,4)*pI(1,K); % 0 0 1
        epIJK(I,J,K,3) = epIJ(I,J,3)*pI(2,K); % 0 1 0
        epIJK(I,J,K,4) = epIJ(I,J,3)*pI(1,K); % 0 1 1
        epIJK(I,J,K,5) = epIJ(I,J,2)*pI(2,K); % 1 0 0
        epIJK(I,J,K,6) = epIJ(I,J,2)*pI(1,K); % 1 0 1
        epIJK(I,J,K,7) = epIJ(I,J,1)*pI(2,K); % 1 1 0
        epIJK(I,J,K,8) = epIJ(I,J,1)*pI(1,K); % 1 1 1
        
        ind_orig = pIJK(I,J,K,:);
        ind = find(ind_orig);
              
        pMI3d(I,J,K) = ...
            (sum(pIJK(I,J,K,ind) .* (log2(pIJK(I,J,K,ind) ...
            ./ epIJK(I,J,K,ind)) + const2)))/nrows;
        
        % This is the unique part of the matrix.
        pMI3d_unique(I,J,K) = pMI3d(I,J,K);
        
        % We store also all the permutations. This will be useful in
        % calculating the ZPX3 matrix.
        pMI3d(I,K,J) = pMI3d(I,J,K);
        pMI3d(J,I,K) = pMI3d(I,J,K);
        pMI3d(J,K,I) = pMI3d(I,J,K);
        pMI3d(K,I,J) = pMI3d(I,J,K);
        pMI3d(K,J,I) = pMI3d(I,J,K);
        
        end        
    end
end
zero_ind = pMI3d == 0;
pMI3d(zero_ind) = NaN;
% pMI = pMI+pMI';
% for m = 1:ncols
% pMI(m,m)=NaN;
% end

end


function [JP] = JointProb2(X)
% % pIJ(1) == 1,1
% % pIJ(2) == 1,0
% % pIJ(3) == 0,1
% % pIJ(4) == 0,0
JP = zeros(4,1);
A1 = [1 1];
A2 = [1 0];
A3 = [0 1];
A4 = [0 0];
JP(1) = sum(X(:,1)==A1(1) & X(:,2)==A1(2));
JP(2) = sum(X(:,1)==A2(1) & X(:,2)==A2(2));
JP(3) = sum(X(:,1)==A3(1) & X(:,2)==A3(2));
JP(4) = sum(X(:,1)==A4(1) & X(:,2)==A4(2));
end


function [JP] = JointProb3(X)
% % pIJK(1) == 0,0,0
% % pIJK(2) == 0,0,1
% % pIJK(3) == 0,1,0
% % pIJK(4) == 0,1,1
% % pIJK(5) == 1,0,0
% % pIJK(6) == 1,0,1
% % pIJK(7) == 1,1,0
% % pIJK(8) == 1,1,1
JP = zeros(8,1);
JP(1) = sum(X(:,1)==0 & X(:,2)==0 & X(:,3)==0);
JP(2) = sum(X(:,1)==0 & X(:,2)==0 & X(:,3)==1);
JP(3) = sum(X(:,1)==0 & X(:,2)==1 & X(:,3)==0);
JP(4) = sum(X(:,1)==0 & X(:,2)==1 & X(:,3)==1);
JP(5) = sum(X(:,1)==1 & X(:,2)==0 & X(:,3)==0);
JP(6) = sum(X(:,1)==1 & X(:,2)==0 & X(:,3)==1);
JP(7) = sum(X(:,1)==1 & X(:,2)==1 & X(:,3)==0);
JP(8) = sum(X(:,1)==1 & X(:,2)==1 & X(:,3)==1);
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


function [pMIP] = MI_to_MIP3d(pMI,ncols)
%--------------------------------------------------------------------------
% MIP calculation
mean_mat = nanmean(pMI(:));
MCA_mat = NaN(ncols,ncols,ncols);

% Here  we calculate the MCA matrix. We recall that the MI3d matrix is
% symmetric in all 3 directions, so only one mean plane is necessary.

mean_plane = squeeze(nanmean(pMI,3));

for m = 1:ncols
    for n = m:ncols
        for o = n:ncols
        MCA_mat(m,n,o)=(mean_plane(m,n)*mean_plane(m,o)*mean_plane(n,o))/mean_mat;
        MCA_mat(m,o,n) = MCA_mat(m,n,o);    
        MCA_mat(n,m,o) = MCA_mat(m,n,o);    
        MCA_mat(n,o,m) = MCA_mat(m,n,o);    
        MCA_mat(o,m,n) = MCA_mat(m,n,o);    
        MCA_mat(o,n,m) = MCA_mat(m,n,o);    
        end
    end
end

% Finally we subtract the MCA matrix from the MI matrix.
pMIP = pMI-MCA_mat;

end


function [pZPX3] = MIP_to_ZPX3d(pMIP,ncols)
%--------------------------------------------------------------------------
% ZPX3 calculation

pZPX3 = NaN(ncols,ncols,ncols);

% Here  we calculate the ZPX3d matrix. We recall that the MIP3d matrix is
% symmetric in all 3 directions, so only one mean plane is necessary.

mean_plane = squeeze(nanmean(pMIP,3));
std_plane = squeeze(nanstd(pMIP,0,3));

for m = 1:ncols
    for n = m:ncols
        for o = n:ncols

        ZPX3_i = (pMIP(m,n,o)-mean_plane(m,n))/std_plane(m,n);
        ZPX3_j = (pMIP(m,n,o)-mean_plane(m,o))/std_plane(m,o);
        ZPX3_k = (pMIP(m,n,o)-mean_plane(n,o))/std_plane(n,o);
        
        pZPX3(m,n,o) = ZPX3_i*ZPX3_j*ZPX3_k;

        pZPX3(m,o,n) = pZPX3(m,n,o);
        pZPX3(n,m,o) = pZPX3(m,n,o);
        pZPX3(n,o,m) = pZPX3(m,n,o);
        pZPX3(o,m,n) = pZPX3(m,n,o);
        pZPX3(o,n,m) = pZPX3(m,n,o);
        
        end    
    end   
end

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




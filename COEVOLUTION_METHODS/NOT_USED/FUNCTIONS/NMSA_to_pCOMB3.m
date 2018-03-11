function [ZPX2_stack,CORR_stack,COMB_stack,MSA_test_stack,...
    MSA_test_pos_stack,min_sum,hits] = ...
    NMSA_to_pCOMB3(nmsa,threshold,covar_vec)
% This function calculates a CORRELATION matrix (bCORR), a ZPX2 matrix 
% (bZPX2), and a matrix that combines the previous two (bCOMB). 
% Usage:
% [bZPX2,bCORR,bCOMB,MSA_test_stack,MSA_test_pos_stack,min_sum] = ...
%       NMSA_to_pCOMB2(nmsa,0.7,1,2);
%
% nmsa is the msa in matalb numeric format.
%
% The recommended value for threshold is 0.7; the calculated weight for  
% each MI map (currently not being used) is 1/k where k is the number of 
% sequences that are more than 70% (for threshold = 0.7) identical to the 
% reference sequence.
%
% 'entropy method' can have values 1 (mean entropy) or 2 (global entropy).
% 'entropy level' can have values 1 (lowest entropy) or 2 (highest
% entropy.
% 
% Place a ~ in the usage syntax for each item you don't want: for example
% if you only want the bCOMB matrix based on the reordered msa with the 
% highest value of global entropy use the following:
% [~,~,bCOMB,~,~,~] = NMSA_to_pCOMB2(nmsa,0.7,2,2);
%
% Don't forget the ";" at the end or the output will go on forever on the
% screen!!
%
% The outputs are:
% bZPX2: best pZPX2 matrix.
% bCORR: best correlation matrix.
% bCOMB: best combination matrix.
% MSA_test_stack: a stack of all the reordered msa in standard format.
% MSA_test_pos_stack: a stack of all the reordered msa in binary format.
% min_sum: 
% 1st column = values of the sum of 1's in the binary matrices. 
% 2nd column = mean entropy values of the same binary matrices. 
% 3rd column = global entropy values of the same binary matrices.
%
%--------------------------------------------------------------------------
% Here we calculate weights and Meff like in DCA; this part is not
% currently used for anything, but helps to have an internal comparison
% with DCA.

[nrows,ncols] = size(nmsa);
W = ones(nrows,1);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:))/ncols;
        dist(j,i) = dist(i,j);
    end
end
dist_threshold = dist > threshold; 
    W = 1./sum(dist_threshold)';
Meff=sum(W);
fprintf('Meff = %f \n', Meff);

MSA_test = nmsa;
MSA_test_bk = MSA_test;

% Here we calculate the initial value of the pMSA_sum

pMSA = false(nrows,ncols);

for i = 2:nrows
    pMSA(i,:) = MSA_test(i,:)~=MSA_test(i-1,:);
end

mH = BNMSA_to_mH(pMSA,nrows,ncols);
% mrH = BNMSA_to_mrH(pMSA,nrows,ncols);
mrH = BNMSA_to_mrH_2(pMSA,nrows,ncols);
% gH = BNMSA_to_gH(pMSA,nrows,ncols);
fprintf('prior binary MSA sum = %f \n', sum(pMSA(:)));
fprintf('prior mean binary MSA entropy = %f \n', mH);
% fprintf('prior global binary MSA entropy = %f \n', gH);
fprintf('prior mean binary MSA rel_entropy = %f \n', mrH);

clear pMSA

%-----------Here we run a full optimization of the binary matrix-----------
%
fprintf('Running full optimization \n');

% Important variables.
% MSA_test_stack: stack of reordered MSAs in standard numeric format.
% MSA_test_pos_stack: stack of reordered binary MSAs.

MSA_test = MSA_test_bk;
MSA_test_stack = zeros(nrows,ncols,nrows);
MSA_test_pos_stack = zeros(nrows,ncols,nrows);
% MSA_test_Q_stack = zeros(nrows,1);
% MSA_test_W_stack = zeros(nrows,nrows);
H_stack = zeros(nrows,2);
ZPX2_stack = zeros(ncols,ncols,nrows);
CORR_stack = zeros(ncols,ncols,nrows);
COMB_stack = zeros(ncols,ncols,nrows);
% DIST_stack = zeros(nrows,nrows,nrows);
NORM_stack = zeros(nrows,2);
hits = zeros(nrows,3);


min_sum = zeros(nrows,1);
order = 1:nrows;
reorder = zeros(1,nrows);

% Here we loop over all possible sequences in the msa as the reference
% sequence in the minimization process.

for j = 1:nrows
 tStart = tic;   
 fprintf('looping over sequence #  %d \n', j);
    
    MSA_test = MSA_test_bk(order,:);

%--------------------------------------------------------------------------
% Here we find the difference between each sequence and the current test
% sequence and we reorder the msa such that there is going to be minimal
% difference between each sequence and the following one.
% In the following loop "MSA_test_stack" stores all the reordered nmsa's in
% the standard 25 numeric symbols
MSA_test = reorder_nmsa_1(MSA_test,nrows,ncols);

%     for k = 1:nrows-1
%     MSA_test_dif = zeros(nrows,1);
%         for i = k:nrows
%         MSA_test_dif(i) = sum(MSA_test(k,:)~=MSA_test(i,:));
%         end
%     [~, sort_MSA_ind] = sort(MSA_test_dif,'ascend');
%     MSA_test_2 = MSA_test(sort_MSA_ind,:);
%     MSA_test = MSA_test_2;
%     end
% 
%--------------------------------------------------------------------------
% The following function is slower,but the code is much more transparent. 
% MSA_test = reorder_nmsa(MSA_test,nrows,ncols);
%--------------------------------------------------------------------------

% Store the reordered msa in the stack
MSA_test_stack(:,:,j) = MSA_test;

% We calculate the distance matrix for each of the reordered msa's

% DIST_stack(:,:,j) = distance_matrix(MSA_test);

% Here we calculate the binary msa.
% In the following loop "MSA_test_pos_stack" stores all the reordered
% binary nmsa's
    
    MSA_test_pos = false(nrows,ncols);
        for i = 2:nrows
            MSA_test_pos(i,:) = MSA_test(i,:)~=MSA_test(i-1,:);
        end
    min_sum(j,1) = sum(MSA_test_pos(:));
    MSA_test_pos_stack(:,:,j) = MSA_test_pos;
    
    NORM_stack(j,1) = norm(+MSA_test_pos);
    NORM_stack(j,2) = norm(+MSA_test_pos,'fro');


% Here we calculate the H_stack and the ZPX2 stack
    % Global entropy
    % gH = BNMSA_to_gH(MSA_test_pos,nrows,ncols);
    mrH = BNMSA_to_mrH_2(MSA_test_pos,nrows,ncols);
    % MI matrix and mean column entropy
    [MI,mH] = BNMSA_to_MI(MSA_test_pos,nrows,ncols);
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
    H_stack(j,1) = mH;
    % H_stack(j,2) = gH;
    H_stack(j,2) = mrH;
    ZPX2_stack(:,:,j) = ZPX2;
    
% Here we calculate the CORR stack    
    CORR = BNMSA_to_CORR(MSA_test_pos,nrows,ncols);
    % MIP matrix
    CORR = MI_to_MIP(CORR,ncols);
    % and finally ZPX2 matrix
    CORR = MIP_to_ZPX2(CORR,ncols);

    CORR_stack(:,:,j) = CORR;
    
% Here we calculate the COMB stack    

    count = (ncols*ncols-ncols)/2;
    lin_ZPX2 = zeros(count,1);
    lin_CORR = zeros(count,1);
    count = 0;
    for ii = 2:ncols
        for jj = ii:ncols
            count = count+1;
        	lin_ZPX2(count) = ZPX2(ii-1,jj);
            lin_CORR(count) = CORR(ii-1,jj);
        end
    end
    lin_ZPX2 = nantozero(lin_ZPX2);
    lin_CORR = nantozero(lin_CORR);

% scale = (l_bZPX2'*l_bCORR)/(l_bCORR'*l_bCORR)
    scale = (lin_ZPX2'*lin_CORR)/(lin_ZPX2'*lin_ZPX2);

% fprintf('best scale CORR/bZPX2 = %f \n', scale);

% Finally we combine the two matrices.

    COMB = ZPX2 + CORR/scale;
    COMB_stack(:,:,j) = COMB;    
    
% Here we calculate the success rate
sorted_CORR = sort_matrix_descend(CORR);
CORR_covar_vec = sorted_CORR(1:50,2:3); 
hits(j,1) = numel(intersect(covar_vec(:,1:2),CORR_covar_vec,'rows'))/2;

sorted_ZPX2 = sort_matrix_descend(ZPX2);
ZPX2_covar_vec = sorted_ZPX2(1:50,2:3); 
hits(j,2) = numel(intersect(covar_vec(:,1:2),ZPX2_covar_vec,'rows'))/2;

sorted_COMB = sort_matrix_descend(COMB);
COMB_covar_vec = sorted_COMB(1:50,2:3); 
hits(j,3) = numel(intersect(covar_vec(:,1:2),COMB_covar_vec,'rows'))/2;


    reorder(1:nrows-1) = order(2:nrows);
    reorder(nrows) = order(1);
    order = reorder;

% Keep track of the time.    
    tElapsed = toc(tStart);
    fprintf('execution time for sequence %d : %6.2f s \n', j,tElapsed);
    
end

% Here we generate the entire stack.
% bZPX2_stack: stack of reordered coevolution matrices.
% H_stack: vector of mean (or global) entropy values for all the binary msa's.
% W_stack: stack of weighted reordered coevolution matrices.

% CORR_stack = zeros(ncols,ncols,nrows);
% H_stack = zeros(nrows,2);
% w_stack = zeros(ncols,ncols,nrows);
    
%    fprintf('Running full optimization at level 2 \n');        

%    for i = 1:nrows
%        fprintf('Matrix # %d for pMSAsum %f \n',i,min_sum(i));
%        pMSA = MSA_test_pos_stack(:,:,i);

% Here we calculate the mean entropy of all columns.
%        H_stack(i,1) = BNMSA_to_mH(pMSA,nrows,ncols);
        
% Here we calculate the global entropy for the entire binary msa. The result 
% of the following function is the same as matlab "entropy" or by 
% calculating the entropy from the partition function.
%        H_stack(i,2) = BNMSA_to_gH(pMSA,nrows,ncols);

% Here we calculate the correlation or the covariance map.
%        CORR_stack(:,:,i) = BNMSA_to_CORR(pMSA,nrows,ncols);
%        CORR_stack(:,:,i) = BNMSA_to_COV(pMSA,nrows,ncols);

%        clear pMSA 
%    end
    
% We first apply weights to each level of the stack.
%     for i =1:nrows        
%         w_stack(:,:,i) = CORR_stack(:,:,i)*W(i);
%     end

%     w_stack = CORR_stack;
    
% Then we calculate mean and std.
%   CORR_std = nanstd(w_stack,1,3);

%     CORR_mean = nanmean(w_stack,3);
    
% Here we calculate the MIP map from the mean CORR map.

%   CORR_mip = MI_to_MIP(CORR_mean,ncols);
    
% We merge minimization and entropy data for each run

    min_sum = [min_sum H_stack NORM_stack];
%    min_sum = [min_sum H_stack MSA_test_Q_stack];
%    min_sum_product = min_sum(:,1).*min_sum(:,2);
%    min_sum = [min_sum min_sum_product];
    z_min_sum = zscore(min_sum);
    min_sum = [min_sum z_min_sum];
    
% Here we select the best level in the stack based on its entropy and we
% calculate the ZPX2 matrix for that level

%    best_ind = min_sum(:,3) == min(min_sum(:,3));
% if entropy_method == 1 && entropy_level == 1
%         best_ind = min_sum(:,2) == min(min_sum(:,2));
%    fprintf('best map is from msa with the lowest mean entropy \n');
% 
% elseif entropy_method == 1 && entropy_level == 2
%         best_ind = min_sum(:,2) == max(min_sum(:,2));
%    fprintf('best map is from msa with the highest mean entropy \n');
% 
% elseif entropy_method == 2 && entropy_level == 1
%         best_ind = min_sum(:,3) == min(min_sum(:,3));
%    fprintf('best map is from msa with the lowest global entropy \n');
% else
%         best_ind = min_sum(:,3) == max(min_sum(:,3));
%    fprintf('best map is from msa with the highest global entropy \n');
% end

%    bpNMSA = MSA_test_stack(:,:,best_ind);

%     % MI matrix
%     bMI = BNMSA_to_MI(MSA_test_pos_stack(:,:,best_ind),nrows,ncols);
%     % MIP matrix
%     bMIP = MI_to_MIP(bMI,ncols);
%     % and finally ZPX2 matrix
%     bZPX2 = MIP_to_ZPX2(bMIP,ncols);

% Here we use the best level also for the correlation. This line overrides
% previous assignment of CORR if there was any.

%      bCORR = BNMSA_to_CORR(MSA_test_pos_stack(:,:,best_ind),nrows,ncols);
     
%    CORR_mean = w_stack(:,:,best_ind);
%    CORR_mean_mip = MI_to_MIP(CORR_mean,ncols);
%    CORR_mean_zpx2 = MIP_to_ZPX2(CORR_mean_mip,ncols);
    
%    CORR_mean = CORR_mean_mip;

% Here we merge MI and correlation matrix by linear regression. This step
% is questionable: we would probably need a "local" scaling for different
% sections of the matrices.

% count = (ncols*ncols-ncols)/2;
% l_bZPX2 = zeros(count,1);
% l_bCORR = zeros(count,1);
% count = 0;
% for i = 2:ncols
%     for j = i:ncols
%         count = count+1;
%         l_bZPX2(count) = bZPX2(i-1,j);
%         l_bCORR(count) = bCORR(i-1,j);
%     end
% end
% l_bZPX2 = nantozero(l_bZPX2);
% l_bCORR = nantozero(l_bCORR);
% 
% % scale1 = (l_bZPX2'*l_bCORR)/(l_bCORR'*l_bCORR)
% scale2 = (l_bZPX2'*l_bCORR)/(l_bZPX2'*l_bZPX2);
% 
%     fprintf('best scale CORR/bZPX2 = %f \n', scale2);
% 
% % Finally we combine the two matrices.
% 
% % COMB = bZPX2 + CORR_mean;
% % COMB = bZPX2 + scale1*CORR_mean;
% 
%  bCOMB = bZPX2 + bCORR/scale2;
%  
% %     COMB = MI_to_MIP(COMB,ncols);
% %     COMB = MIP_to_ZPX2(COMB,ncols);


%     fprintf('best COMB map is # %d corresponding to:\n', find(best_ind));
%     fprintf('a. posterior binary MSA sum = %f \n', min_sum(best_ind,1));
%     fprintf('b. posterior mean binary MSA entropy = %f \n', min_sum(best_ind,2));
%     fprintf('c. posterior global binary MSA entropy = %f \n', min_sum(best_ind,3));

% Here we filter the COMB matrix taking only the top # eigenvectors.
% mCOMB = COMB;
% nan_COMB = isnan(COMB); 
% mCOMB(nan_COMB) = 0;
% [V_COMB, D_COMB] = spectral_2(mCOMB,100);
% eCOMB = V_COMB*D_COMB*V_COMB';
% eCOMB(nan_COMB) = NaN;
%---End of full optimization-----------------------------------------------

clear pMSA MSA_test MSA_test_2 sort_MSA_ind

end


function [W] = nmsa_to_weights(nmsa,threshold,nrows,ncols)
%--------------------------------------------------------------------------
% Weights calculation for the partition function
W = ones(nrows,1);
dist = zeros(nrows,nrows);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:))/ncols;
        dist(j,i) = dist(i,j);
    end
end
dist_threshold = dist > threshold; 
    W = 1./sum(dist_threshold)';
end


function [mH] = BNMSA_to_mH(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% Binary Entropy calculation.
bnmsa = + bnmsa;
const = log2(nrows);
H = zeros(1,ncols);

pI = zeros(2,ncols);
for I = 1:ncols
        pI(1,I) = sum(bnmsa(:,I));
        pI(2,I) = nrows - pI(1,I);
        ind = find(pI(:,I));
        H(I) = -(sum(pI(ind,I).*(log2(pI(ind,I)) - const)))/nrows;
end
mH = mean(H);

end


function [H] = BNMSA_to_gH(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% Binary Entropy calculation.
bnmsa = + bnmsa(:);
nrows = nrows*ncols;
const = log2(nrows);

% Even faster
pI = zeros(2,1);
        pI(1) = sum(bnmsa);
        pI(2) = nrows - pI(1);
        ind = find(pI);
        H = -(sum(pI(ind).*(log2(pI(ind)) - const)))/nrows;
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


function [ordered] = reorder_nmsa_2(nmsa,nrows,ncols)
%--------------------------------------------------------------------------
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)==nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
    dist(i,i) = 0;
end

ordered = zeros(nrows,ncols);

% The 1st sequence is set by default
ordered(1,:) = nmsa(1,:);
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


function [dist] = distance_matrix(nmsa)
[nrows,~] = size(nmsa);
for i = 1:nrows
    for j = i:nrows
        dist(i,j) = sum(nmsa(i,:)~=nmsa(j,:));
        dist(j,i) = dist(i,j);
    end
end
end


function [pCORR] = BNMSA_to_CORR(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% CORR calculation; the following are the rules used in generating the
bnmsa = + bnmsa;
    pCORR = corr(bnmsa);
for m=1:ncols
    pCORR(m,m) = NaN;
end

end


function [pCOV] = BNMSA_to_COV(bnmsa,nrows,ncols)
%--------------------------------------------------------------------------
% CORR calculation; the following are the rules used in generating the
bnmsa = + bnmsa;
    pCOV = cov(bnmsa);
for m=1:ncols
    pCOV(m,m) = NaN;
end

end


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




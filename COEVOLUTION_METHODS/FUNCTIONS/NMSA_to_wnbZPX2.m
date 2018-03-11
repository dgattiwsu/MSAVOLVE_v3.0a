function [ZPX2,W,MIP,MI,ordered,bin_ordered,z_ordered] = ...
    NMSA_to_wnbZPX2(nmsa,sort_method,threshold,lambda,gaps_method,nsymbols)
% This function differs from 'nb2_ZPX2' as it implements a weighing scheme
% to account for the similarity between sequences and also uses a
% psudocount ('lambda'). If threshold = 1, all sequences receive a weight
% of 1, and if lambda is set to zero no pseudocount is applied. The flag
% 'gaps_method' can be set to 'INCLUDE' or 'EXCLUDE' to include or exclude
% gaps in the calculation. 'nsymbols' is the number of symbols used in
% scaling the lambda parameter for the pseudocount. Recommended use: [ZPX2]
% = NMSA_to_wnbZPX2(nmsa,'DISTANCE',0.8,0.5,'INCLUDE',21); or [ZPX2] =
% NMSA_to_wnbZPX2(nmsa,'DISTANCE',0.8,0.5,'EXCLUDE',20); It becomes
% identical to NMSA_to_nbZPX2 if run as: [ZPX2] =
% NMSA_to_wnbZPX2(nmsa,'DISTANCE',1.0,0.0,'INCLUDE',21);

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
bin_ordered = nmsa_to_binmsa_21q(nmsa);
dist = bin_ordered*bin_ordered';
[W,Meff] = nmsa_to_W(nmsa,dist,threshold);

switch sort_method
%--------------------------------------------------------------------------
    case 'DISTANCE'
    % Here we reorder based on the distance matrix.
    triu_dist = triu(dist,1);
    [~,max_ind] = max(triu_dist(:));
    [firstrow,firstcol] = ind2sub([nrows,nrows],max_ind);

    nmsa1 = nmsa(firstrow:nrows,:);
    nmsa2 = nmsa(1:firstrow-1,:);
    nmsa = [nmsa1;nmsa2];

    ordered = reorder_nmsa_1(nmsa,nrows,ncols);
%--------------------------------------------------------------------------
    case 'EIGEN'
    % Here we reorder based on the eigenvector corresponding to the  
    % largest eigenvalue of the distance matrix.
    [V,~] = eigs(dist,1,'lm');
    % [V,D] = eig(dist);
    % for i=1:nrows
    %     switch_sign = sign(mean(V(:,i)));
    %     V(:,i) = switch_sign*V(:,i);
    % end
    % V = sum(V,2);
    % switch_sign = sign(mean(V));
    % V = switch_sign*V;
    [~,ind] = sort(abs(V),'descend');
    ordered = nmsa(ind,:);
%--------------------------------------------------------------------------
end

% Here we convert the msa from matlab numeric format to binary
% differential.
bin_ordered = false(nrows,ncols);
    for i = 2:nrows
        bin_ordered(i,:) = ordered(i,:)~=ordered(i-1,:);
    end
    
% Here we set to zero all the indices of the reordered nmsa that are 0 in 
% the binary msa. 
z_ordered = ordered;
z_ordered(~bin_ordered) = 0;

% Here we convert to 0 all the positions that have a value of 21.
% This is equivalent to assuming the each time there is a gap in the real
% msa, the position does not change in the reordered msa.

switch gaps_method
    case 'EXCLUDE'
        ind_21 = z_ordered == 21;
        z_ordered(ind_21) = 0;
    case 'INCLUDE'
end

    % MI matrix
    [MI] = NMSA_to_MI(z_ordered,W,lambda,nsymbols);
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % and finally ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
%    fprintf('a. binary MSA sum = %f \n', sum(sum(bin_ordered)));
%    fprintf('b. posterior mean binary MSA entropy = %f \n', pH);
            
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


function [scale] = scale_matrices(mat1,mat2)
% This function returns the best scale for mat1/mat2.
[rows,cols] = size(mat1);
template = ones(rows,cols);
upper = triu(template,1);
triu_ind = template == upper;
mat1 = mat1(triu_ind);
mat2 = mat2(triu_ind);
scale = mat1\mat2;
% scale1 = mat1\mat2;
% scale1 = (mat1'*mat2)/(mat1'*mat1);
% scale2 = mat2\mat1;
% scale2 = (mat1'*mat2)/(mat2'*mat2);
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


function [ MI_mat ] = NMSA_to_MI(nmsa,W,lambda,nsymbols)
% This function produce a MI matrix starting from a MSA in matlab numeric
% format (.nmsa). The function uses  the "MutualInformation" function from 
% Dwinnel package. It differs from 'NMSA_to_gcMI' because it does not apply 
% a correction for gaps. It is the same function as NMSA_to_ngcMI, where
% ngc = no gap correction.

[~,ncols] = size(nmsa);
MI_mat = zeros(ncols,ncols);

for i = 1:ncols
    for j = i:ncols
    MI_mat(i,j) = MutualInformation(nmsa(:,i),nmsa(:,j),W,lambda,nsymbols);
    MI_mat(j,i) = MI_mat(i,j);    
    end
    MI_mat(i,i)=NaN;
end

end


function I = MutualInformation(X,Y,W,lambda,nsymbols)

% MutualInformation: returns mutual information (in bits) of the 'X' and
% 'Y'. Modified from Will Dwinnell's original.

    I = Entropy(X,W,lambda,nsymbols,1) + ...
        Entropy(Y,W,lambda,nsymbols,1) - ...
        JointEntropy([X Y],W,lambda,nsymbols,2);
%     I = sum(Entropy([X Y],W,lambda,nsymbols,1)) - ...
%         JointEntropy([X Y],W,lambda,nsymbols,2);

end


function H = Entropy(X,W,lambda,nsymbols,power)

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
    loq = lambda/(nsymbols^power);
    	
    % Calculate sample frequencies
    for symbol = 1:local_symbols
        Count_vector = (X(:,Column) == Alphabet(symbol)).*W;
        Count(symbol) = loq + sum(Count_vector);
    end
	
    % Calculate sample class probabilities
    l_Meff = lambda + sum(Count);
    P = Count / l_Meff;

    % Calculate entropy in bits
    % Note: floating point underflow is not an issue since we are
    % dealing only with the observed alphabet.
    
    H(Column) = -sum(P .* log2(P));
    
end

end


function H = JointEntropy(X,W,lambda,nsymbols,power)
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
H = Entropy(VectorX,W,lambda,nsymbols,power);
end


function H = JointEntropy2(X,W)

% Much more transparent than the previous one, but also much slower.

% Establish size of data
 [nrows,ncols] = size(X);

% Housekeeping

    % Assemble observed alphabet
    Alphabet = unique(X,'rows');
	
    % Housekeeping
    nsymbols = length(Alphabet);
    Count = zeros(nrows,ncols);
    Frequency = zeros(nsymbols,1);
    	
    % Calculate sample frequencies
    for symbol = 1:length(Alphabet)
        for col = 1:ncols
        Count(:,col) = (X(:,col) == Alphabet(symbol,col));
        end
        Count_sum = sum(Count,2);
        Count_ind = Count_sum == ncols;
        Count_ind = Count_ind.*W;
        Frequency(symbol) = sum(Count_ind);
    end
	
    % Calculate sample class probabilities
    P = Frequency / sum(Frequency);
	
    % Calculate entropy in bits
    % Note: floating point underflow is not an issue since we are
    % dealing only with the observed alphabet
    H = -sum(P .* log2(P));

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



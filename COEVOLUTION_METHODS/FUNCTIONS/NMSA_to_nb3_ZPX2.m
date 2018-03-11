function [ZPX2,MIP,MI,ordered_ind,W,ordered_W,ordered,bin_ordered,z_ordered] = ...
    NMSA_to_nb3_ZPX2(nmsa,dist_method,sort_method,...
    threshold,lambda,nsymbols,gaps_method)
% This is the most complete of the 'nb' functions. It implements a
% weighing scheme to account for the similarity between sequences and also
% uses a psudocount ('lambda'). If threshold = 1, all sequences receive a
% weight of 1, and if lambda is set to zero no pseudocount is applied. This
% function can also use a completely automatic way (sort_method =
% 'WEIGHTS') of determining weights based on the distance matrix and
% reorders the msa based on those weights. In this case the value of
% threshold is ignored. In this function there is no default syntax and all
% the input parameters must be entered. 'nsymbols' (20 or 21) is only a
% divider that scales the lambda of the pseudocount. It is difficult to
% recommend a best syntax for every sequence: simulations can be very
% helpful for this purpose. However, it is not recommended to use together
% 'DISTANCE' with a similarity threshold < 1.0. Possible types of syntax
% are: 
% [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,'GAPS','WEIGHTS',1.0,0.5,21,'INCLUDE'); 
% [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,'NOGAPS','EIGEN',0.6,1.0,20,'EXCLUDE'); 
% [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,'NOGAPS','DISTANCE',1.0,0.5,21,'INCLUDE'); 
% [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,'GAPS','DISTANCE',1.0,0.5,21,'EXCLUDE');
% It becomes essentially identical to NMSA_to_nbZPX2 if run as: 
% [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,'GAPS','DISTANCE',1.0,0.0,21,'INCLUDE');
%
% If unsure about what syntax to use probably the most consistent results
% are obtained with:
% [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,'NOGAPS','DISTANCE',1.0,0.5,21,'EXCLUDE'); 

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

    % Here we find which sequence will be at the top in the reordered msa.
        triu_dist = triu(mdist,1);
        [~,max_ind] = max(triu_dist(:));
        [firstrow,firstcol] = ind2sub([nrows,nrows],max_ind);

    % Choose as 1st the sequence with the smallest number of gaps, if the
    % two sequences are of different length.
    
        if dist(firstcol,firstcol) > dist(firstrow,firstrow)
            firstrow = firstcol;
        end
        
    % Now that we have selected the 1st sequence, the merged distance
    % matrix becomes the final distance matrix.
    
        dist = mdist;
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;

    % Here we find which sequence will be at the top in the reordered msa.
        triu_dist = triu(dist,1);
        [~,max_ind] = max(triu_dist(:));
        [firstrow,~] = ind2sub([nrows,nrows],max_ind);
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

%--------------------------------------------------------------------------
% Here we reorder the indices of the msa and of the weights vector. We have
% three options: we can use the distance matrix as such, its eigen
% decomposition, or the columns sums.

switch sort_method
%--------------------------------------------------------------------------
    case 'DISTANCE'

        % The index of the 1st sequence is the one we found before.
        ordered_ind = zeros(nrows,1);
        ordered_ind(1) = firstrow;
        
        ind = firstrow;

        for i = 2:nrows
            % zero the entire test column
            dist(:,ind) = 0;
            % find the index of the closest sequence
            [~,newind] = max(dist(ind,:));
            ordered_ind(i) = newind;
            % the next index will be the index of the sequence just selected
            ind = newind;
        end        
%--------------------------------------------------------------------------
    case 'EIGEN'
    % Here we reorder based on the eigenvector corresponding to the  
    % largest absolute eigenvalue of the distance matrix.
    % [V,~] = eigs(dist,1,'lm');
    [V,D]=eig(dist);
    [sorted_D,ind]=sort(abs(diag(D)),'descend');
    top_V = V(:,ind(1:nrows));
    top_D = sorted_D(1:nrows);

    top_D_mat = zeros(nrows);

    for i = 1:nrows
    top_D_mat(i,i) = top_D(i,1);
    end

    for i=1:nrows
    switch_sign = sign(mean(top_V(:,i)));
    top_V(:,i) = switch_sign*top_V(:,i);
    end
    % V = sum(V,2);
    % switch_sign = sign(mean(V));
    % V = switch_sign*V;
    [~,ordered_ind] = sort(abs(top_V(:,1)),'descend');    
%--------------------------------------------------------------------------
    case 'WEIGHTS'
    % Here we reorder based on the sum of the columns of the distance
    % matrix.        
    dist_sum = sum(dist)';
    [~,ordered_ind] = sort(dist_sum,'descend'); 
    W = ncols./dist_sum;
%--------------------------------------------------------------------------    
end

% Here we reorder the msa and weights vector based on the reordered
% indices.

    ordered = nmsa(ordered_ind,:);
    ordered_W = W(ordered_ind,1);

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

% Here we convert to 0 all the positions that have a value of 21.
% This is equivalent to assuming the each time there is a gap in the real
% msa, the position does not change in the reordered msa.

switch gaps_method
    case 'EXCLUDE'
        ind_21 = z_ordered == 21;
        z_ordered(ind_21) = 0;
        post_binsum = numel(find(z_ordered));
        fprintf('Gap Corrected Binary Sum = %f \n', post_binsum);        
    case 'INCLUDE'
end

    % MI matrix
    [MI] = NMSA_to_MI(z_ordered,ordered_W,lambda,nsymbols);
    % MIP matrix
    MIP = MI_to_MIP(MI,ncols);
    % ZPX2 matrix
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
                
end


function [ MI_mat ] = NMSA_to_MI(nmsa,W,lambda,nsymbols)

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

    % Calculate sample class probabilities: hard to see which of the two
    % lines below is more correct theoretically. Using 'length(X)' seems to
    % require more uniform values of the QUIC subroutine 'lambda'.
    % l_Meff = lambda + sum(Count);
    l_Meff = lambda + length(X);    
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



function [pcorr_nulsp,e_nmsa,W,pcorr_colsp,pcorr_vsum,pcorr_orig] = ...
    NMSA_to_PCORR_property(nmsa,property,dist_method,threshold,...
    psc_lambda,nsymbols,slice_scale,col_scale,proj_space,nprocs)
% This function calculates a matrix of partial correlation coefficients for
% the columns of the msa after it has been converted to a non-integer
% 'property'format: this can be a frequency or property based
% representation of each residue. The idea is that for every pair of
% columns, we are going to get the columns projection on the subspace of
% the matrix that does not contain those two columns. 'proj_space' can be
% 'ALL' for the column and left-null spaces or 'NULL' for only the
% left-null space. 'slice_scale' determines whether the columns are divided
% by the number of observations. This number can be the sum of the values
% in each slice ('SUM'), the total number of rows ('NROWS'), the sum of all
% the row weights ('WEIGHTS'), or 1 ('NONE'). In principle, if we want to
% obtain the exact value of the correlation/partial correlation,
% 'slice_scale' should be set to 'NONE'. 'col_scale' ('GAPS'|'NONE')
% determines whether each slice is scaled by the percentage of non-gaps
% rows.

%--------------------------------------------------------------------------

[nrows,ncols] = size(nmsa);
ncols_o = ncols;
e_nmsa = zeros(nrows,ncols);

% First we replace unusual symbols possibly present in the nmsa

ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

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
        dist = mdist + mdist' + eye(nrows);
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;
end

% Here we calculate the weights for the similarity between sequences.

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff=round(sum(W));
end

fprintf('Meff = %d \n', Meff);

% Here we create a 3d matrix in which every layer has the dimensions of the
% nmsa and represents one of 20 symbols (if gaps are excluded). Then every
% row of each layer is scaled by the weight of that sequence.

nmsa_3 = zeros(nrows,ncols,nsymbols);
for i = 1:nsymbols    
    nmsa_3(:,:,i) = nmsa == i;        
end

% Here we add a total extra fraction of 1 aa equal to the pseudocount
% value to the columns of each slice that count 0; the fraction is added
% only to the rows which are not gaps. Then every row of each slice is
% scaled by the weight of that sequence, based on the similarity threshold.
% Finally we determine the mean of each column counting only the rows
% in which there are no gaps, and we subtract this mean from each column.

for i = 1:ncols
    slice = squeeze(nmsa_3(:,i,:));
    slice_zeros = slice == 0;
    row_gap_ind = nmsa(:,i) == 21;
    row_nongap_ind = nmsa(:,i) ~= 21;
    slice_zeros(row_gap_ind,:) = false;
    col_sum = sum(slice,1);
    col_nongap_ind = col_sum ~= 0;
    slice_zeros(:,col_nongap_ind) = false;
    loq = psc_lambda/sum(slice_zeros(:));
    slice(slice_zeros) = loq;
    slice = slice.*W(:,ones(1,nsymbols));
    % The following line replaces dividing by the no. of observations.
    switch slice_scale
        case 'SUM'
        slice = slice/sum(slice(:));
        case 'NROWS'
            slice = slice/nrows;
        case 'WEIGHTS'
            slice = slice/sum(W);
        case 'NONE'
    end    
    slice_mean = mean(slice(row_nongap_ind,:));
    slice_mean_mat = zeros(nrows,nsymbols);
    row_non_gap_ind = find(row_nongap_ind);
    for k = 1:length(row_non_gap_ind)
        k_ind = row_non_gap_ind(k);
        slice_mean_mat(k_ind,:) = slice_mean;
    end
    slice = slice - slice_mean_mat;
    
    % Here we scale each slice by the number of gaps in the original column
    % of the msa.
    sumW = sum(W);
    if strcmp(col_scale,'GAPS')
        slice = slice*(sum(row_nongap_ind.*W)/sumW);
    end        
    nmsa_3(:,i,:) = slice;

    if size(property,2) > 1
    e_nmsa(:,i) = slice*property(:,i);
    else    
    e_nmsa(:,i) = slice*property;
    end
end

%--------------------------------------------------------------------------

% Set the pool of cpu's for parallel processing
        
        if matlabpool('size') > 0
            matlabpool close
        end
        
        if nprocs > 1
        matlabpool(nprocs)
        end

% Here we calculate directly the correlation and/or partial correlation 
% matrix.
        
pcorr_nulsp = zeros(ncols_o,ncols_o);
pcorr_colsp = zeros(ncols_o,ncols_o);
pcorr_vsum = zeros(ncols_o,ncols_o);

switch proj_space
    case 'ALL'

        % Here we establish the total number of cycles required to complete 
        % the calculation: this is used to determine the time left to
        % completion.
        total_count = ncols_o*(ncols_o - 1)/2;
        cycle_count = 0;
        time = 0;
        
        row = true(1,ncols_o);
                                
        for i = 1:ncols
            
            tic
            v1 = e_nmsa(:,i);
            
            parfor j = i:ncols_o
                
                v2 = e_nmsa(:,j);                
                row1 = row;

                % Here we select all the columns except the two for which
                % we want to measure the correlation coefficient.
                row1([i j]) = false;
                s_nmsa = e_nmsa(:,row1);

                % We use 'svd' to identify the 'column space' and the 'left
                % null space' of the matrix that does not contain the two
                % columns v1 and v2.
                
                [U,S,~] = svd(s_nmsa);

                % Alternatively we could take the transpose first and then
                % work on V rather than U. In both cases, the column/row
                % space bases are exactly the same, but the null/left null
                % space bases come out different (although both perfectly
                % orthonormal).
                
                % [~,S,V] = svd(s_nmsa');
                               
                % We determine the rank 'r'.
                s = diag(S);
                tol = max(size(s_nmsa))*eps(max(s));
                r = sum(s > tol);
                                
                % Here we select only the first 'r' columns of U (or V if 
                % we took the transpose).
                C = U(:,1:r);

                P1 = C*C';
                cv1 = P1*v1;
                cv2 = P1*v2;
                pcorr_colsp(i,j) = corr(cv1,cv2);
                % pcorr_colsp(i,j) = (cv1'*cv2);
                
                % In principle we should calculate all the dot
                % products after multiplying each vector by the projection
                % matrix P2:
                %   cov = (P1*v1)'*(P1*v2). 
                % If we take the transpose of the 1st term in the product,
                % that becomes:
                %   cov = v1'*(P1'*P1)*v2
                % Since the square of the projection matrix is itself we
                % can simplify further and then extract the norm directly:

                % pcorr_colsp(i,j) = (v1'*P1*v2);

                % We get also the perpendicular components: these are the
                % projection vectors onto N(s_nmsa'). The projection matrix
                % is I-P1.
                % [p_nrows,p_ncols] = size(P1);
                % P2 = eye(p_nrows,p_ncols);
                % P2 = P2-P1;

                % Alternatively we could use the columns from U (or V if we
                % initially took the transpose) that correspond to
                % N(s_nmsa').
                
                N = U(:,r+1:end);
                P2 = N*N';
                nv1 = P2*v1;
                nv2 = P2*v2;
                pcorr_nulsp(i,j) = corr(nv1,nv2);
                % pcorr_nulsp(i,j) = nv1'*nv2;

                % Same as the following line:
                % pcorr_nulsp(i,j) = (v1'*P2*v2);

                % Finally we check that the two perpendicular projections
                % sum up to the original vector.
                rv1 = nv1 + cv1;
                rv2 = nv2 + cv2;
                pcorr_vsum(i,j) = corr(rv1,rv2);
                % pcorr_vsum(i,j) = rv1'*rv2;
 
            end
            
            cycle_time = toc;
            time = time + cycle_time;  
            cycle_count = cycle_count + (ncols_o - i);
            progress = cycle_count/total_count;
            count_left = total_count - cycle_count;
            fprintf('Percent completion = %5.3f \n', progress);
            time_left = count_left*(time/cycle_count)/60;
            fprintf('Time elapsed: %6.1f minutes \n',time/60);            
            fprintf('Est. time to completion: %6.1f minutes \n',time_left);
                       
        end
        
        for i = 1:ncols_o            
            for j = i:ncols_o
                pcorr_colsp(j,i) = pcorr_colsp(i,j);
                pcorr_nulsp(j,i) = pcorr_nulsp(i,j);
                pcorr_vsum(j,i) = pcorr_vsum(i,j);
            end
        end

    pcorr_orig = corr(e_nmsa);
    % pcorr_orig = e_nmsa'*e_nmsa;
        

    case 'NULL'
        
        % Here we establish the total number of cycles required to complete 
        % the calculation: this is used to determine the time left to
        % completion.
        total_count = ncols_o*(ncols_o - 1)/2;
        cycle_count = 0;
        time = 0;
        
        row = true(1,ncols_o);
        
        for i = 1:ncols_o
            
            tic
            v1 = e_nmsa(:,i);
                        
           parfor j = i:ncols_o
                
                v2 = e_nmsa(:,j);
                row1 = row;
        
                % Here we select all the columns except the two for which
                % we want to measure the correlation coefficient.
                
                row1([i j]) = false;
                s_nmsa = e_nmsa(:,row1);
                                
                % We use 'svd' to identify the 'column space' and the 'left
                % null space' of the matrix that does not contain the two
                % columns v1 and v2.
               
                [U,S,~] = svd(s_nmsa);

                % We determine the rank 'r'.
                s = diag(S);
                tol = max(size(s_nmsa))*eps(max(s));
                r = sum(s > tol);
                
                % Here we calculate the projection vectors of the two
                % selected vectors from nmsa onto N(s_nmsa').
                                
                N = U(:,r+1:end);              
                P2 = N*N';
                                 
                % Here we calculate a cross-covariance matrix slice by
                % slice. In principle we should calculate all the dot
                % products after multiplying each vector by the projection
                % matrix P2:
                %   cov = (P2*i_slice)'*(P2*j_slice). 
                % If we take the transpose of the 1st term in the product,
                % that becomes:
                %   cov = i_slice'*(P2'*P2)*j_slice
                % Since the square of the projection matrix is itself we
                % can simplify further and then extract the norm directly:
                % pcorr_nulsp(i,j) = (v1'*P2*v2);

                nv1 = P2*v1;
                nv2 = P2*v2;
                pcorr_nulsp(i,j) = corr(nv1,nv2);
                                                                      
           end
           
            cycle_time = toc;
            time = time + cycle_time;  
            cycle_count = cycle_count + (ncols_o - i);
            progress = cycle_count/total_count;
            count_left = total_count - cycle_count;
            fprintf('Percent completion = %5.3f \n', progress);
            time_left = count_left*(time/cycle_count)/60;
            fprintf('Time elapsed: %6.1f minutes \n',time/60);            
            fprintf('Est. time to completion: %6.1f minutes \n',time_left);
            
        end
        
        for i = 1:ncols_o            
            for j = i:ncols_o
                pcorr_nulsp(j,i) = pcorr_nulsp(i,j);
            end
        end

end

    if nprocs > 1
        matlabpool close
    end
    
end





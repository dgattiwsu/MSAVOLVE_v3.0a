function [recomb_MSA,history,recomb_history,COV,recomb_COV,glob_COV] = ...
    hor_transfer_6s_3D_cov_hist( MSA,history,recomb_history,COV,...
    recomb_COV,glob_COV,hrf,nzones,nhrt )

% This function produces a recombination of a MSA based on:
% nzones = the number of fragments of a sequence that can be exchanged
% hrf = the sequence numbers corresponding to the edges of the
% recombination zones.
% nhrt = the number of times the recombination is carried out.
% In this algorithm the recombination is rather seen as a wave that extends
% a certain fragment to a number of sequences.

% Please, REMEMBER! the 'nseq' inside this function is not the same value
% as 'nseq' in the main program. Here 'nseq' has a different value in each 
% of the three levels.
[nseq,npos] = size(MSA);

% Note that here 'nedges' is nzones + 1 (not nzones - 1) because we are
% including also the 1st and last position of the MSA.
nedges = nzones + 1;

% Here we first store a copy of the starting history and COV. They will be
% used to recover the contribution of just the recombination process to the
% history and COV matrices.

MSA1 = MSA;

for i=1:nhrt
    MSA0 = MSA1;    
    % Here we select the number of rows that will become identical in a
    % certain fragment.
    n = randi(nseq);
    % Here we pick the upstream and downstream edge of the recombining
    % fragment
    ind = randi(nedges,1,2);
    ind = sort(ind);
    % Here we pick which rows (sequences) will be affected by the
    % recombination wave.
    ind2 = randi(nseq,n,1);
    % Here we find the sequence numbers corresponding to the upstream and 
    % downstream edge of the recombining fragment.
    a = hrf(ind(1));
    b = hrf(ind(2));
    % Here in all the rows affected by the wave the recombining fragment
    % becomes identical to that of one of the rows picked at random.
    r = randi(n);
    for j = 1:n    
    MSA1(ind2(j),a:b) = MSA1(ind2(r),a:b);
    end

    history_part = MSA1 ~= MSA0;
    % CRITICAL: make sure history arrays are numerical.
    history_part = + history_part;
    
    for k = 1:nseq
        COV_part = zeros(npos,npos,nseq);
        for m = 1:npos
            for j = m:npos    
                COV_part(m,j,k) = history_part(k,m)*history_part(k,j);
            end
        end
        COV = COV + COV_part;    
        glob_COV = glob_COV + COV_part;    
        recomb_COV = recomb_COV + COV_part;    
    end
    
% CRITICAL: make sure history arrays are numerical.    
history = + history;
recomb_history = + recomb_history;
history = history + history_part;
recomb_history = recomb_history + history_part;
end

recomb_MSA = MSA1;

end


% Copyright (c) 2012, Sharon H. Ackerman, Elisabeth Tillier, Domenico L. Gatti
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
% 
% MSAvolve v.1.0: This program simulates the evolution of a protein from 
% a hypothetical ancestor, and returns a multiple sequence alignment (MSA).
% 
%% 'run' section.
% Here we decide if it is a completely new run or a restart from a crush.
% This is a supercell and will run everything down to the assignment of all
% the array. If you are running in cell mode, skip this cell and start from
% the 'random number generator' section. 
if NEW_RUN
    start_cycle = 1;
    
%% 'random number generator' section.

% Here we initialize the random number generator such that each time we get a 
% different result. 

if isempty(RANDOM_SEED_FILE) || isempty(RANDOM_SEED_VAR)
    s = rng('shuffle');
else 
% Here we recall the stored startpoint for the random number generator to
% compare two runs exactly.
    load(RANDOM_SEED_FILE,RANDOM_SEED_VAR); 
    rng(eval(RANDOM_SEED_VAR));
end

 %% 'path' section
 
 % Here we can add some paths to m-files that will be heavily used. 
 % For example:
 
%  addpath('DWINNEL_MI');
%  addpath('MATLAB_LOCAL');
 
 %% 'Reference msa' section.

% If we use a reference MSA we first read a MSA in Clustal or MSF format and 
% convert it to the standard Matlab numeric format.
% REF_smsa is a structure containing both headers and sequences.
% REF_nmsa is the MSA in numeric format.
% REF_cmsa is a form of the letter format without headers

if MAKE_MSA

% Make sure the total number of sequences can be evenly divided by the 
% number of branches    
branch_no_seq = round(no_seq/no_branches);
no_seq = branch_no_seq*no_branches;

% Make sure the largest crossover points does not exceed the number of 
% positions.
if max(crossover_points)>no_pos
disp('ERROR: MSAvolve cannot continue because')    
fprintf('crossover point %d is larger \n than the number of positions. \n',...
    max(crossover_points));
return
end    

for i = 1:no_seq
    for j = 1:no_pos
    REF_nmsa(i,j) = randi(20);
    end
end
clear no_pos

REF_cmsa = int2aa(REF_nmsa);

for i = 1:size(REF_nmsa,1)
    REF_smsa(1,i).Header = ['Sequence ' num2str(i)];
    REF_smsa(1,i).Sequence = int2aa(REF_nmsa(i,:));    
end

else

    if fasta
        [REF_smsa,REF_nmsa] = faln_to_nmsa(ext_MSA);
    else
        [REF_smsa,REF_nmsa] = aln_to_nmsa(ext_MSA);
    end
        
REF_cmsa = int2aa(REF_nmsa);

end

[nseq,npos] = size(REF_nmsa);

% Then we need to set a table for the background and for a "biased"
% occurrence of aa's based on the profile of the experimnetal MSA. 
% We can use values from the entire database or just
% from a set of proteins with which we are familiar. In this latter case it
% is convenient to take the 'null' and 'match' emissions from a hmm 
% model. We can derive the emissions using the external program HMMER v3.0
% or we can calculate them using the BioInformatics Toolbox of Matlab:

REF_hmm_model = nmsa_to_hmm_model(REF_nmsa);
bg_frequencies = REF_hmm_model.NullEmission';
t_emissions = REF_hmm_model.MatchEmission';

% Alternatively we could use established background probabilities derived
% from a non redundant database of proteins. For convenience we store the 
% values as a column vector. The frequencies below refer to the following
% symbols: A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V

% bg_frequencies = [ 0.0730 0.0520 0.0430 0.0500 0.0250 0.0400 0.0610 ...
%    0.0720 0.0230 0.0530 0.0890 0.0640 0.0230 0.0420 0.0520 0.0730 ...
%    0.0560 0.0130 0.0330 0.0630 ]';

% Next we calculate the profile of the MSA including gaps: the 21th row of 
% the profile refers to the gaps which in matlab are however represented by
% number 25. 

[REF_profile_gaps,REF_symbols_gaps] = seqprofile(REF_smsa,'Gaps',gaps_count);
REF_bg_frequencies_gaps = sum(REF_profile_gaps,2)/npos;

% Here we get a measure of relative entropy taking into considerations the
% presence of gaps in the reference msa.
% [ rH_ai,rH_i,D_ai,D_i,Pm_ai ] = ...
%    kldiv_gaps_nmsa( REF_nmsa,REF_bg_frequencies_gaps,REF_profile_gaps );

% Next we calculate the profile of the MSA without including gaps.
[REF_profile,REF_symbols] = seqprofile(REF_smsa,'Gaps','none');

% Here we determine the level of conservation using the traditional 
% definitions of entropy and relative entropy.

REF_entropy = Entropy(REF_nmsa);
[REF_rel_entropy_aa,REF_rel_entropy] = rel_entropy_nmsa(REF_nmsa,...
    bg_frequencies,REF_profile);

% Here we find the indices of the positions that are fully conserved. These
% positions must be removed from any random assignment of covarying
% positions.

REF_conserved_ind = find(REF_entropy == 0); 

% Here we calculate the relative entropy using as profile the emission
% probabilities from the hmm model.

[REF_hmm_rel_entropy_aa,REF_hmm_rel_entropy] = rel_entropy_nmsa(REF_nmsa,...
    bg_frequencies,t_emissions);

% Unflag the following part to plot the relative entropy values.
%
% REF_seq_cons=figure; clf; 
% set(REF_seq_cons,'Units','normalized','Position',[0 0.5 0.4 0.5],'Name',...
%     'REF Sequence Conservation');
%
% subplot(2,1,1);
% imagesc(REF_rel_entropy_aa,'DisplayName','REF_rel_entropy_aa');figure(gcf)
% xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
% ylabel('REF AAs Relative Entropy','FontSize',12,'FontWeight','n');
% set(gca,'Xlim',[1 npos],'YTickLabel',{'A','R','N','D','C','Q',...
%     'E','G','H','I','L','K',...
%     'M','F','P','S','T','W','Y','V'},...
%     'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20],...
%     'YDir','reverse',...
%     'TickDir','in',...
%     'Layer','top');
%
% subplot(2,1,2);
% stairs(REF_rel_entropy,'-r','DisplayName','REF_rel_entropy_sum');figure(gcf)
% xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
% ylabel('REF Total Relative Entropy','FontSize',14,'FontWeight','n'); 
% set(gca,'Xlim',[1 npos]);

% An alternative way to measure conservation is to calculate
% the consensus sequence of the experimental msa and the
% consensus score for each position. The scores are the 
% average euclidean distance between the scored symbol and the 
% M-dimensional consensus value. Therefore, the higher the level 
% of conservation, the lower the consensus score. For example for a fully 
% conserved position the consensus score is 0. 

[REF_s_consensus,REF_s_consensus_score] = seqconsensus(REF_cmsa);

% Here we plot together the two measures of conservation.

% REF_seq_cons_2=figure; clf; 
% set(REF_seq_cons,'Units','normalized','Position',[0 0.2 0.5 0.8],'Name',...
%    'REF Sequence Conservation');
%
% subplot(2,1,1);
% stairs(REF_rel_entropy,'-r','DisplayName','REF_rel_entropy_sum');figure(gcf)
% xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
% ylabel('REF Total Relative Entropy','FontSize',14,'FontWeight','n'); 
% set(gca,'Xlim',[1 npos]);
%
% subplot(2,1,2);
% stairs(REF_s_consensus_score,'-b','DisplayName',...
%     'REF_consensus_score');figure(gcf)
% xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
% ylabel('REF Consensus Score','FontSize',14,'FontWeight','n'); 
% set(gca,'Xlim',[1 npos],'Ylim',[-1 (max(REF_s_consensus_score)+1)]);

if GAPS
% Here we assign a mutation probability at each position of the sequence: 
% we get the probabilities from the profile of the reference msa calculated
% including gaps.

symbols = ([1:20,25])';
mut_prob_PD = cell(npos,1);
mut_freq_PD = zeros(21,npos);

for i = 1:npos
    mut_prob_PD{i} = fitdist(symbols,'kernel','width',0.01,'frequency',...
        round(REF_profile_gaps(:,i)*1e6));
end

for i= 1:npos
nfreq = size(mut_prob_PD{i,1}.InputData.freq,1);
mut_freq_PD(1:nfreq,i) = mut_prob_PD{i,1}.InputData.freq;
end

REF_mut_prob_PD = mut_prob_PD;
REF_mut_freq_PD = mut_freq_PD;


bg_prob_PD = fitdist(symbols,'kernel','width',0.01,'frequency',...
         round(REF_bg_frequencies_gaps*1e6));

clear symbols i

else
    
% Here we get the probabilities from the match emissions of the hmm model.

symbols = (1:1:20)';
mut_prob_PD = cell(npos,1);
mut_freq_PD = zeros(20,npos);

for i = 1:npos
    mut_prob_PD{i} = fitdist(symbols,'kernel','width',0.01,'frequency',...
        round(t_emissions(:,i)*1e6));
end

for i= 1:npos   
mut_freq_PD(:,i) = mut_prob_PD{i,1}.InputData.freq;
end

REF_mut_prob_PD = mut_prob_PD;
REF_mut_freq_PD = mut_freq_PD;

bg_prob_PD = fitdist(symbols,'kernel','width',0.01,'frequency',...
         round(bg_frequencies*1e6));

clear symbols i
    
end

    
%% 'ancestor preparation' section.

% The ancestor is actually an m (sequences) by n (aa's) array of sequences 
% all identical to the ancestor. These sequences will become the main branches 
% of the phylogenetic tree. If we are using a reference msa, we can determine 
% the number of branches in that msa from the distance matrix or from the 
% covariance matrix. We calculate the distance between all the sequences 
% with the p-distance method, where p is the proportion of sites at which the 
% two sequences differ. p is close to 1 for poorly related sequences, and is 
% close to 0 for similar sequences. We can calculate a covariance matrix of the msa after 
% converting the nmsa into a binary form in which each sequence of length n 
% is treated as a vector of size 20*n with 0 and 1. 
% VERY IMPORTANT: in this case the variables are considered the various
% sequences and the observations are the different positions in each
% sequence. Therefore, the covariance matrix that we will obtain will be 
% a 'sequence correlation' matrix.
% For this reason, before calculating the covariance matrix we take the 
% transpose of the binary matrix. For example, if the binary matrix is of
% dimensions [348,7000]: its transpose will be [7000 obs , 348 vars] and
% the covariance matrix will be [348,348] matching the distance matrix. 
% If we were not taking the transpose, the different positions in each 
% sequence would be the 'variables' and the different sequences would be 
% the 'observations'. In this case the covariance matrix of dimensions
% [7000,7000] would be a 'positional correlation' matrix.
 
 if DISTANCE
 REF_dist = seqpdist(REF_cmsa,'SquareForm',true,'Method','p-distance');
 REF_cov_mat = 1-REF_dist;
 else
 REF_binmsa = nmsa_to_binmsa(REF_nmsa);
 REF_cov_mat = cov(REF_binmsa',1);  
 end

% Although not used here immediately the covariance matrix provides a
% way to define an overall similarity score between all the sequences using
% the definition of multidimensional mutual information from:
% Wang,B. and Shen,Y. "A method on calculating high-dimensional mutual 
% information and its application to registration of multiple ultrasound 
% images" Ultrasonics 44 (2006) e79?e83.
% Here we used the covariance matrix instead of the mutual information
% matrix to express the similarity between sequences.
 
 REF_sim_score = multidim_MI(REF_cov_mat);

% If we want to calculate manually the covariance matrix without using the
% 'cov' function from Matlab Statistics ToolBox we can use the following:

% data = REF_binmsa';
% [obs,vars] = size(data);
% for i = 1:vars
%   centered_data(:,i) = data(:,i) - mean(data(:,i));
% end
% cov_data = centered_data' * centered_data;
% cov_data = cov_data/obs; 

% Next we carry out a spectral analysis of the covariance matrix:

[REF_pc,REF_ev] = spectral(REF_cov_mat);

% IMPORTANT: the columns in the eigevector matrix are the principal 
% components. The elements (rows) of the eigenvectors are the principal
% components 'coefficients' or 'loadings'. The meaning is the following:
% the rows of the eigenvectors are the contributions or 'coefficients' of 
% each variable (each sequence) to the principal component. In our case
% the variables are the individual sequences, and so the coefficents of the
% eigenvectors of the covariance matrix represent the contribution of each
% sequence to the direction of the principal components in the nseq-space.
% We also notice that according to this formalism the individual aa's in
% each sequence represent the 'observations'.
% Therefore a scatter plot of the first two principal components reveals if
% there are clusters of sequences. A similar result would be obtained 
% using the distance matrix instead of the covariance matrix; in fact the 
% distance matrix is also positive-semidefinite and therefore is also 
% a form of covariance matrix. Typically, without even plotting we can get 
% the clusters from the 1st two PC's.

REF_Clusters = clusterdata(REF_pc(:,1:DIMENSIONS),'maxclust',MAX_CLUSTERS);

% Here we override the cluster determination using only a phylogenetic tree.

if TREE_ONLY
    smsa_dist = seqpdist(REF_smsa,'Method',distance_method);
    if tree_method == 'seqneighjoin'
        smsa_tree = seqneighjoin(smsa_dist,neighjoin_method);
    else
        smsa_tree = seqlinkage(smsa_dist,tree_method);
    end
    [REF_Clusters,NODE_Clusters,~] = ...
     cluster(smsa_tree,[],'criterion','silhouette','MaxClust',MAX_CLUSTERS);

    h = plot(smsa_tree,'Orientation','top','TerminalLabels','false'); 
    set(h.BranchLines(NODE_Clusters==1),'Color','b')
    set(h.BranchLines(NODE_Clusters==2),'Color','r')
    set(h.BranchLines(NODE_Clusters==3),'Color','g')
    set(h.BranchLines(NODE_Clusters==4),'Color','c')
    set(h.BranchLines(NODE_Clusters==5),'Color','y')
    set(h.BranchLines(NODE_Clusters==6),'Color','k')
    set(h.BranchLines(NODE_Clusters==7),'Color',[1.0,0.5,0.0])
    set(h.BranchLines(NODE_Clusters==8),'Color','m')
    set(h.BranchLines(NODE_Clusters>8),'Color',[0.3,0.6,0.2])
end

% The number of main branches in the simulated tree are set accordingly.

nbranches = max(REF_Clusters);

% If we are making our own reference MSA we can override the number of
% branches assigned automatically.

if MAKE_MSA
    branch_no_seq = round(no_seq/no_branches);
    nbranches = no_branches;
    REF_Clusters = [];
    for i = 1:nbranches
        temp_clusters = ones(branch_no_seq,1)*i;
        REF_Clusters = [REF_Clusters;temp_clusters];
    end
clear branch_no_seq temp_clusters  
end

% REF_cov_mat_lin=nonzeros(triu(REF_cov_mat,1));

% We will re-use this figure later comparing the results with the same
% analysis carried out with the simulated MSA. Thus, no need to plot it now
% unless you want to check the assignment of the main branches.

%     if CHECK_BRANCHES
%         [~] = make_REF_MSA_comp_1(REF_nmsa,DISTANCE,DIMENSIONS,...
%         MAX_CLUSTERS,REF_rel_entropy,cov_coord,TREE_ONLY,...
%         distance_method,tree_method,neighjoin_method);    
%     % Here we pass control to the keyboard. If you are satisfied with the
%     % automatic assignment of branches based on clusters, type 'return' and 
%     % then enter. If you want to quit the program type 'dbquit' and then
%     % enter.
%     keyboard
%     end

% We can find out which sequences are in the clusters identified by the
% spectral analysis.

REF_clust1 = find(REF_Clusters == 1);
REF_clust2 = find(REF_Clusters == 2);
REF_clust3 = find(REF_Clusters == 3);
REF_clust4 = find(REF_Clusters == 4);
REF_clust5 = find(REF_Clusters == 5);
REF_clust6 = find(REF_Clusters == 6);
REF_clust7 = find(REF_Clusters == 7);
REF_clust8 = find(REF_Clusters == 8);

% We could also get the msa's corresponding to the clusters...
% REF_msa_clust1 = REF_msa(REF_clust1);
% REF_msa_clust2 = REF_msa(REF_clust2);
% for i = 1:size(REF_msa_clust1,2)
%     REF_msa_clust1(1,i).Header = REF_msa_clust1(1,i).Header(8:end);
% end
% for i = 1:size(REF_msa_clust2,2)
%     REF_msa_clust2(1,i).Header = REF_msa_clust2(1,i).Header(8:end);
% end

% ... and print them out in clustal or msf format.
% multialignwrite('REF_cluster_1.aln',REF_msa_clust1,...
%     'Header','REF Cluster 1','WriteCount','true');
% multialignwrite('REF_cluster_2.aln',REF_msa_clust2,...
%     'Header','REF Cluster 2','WriteCount','true');
% multialignwrite('REF_cluster_1.msf',REF_msa_clust1,...
%     'Header','REF Cluster 1','WriteCount','true');
% multialignwrite('REF_cluster_2.msf',REF_msa_clust2,...
%     'Header','REF Cluster 2','WriteCount','true');

% If the histogram shows a peak at the highest values of
% covariance it means there is a large number of sequences that are too
% similar to each other. We could use a covariance cutoff to remove these
% sequences, but we recommend using the program 't_coffee' to eliminate the
% redundancy from the original data set. For example to keep only sequences
% that are no more identical to any other sequence than 90% use the
% following:
% "t_coffee -other_pg seq_reformat -in complete.aln -action +trim 
%   _aln_%%90_ +upper > trimmed.aln"
% For t_coffee installation see http://www.tcoffee.org/
% Alternatively you can use the utility 'trim_nmsa_by_threshold.m'.

if GAPS

symbols = ([1:20,25])';
REF_branch_mut_prob_PD = cell(npos,nbranches);
REF_branch_mut_freq_PD = zeros(21,npos,nbranches);
REF_branch_t_emissions = zeros(20,npos,nbranches);
REF_branch_profile = zeros(21,npos,nbranches);
REF_branch_sim_score = zeros(nbranches,1);
REF_branch_sim_score_true = zeros(nbranches,1);
REF_branch_size = zeros(nbranches,1);

for k = 1:nbranches
REF_branch_ind = find(REF_Clusters == k);
% We also measure the size of each branch for selection at the end of the
% evolution.
REF_branch_size(k,1) = size(REF_branch_ind,1);
if REF_branch_size(k,1) == 1
    disp('ERROR: there is only one sequence in one of the clusters.');
    disp('We cannot calculate a hmm model from just one sequence!');
    disp('Rerun reducing the number of clusters and/or dimensions.');
    return
end
REF_branch_profile(:,:,k) = ...
    seqprofile(REF_cmsa(REF_branch_ind,:),'Gaps',gaps_count);
        
    for i = 1:npos
        REF_branch_mut_prob_PD{i,k} = ...
            fitdist(symbols,'kernel','width',0.01,'frequency',...
            round(REF_branch_profile(:,i,k)*1e6));
    end

    for i= 1:npos
        nfreq = size(REF_branch_mut_prob_PD{i,k}.InputData.freq,1);
        REF_branch_mut_freq_PD(1:nfreq,i,k) = ...
            REF_branch_mut_prob_PD{i,k}.InputData.freq;
    end
end    

else
    
symbols = (1:1:20)';
REF_branch_mut_prob_PD = cell(npos,nbranches);
REF_branch_mut_freq_PD = zeros(20,npos,nbranches);
REF_branch_t_emissions = zeros(20,npos,nbranches);
REF_branch_profile = zeros(20,npos,nbranches);
REF_branch_sim_score = zeros(nbranches,1);
REF_branch_sim_score_true = zeros(nbranches,1);
REF_branch_size = zeros(nbranches,1);

for k = 1:nbranches
REF_branch_ind = find(REF_Clusters == k);
% We also measure the size of each branch for selection at the end of the
% evolution.
REF_branch_size(k,1) = size(REF_branch_ind,1);
if REF_branch_size(k,1) == 1
    disp('ERROR: there is only one sequence in one of the clusters.');
    disp('We cannot calculate a hmm model from just one sequence!');
    disp('Rerun reducing the number of clusters and/or dimensions.');
    return
end
REF_branch_hmm_model = nmsa_to_hmm_model(REF_nmsa(REF_branch_ind,:));
REF_branch_t_emissions(:,:,k) = REF_branch_hmm_model.MatchEmission';

% If too many clusters we may have that some aa's are not represented. This
% will produce a mismatch with the arrays that expect 20 aa.

if any(REF_branch_hmm_model.MatchEmission(:) == 0)
emis_zero_ind = (REF_branch_t_emissions(1:20,:,k) == 0);
layer = REF_branch_t_emissions(1:20,:,k);
layer(emis_zero_ind) = eps;
REF_branch_t_emissions(1:20,:,k) = layer;
    disp('CAUTION: there is one cluster in which not all aa are represented.');
    disp('We will continue with a numerical fix, but we recommend that you');
    disp('reduce the number of clusters or try different number of dimensions.');
%     return
end

REF_branch_profile(:,:,k) = ...
    seqprofile(REF_cmsa(REF_branch_ind,:),'Gaps','none');
        
    for i = 1:npos
        REF_branch_mut_prob_PD{i,k} = ...
            fitdist(symbols,'kernel','width',0.01,'frequency',...
            round(REF_branch_t_emissions(:,i,k)*1e6));
    end

    for i= 1:npos
        REF_branch_mut_freq_PD(:,i,k) = ...
            REF_branch_mut_prob_PD{i,k}.InputData.freq;
    end
end

end % end of GAPS

    
if DISTANCE
 REF_branch_dist = seqpdist(REF_cmsa(REF_branch_ind),...
     'SquareForm',true,'Method','p-distance');
 REF_branch_cov_mat = 1-REF_branch_dist;
else
 REF_branch_binmsa = nmsa_to_binmsa(REF_nmsa(REF_branch_ind));
 REF_branch_cov_mat = cov(REF_branch_binmsa',1);
end

% Here we calculate a global similarity score for each branch and we store 
% a back-up copy.

REF_branch_sim_score(k) = multidim_MI(nantozero(REF_branch_cov_mat));
REF_branch_sim_score_true(k) = REF_branch_sim_score(k);

clear REF_branch_ind REF_branch_hmm_model 
clear REF_branch_dist REF_branch_binmsa REF_branch_cov_mat

% Here we introduce a multiplier that will be used to scale the simulated
% msa if we choose uneven branches.

if isempty(msa_size)
    msa_size_mult = 1;
else
    msa_size_mult = msa_size/nseq;
end

if UNEVEN_BRANCHES

% Here we scale the size of the branches such that the final msa has
% approximately the size of the reference msa.

 MSA_branch_size = ncopies_1*ncopies_2;

for i = 1:nbranches
     MSA_branch_mult = ...
         sqrt(msa_size_mult*REF_branch_size(i)/MSA_branch_size);
     bncopies_1(i) = round(ncopies_1*MSA_branch_mult);     
     bncopies_2(i) = round(ncopies_2*MSA_branch_mult);
end
s_MSA_branch_size = (bncopies_1.*bncopies_2)';
REF_MSA_branch_size_ratio = msa_size_mult*REF_branch_size./s_MSA_branch_size; 
nselect = sum(s_MSA_branch_size);

else
    
for i = 1:nbranches    
     bncopies_1(i) = ncopies_1;     
     bncopies_2(i) = ncopies_2;
end    
s_MSA_branch_size = (bncopies_1.*bncopies_2)';
REF_MSA_branch_size_ratio = REF_branch_size./s_MSA_branch_size; 
nselect = sum(s_MSA_branch_size);
 
end

if RECOMB_SCALING 
% Here we scale the branch similarity scores such that the sum of the 
% scores is equal to the number of branches. ...Just trust me on this.

 REF_branch_sim_score = ...
     REF_branch_sim_score_true*nbranches/sum(REF_branch_sim_score_true);

else
    
% Here we set to 1 the score effectively used during the evolution process.
 REF_branch_sim_score = ones(k,1);

end     % End of RECOMB_SCALING

%% 'recombination' section
if isempty(crossover_points)
    hrf = [1 npos];
    nzones = numel(hrf)-1;
    rec_cycles = 0;
    nhrt_range_1 = [0,0];
    nhrt_range_2 = [0,0];
    nhrt_range_3 = [0,0];
    nhrt_range_4 = [0,0];
    nhrt_range_5 = [0,0];
else
hrf = [1 crossover_points npos];
end

if ONE_RAND_RECOMB_ZONES
    RAND_RECOMB_ZONES = 0;
    nedges = nzones - 1;
    hrf = randi(npos,1,nedges);
    hrf = sort(hrf,2);
    hrf = [1 hrf npos];
end

%% 'covarions' section

% This is an extremely tricky section at the core of MSAvolve. There
% is a lot of redundancy in the checks and the code could be slimmer, but
% it works well in its current version: don't change it.

% If SET_COVARIONS is false and RAND_COVARIONS is true this whole section 
% is skipped and covarions are assigned randomly at each evolution trial.

if SET_COVARIONS
    RAND_COVARIONS = 0;

% Here we randomly pick a certain percentage (e.g.10%) of the positions 
% in the sequence of the ancestor to be covarions by calling the function 
% that sets the covarying columns in the msa.

if ENTROPY
    REL_ENTROPY = 0;
    METHOD = 0;
    disp('Entropy based assignment of covarions');
    % [ covar_vec_low,covar_vec_mid,covar_vec_high ] = ...
    % pick_covarions_3(REF_nmsa,npos,fcov,cov_gaps);
    [ covar_vec_low,covar_vec_mid,covar_vec_high ] = ...
    pick_covarions_5(REF_nmsa,REF_entropy,npos,fcov,cov_gaps);
end

if REL_ENTROPY
    ENTROPY = 0;
    METHOD = 0;
    disp('Relative Entropy based assignment of covarions');
    % [ covar_vec_low,covar_vec_mid,covar_vec_high ] = ...
    % pick_covarions_4(REF_nmsa,bg_frequencies,npos,fcov,cov_gaps);
    [ covar_vec_low,covar_vec_mid,covar_vec_high ] = ...
    pick_covarions_5(REF_nmsa,REF_rel_entropy,npos,fcov,cov_gaps);
end

% Choose which range of (rel)entropy will be used: low,mid,high. Remember
% that high entropy corresponds to low relative entropy and low entropy to
% high relative entropy. Thus, your choice depends on which kind of entropy
% is being used to select the covarions. 

if (ENTROPY || REL_ENTROPY)
    switch cov_entropy_range
        case 'low'
            covar_vec = covar_vec_low;
        case 'mid'
            covar_vec = covar_vec_mid;
        case 'high'
            covar_vec = covar_vec_high;
        case 'low-mid'
           merge_covar_vec = [covar_vec_low covar_vec_mid];
           ncov = size(merge_covar_vec,1);
           covar_vec = zeros(ncov,2);
           mix_ind = zeros(ncov,2);
           for i = 1:ncov
           mix = randperm(4);
           mix_ind = mix(1:2);
           covar_vec(i,:) = merge_covar_vec(i,mix_ind);
           end
        case 'low-high'
           merge_covar_vec = [covar_vec_low covar_vec_high];
           ncov = size(merge_covar_vec,1);
           covar_vec = zeros(ncov,2);
           mix_ind = zeros(ncov,2);
           for i = 1:ncov
           mix = randperm(4);
           mix_ind = mix(1:2);
           covar_vec(i,:) = merge_covar_vec(i,mix_ind);
           end
        case 'mid-high'
           merge_covar_vec = [covar_vec_mid covar_vec_high];
           ncov = size(merge_covar_vec,1);
           covar_vec = zeros(ncov,2);
           mix_ind = zeros(ncov,2);
           for i = 1:ncov
           mix = randperm(4);
           mix_ind = mix(1:2);
           covar_vec(i,:) = merge_covar_vec(i,mix_ind);
           end
        case 'low-mid-high'
           merge_covar_vec = [covar_vec_low covar_vec_mid covar_vec_high];
           ncov = size(merge_covar_vec,1);
           covar_vec = zeros(ncov,2);
           mix_ind = zeros(ncov,2);
           for i = 1:ncov
           mix = randperm(6);
           mix_ind = mix(1:2);
           covar_vec(i,:) = merge_covar_vec(i,mix_ind);
           end
        case 'mix'
           merge_covar_vec = [covar_vec_low covar_vec_mid covar_vec_high];
           mix_covar_vec = merge_covar_vec;
           ncov = size(merge_covar_vec,1);
           covar_vec = zeros(ncov,2);
           for i = 1:6
               mix_ind = randperm(ncov);
               mix_covar_vec(:,i) = merge_covar_vec(mix_ind,i);
           end          
           for i = 1:ncov
           mix_ind = randperm(6);
           mix_covar_vec(i,:) = merge_covar_vec(i,mix_ind);
           end           
           covar_vec = mix_covar_vec(:,1:2);
    end
covar_vec_2 = sortrows(covar_vec,1);
covar_vec_2 = sort(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec_2 = sortrows(covar_vec_2,2);
covar_vec = sortrows(covar_vec_2,1);
ncovar_vec = size(covar_vec,1);
covar_vec_bk = covar_vec;
cov_coord = covar_vec;
ncov_coord = size(cov_coord,1);
clear covar_vec_2
end

% We can select the covarions based on the identification provided by any
% of the methods available to detect coevolution. This is useful to test 
% how good the statistics of that method are with a set of simulated MSAs.
% Possible methods are:
% 'MI','NMI','MIP','ZPX','ZPX2','ZRES',ZRES2','ZMI','ZNMI','OMES','McBASC',
% 'ELSC','fodorSCA','ssemSCA','rsemSCA,'ramaSCA','logR','DCA','dbZPX2',
% 'nbZPX2','dgbZPX2','none'.
% ZRES is the original algorithm by Little and Chen, ZRES2 produces the 
% same result using the ZPX algorithm.
% The function 'get_nmsa_covar_vec' is also independently useful to get the
% coevolution matrix and the covarions vector from any msa in matlab
% numeric format using any of the listed methods. Usage is:
% [coevol_matrix,coevol_vector] = get_nmsa_covar_vec(nmsa,fcov,'method')

if METHOD
    disp('Method based assignment of covarions');    
    [~,covar_vec_coev_method] = ...
        get_nmsa_covar_vec(REF_nmsa,fcov,cov_method);
    covar_vec = unique(covar_vec_coev_method,'rows');
    if isempty(covar_vec)
    err('No method selected...stopping execution.')
    return
    else
    covar_vec_2 = sortrows(covar_vec,1);
    covar_vec_2 = sort(covar_vec_2,2);
    covar_vec_2 = sortrows(covar_vec_2,1);
    covar_vec_2 = sortrows(covar_vec_2,2);
    covar_vec = sortrows(covar_vec_2,1);
    ncovar_vec = size(covar_vec,1);
    covar_vec_bk = covar_vec;
    cov_coord = covar_vec;
    ncov_coord = size(cov_coord,1);
    clear covar_vec_2
    end

% Automatic activation of MULTIPLETS if any position appears in more
% than one pair.
%
%     u_covarions = unique(covar_vec);
%     for i = 1:numel(u_covarions)
%         u_covarions_ind = covarions == u_covarions(i);
%         mult_covarions = sum(u_covarions_ind(:));        
%         if mult_covarions > 1
%             MULTIPLETS = 1;
%             break
%         end
%     end
end

if VECTOR
    disp('User-list based assignment of covarions'); 
    covarions = unique(covarions,'rows');
    u_covarions = unique(covarions);

% Automatic activation of MULTIPLETS if any position appears in more
% than one pair.
%    
%     for i = 1:numel(u_covarions)
%         u_covarions_ind = covarions == u_covarions(i);
%         mult_covarions = sum(u_covarions_ind(:));        
%         if mult_covarions > 1
%             MULTIPLETS = 1;
%             break
%         end
%     end
    
    [ncov,ncovcols] = size(covarions);
    for i = 1:ncov
    temprow = nonzeros(covarions(i,:))';
    tempcols = size(temprow,2);
    temprow = [temprow zeros(1,ncovcols - tempcols)];
    covarions(i,:) = temprow;
    end
    clear temprow tempcols

    covar_vec = covarions;
    covar_vec_bk = covarions;
    cov_coord = covarions;
    ncov_coord = size(cov_coord,1);

    % If there are more than two columns in the covarion matrix, we must 
    % find the coordinates at which to evaluate covariation in the 
    % statistical tests.
    
    if size(covarions,2)>2

    % Automatic activation of MULTIPLETS if the 'covarions' matrix has more 
    % than 2 columns.
    %   MULTIPLETS = 1;

        u_covarions = unique(covarions);    
        m = 1;
        for i = 1:numel(u_covarions)
            for j = 1:numel(u_covarions)           
            m = m+1;
            pair(m,1) = u_covarions(i);
            pair(m,2) = u_covarions(j);
            end
        end
        pair_2 = sortrows(pair,1);
        pair_2 = sort(pair_2,2);
        pair_2 = sortrows(pair_2,1);
        pair_2 = sortrows(pair_2,2);
        pair = sortrows(pair_2,1);
        pair = unique(pair,'rows');
        same_ind = pair(:,1) == pair(:,2);
        pair(same_ind,:) = [];
        zero_ind = pair(:,1) == 0 | pair(:,2) == 0;
        pair(zero_ind,:) = [];
        cov_coord = pair;    
        ncov_coord = size(cov_coord,1);
    end    
end

% Reorder covar_vec with the first elements of the pairs in ascending order.
% IMPORTANT! This step is critical, because if the numbers are not
% ascending along the rows of the covar_vec matrix, reflecting the order of 
% the columns in the msa, the alphabet is inverted. This is the last point
% we can reorder, because if MULTIPLETS is turned on, no more reordering
% can be done before the PD's are calculated.
covar_vec_2 = sortrows(covar_vec,1);
covar_vec_2 = sort(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec_2 = sortrows(covar_vec_2,2);
covar_vec = sortrows(covar_vec_2,1);
clear covar_vec_2

% Make sure again covar_vec does not carry any leading zeros.

[ncov,ncovcols] = size(covar_vec);
for i = 1:ncov
    temprow = nonzeros(covar_vec(i,:))';
    tempcols = size(temprow,2);
    temprow = [temprow zeros(1,ncovcols - tempcols)];
    covar_vec(i,:) = temprow;
end
clear temprow tempcols
covar_vec = unique(covar_vec,'rows');
[ncov,ncolcov] = size(covar_vec);

%--------------------------------------------------------------------------
if MULTIPLETS
% Here we determine a larger covarion matrix where all the possible
% multiplets based on individual pairs are listed.
    
covar_vec_unique = unique(covar_vec);
temp_mat_lin = zeros(numel(covar_vec_unique),multiples);
for i = 1:numel(covar_vec_unique)
    pair_ind = covar_vec == covar_vec_unique(i);
    pair_ind = any(pair_ind,2);
    temp_mat = covar_vec(pair_ind,:);
    unique_covs = unique(temp_mat(:))';
    nunique = numel(unique_covs);
    temp_mat_lin(i,1:nunique) = unique_covs; 
    clear pair_ind
end
dim1 = any(temp_mat_lin,2);
dim2 = any(temp_mat_lin,1);
temp_mat_lin = temp_mat_lin(dim1,dim2);
temp_mat_lin = unique(temp_mat_lin,'rows');
[n_temp_mat_lin,ncols_temp_mat_lin] = size(temp_mat_lin);
temp_mat_lin = ...
    [temp_mat_lin zeros(n_temp_mat_lin,multiples-ncols_temp_mat_lin)];
clear dim1 dim2 ncols_temp_mat_lin

temp_mat_lin2 = zeros(numel(covar_vec_unique),multiples);
for i = 1:numel(covar_vec_unique)
    for j = 1:n_temp_mat_lin
    pair_ind(j) = sum((temp_mat_lin(j,:) == covar_vec_unique(i)));
    end
    pair_ind = logical(pair_ind);
    temp_mat = temp_mat_lin(pair_ind,:);
    unique_covs = unique(nonzeros(temp_mat(:)))';
    nunique = numel(unique_covs);
    temp_mat_lin2(i,1:nunique) = unique_covs; 
    clear pair_ind
end
temp_mat_lin2 = unique(temp_mat_lin2,'rows');
temp_mat_lin3 = [[covar_vec, zeros(ncov,multiples-ncolcov)];...
    temp_mat_lin; temp_mat_lin2];
temp_mat_lin3 = unique(temp_mat_lin3,'rows');
covar_vec = temp_mat_lin3;

[ncov,ncovcols] = size(covar_vec);
for i = 1:ncov
    temprow = nonzeros(covar_vec(i,:))';
    tempcols = size(temprow,2);
    temprow = [temprow zeros(1,ncovcols - tempcols)];
    covar_vec(i,:) = temprow;
end
dim1 = any(covar_vec,2);
dim2 = any(covar_vec,1);
covar_vec = covar_vec(dim1,dim2);

clear dim1 dim2 n_temp_mat_lin 
clear temprow tempcols

end   % End of MULTIPLETS
%--------------------------------------------------------------------------

[ncov,ncovcols] = size(covar_vec);

% Covarying PD of each pair or multiple in the entire MSA.
REF_cov_prob_PD = cell(ncov,1);
REF_cov_alphabet = zeros(nseq,ncovcols+1,ncov);

for i = 1:ncov
nmsa = REF_nmsa(:,nonzeros(covar_vec(i,:)));
if cov_gaps == 0;
    nmsa25 = nmsa == 25;
    nmsa25 = sum(nmsa25,2);
    nogaps = nmsa25 == 0;
    nmsa = nmsa(nogaps,:);
end
alphabet = JointProbDistr_9(nmsa);
    [n_alphabet,alphabet_cols] = size(alphabet);
    REF_cov_alphabet(1:n_alphabet,1:alphabet_cols,i) = alphabet;
    REF_cov_prob_PD{i} = ...
        fitdist((1:1:n_alphabet)','kernel',...
        'width',0.01,'frequency',...
        round(REF_cov_alphabet(1:n_alphabet,alphabet_cols,i)*1e6));    
    clear nmsa nmsa25 nogaps alphabet n_alphabet
end

% Covarying PD of each pair or multiple in the main branches.
REF_branch_cov_prob_PD = cell(ncov,nbranches);
REF_branch_cov_alphabet = zeros(nseq,ncovcols+1,ncov,nbranches);

for k = 1:nbranches
REF_branch_ind = find(REF_Clusters == k);
REF_nmsa_branch = REF_nmsa(REF_branch_ind,:);    

for i = 1:ncov
nmsa = REF_nmsa_branch(:,nonzeros(covar_vec(i,:)));
if cov_gaps == 0;
    nmsa25 = nmsa == 25;
    nmsa25 = sum(nmsa25,2);
    nogaps = nmsa25 == 0;
    nmsa = nmsa(nogaps,:);
end
alphabet = JointProbDistr_9(nmsa);
    [n_alphabet,alphabet_cols] = size(alphabet);
    REF_branch_cov_alphabet(1:n_alphabet,1:alphabet_cols,i,k) = alphabet;
    REF_branch_cov_prob_PD{i,k} = ...
        fitdist((1:1:n_alphabet)','kernel',...
        'width',0.01,'frequency',...
        round(REF_branch_cov_alphabet(1:n_alphabet,alphabet_cols,i,k)*1e6));    
    clear nmsa nmsa25 nogaps alphabet n_alphabet
end
clear REF_nmsa_branch
end

    if CHECK_BRANCHES
        [~] = make_REF_MSA_comp_1(REF_nmsa,DISTANCE,DIMENSIONS,...
        MAX_CLUSTERS,REF_rel_entropy,cov_coord,TREE_ONLY,...
        distance_method,tree_method,neighjoin_method);    
    % Pass control to the keyboard
    keyboard
    end
               
else  % Else for SET_COVARIONS statement.
    
% Here we set the maximum number of covarions to initialize all the
% arrays in the next section.
    
    RAND_COVARIONS = 1;
    ncov_coord = round(npos*fcov/100);
    gncov = ncov_coord;
    covar_vec = zeros(gncov,2);
    covar_vec_bk = covar_vec;
    cov_coord = covar_vec;
    
end   % End of SET_COVARIONS statement.

% The following line is critical: don't remove it.
covar_vec_large = covar_vec;

%% 'random msa' section
% Here we generate random MSA's simply using the probability distributions 
% of the different aa's at all non-covarying positions and the probability
% distributions of the covarions at all covarying positions. There is no
% evolution, and therefore no recombination, involved in this process.

if RANDOM_MSA
all_branch_size = sum(s_MSA_branch_size);
RAND_MSA_select = zeros(all_branch_size,npos,random_msa_no);
RAND_MSA_rel_entropy = zeros(random_msa_no,npos);
RAND_MSA_entropy = zeros(random_msa_no,npos);
corr_REF_RAND_MSA_rel_entropy = zeros(random_msa_no,1);
corr_REF_RAND_MSA_entropy = zeros(random_msa_no,1);
RAND_MSA_sim_score = zeros(random_msa_no,1);
RAND_start_time = tic;

    if WITH_COVARIONS
        for i = 1:random_msa_no
        branch_start = 0;
        for k = 1:nbranches
                mut_prob_PD = REF_branch_mut_prob_PD(:,k);
                mut_freq_PD = REF_branch_mut_freq_PD(:,:,k);
                cov_prob_PD = REF_branch_cov_prob_PD(:,k);
                cov_alphabet = REF_branch_cov_alphabet(:,:,:,k);

            branch_size = s_MSA_branch_size(k);
            branch_end = branch_start+branch_size;        
            branch_start = branch_start+1;
            RAND_MSA_select(:,:,i) = ...
            make_ancestor_3B(RAND_MSA_select(:,:,i),branch_start,...
            branch_end,npos,mut_prob_PD,mut_freq_PD,random_cycles);
%             RAND_MSA_select(:,:,i) = ...
%             make_random_covarions_1(RAND_MSA_select(:,:,i),branch_start,...
%             branch_end,covar_vec,cov_prob_PD,cov_alphabet);
            RAND_MSA_select(:,:,i) = ...        
            make_random_covarions_1B(RAND_MSA_select(:,:,i),branch_start,...
            branch_end,covar_vec,cov_prob_PD,cov_alphabet,random_cov_cycles);
        
            branch_start = branch_end;
        end
        end
    else
        for i = 1:random_msa_no
        branch_start = 0;
        for k = 1:nbranches
                mut_prob_PD = REF_branch_mut_prob_PD(:,k);
                mut_freq_PD = REF_branch_mut_freq_PD(:,:,k);

            branch_size = s_MSA_branch_size(k);
            branch_end = branch_start+branch_size;        
            branch_start = branch_start+1;
            RAND_MSA_select(:,:,i) = ...
            make_ancestor_3B(RAND_MSA_select(:,:,i),branch_start,...
            branch_end,npos,mut_prob_PD,mut_freq_PD,random_cycles);
            branch_start = branch_end;        
        end
        end    
    end

% Here we calculate relative entropies and mean similarity scores.
for i = 1:random_msa_no
    s_RAND_MSA_select = nmsa_to_smsa(RAND_MSA_select(:,:,i));
    [RAND_MSA_profile,~] = seqprofile(s_RAND_MSA_select,'Gaps','none');
    [~,RAND_MSA_rel_entropy(i,:)] = ...
        rel_entropy_nmsa(RAND_MSA_select(:,:,i),...
        bg_frequencies,RAND_MSA_profile);
    RAND_MSA_entropy(i,:) = Entropy(RAND_MSA_select(:,:,i));
    corr_REF_RAND_MSA_rel_entropy(i) = ...
        corr(REF_rel_entropy',RAND_MSA_rel_entropy(i,:)');
    corr_REF_RAND_MSA_entropy(i) = ...
        corr(REF_entropy',RAND_MSA_entropy(i,:)');    
end

if DISTANCE
    for i = 1:random_msa_no
        c_RAND_MSA_select = int2aa(RAND_MSA_select(:,:,i));
        RAND_MSA_dist = ...
            seqpdist(c_RAND_MSA_select,'SquareForm',...
            true,'Method','p-distance');
        RAND_MSA_cov_mat = 1-RAND_MSA_dist;
        RAND_MSA_sim_score(i,1) = multidim_MI(RAND_MSA_cov_mat);
    end
else
    for i = 1:random_msa_no    
        RAND_MSA_binmsa = nmsa_to_binmsa(RAND_MSA_select(:,:,i));
        RAND_MSA_cov_mat = cov(RAND_MSA_binmsa',1);    
        RAND_MSA_sim_score(i,1) = multidim_MI(RAND_MSA_cov_mat);
    end
end

% Check the execution time.
RAND_tElapsed = toc(RAND_start_time);
fprintf('RAND msa calculation: %d minutes \n', RAND_tElapsed/60);
clear RAND_start_time RAND_tElapsed

if random_msa_no > 1
% Here we look at the relationship between the consensus sequence in the
% simulated MSA and the ancestor.
    
RAND_MSA_sim_score_gevPD = fitdist(RAND_MSA_sim_score,'normal');
RAND_MSA_sim_score_gevPD_params(1,1) = mean(RAND_MSA_sim_score);
RAND_MSA_sim_score_gevPD_params(1,2) = std(RAND_MSA_sim_score);
RAND_MSA_sim_score_gevPD_params(1,3:4) = RAND_MSA_sim_score_gevPD.Params;

mean_RAND_MSA_entropy = mean(RAND_MSA_entropy);
mean_RAND_MSA_rel_entropy = mean(RAND_MSA_rel_entropy);
mean_RAND_MSA_cov_mat = mean(RAND_MSA_cov_mat_ALL,3);
end

end % End of RANDOM_MSA

%% 'evolution' section: declaration of variables

% Don't touch this section!

MSA_sim_score = zeros(end_cycle,1);
MSA_cons_anc_sim_score = zeros(end_cycle,1);
MSA_consensus = zeros(end_cycle,npos);
MSA_consensus_score = zeros(end_cycle,npos);
MSA_entropy = zeros(end_cycle,npos);
MSA_joint_entropy = zeros(npos,npos,end_cycle);
MSA_rel_entropy = zeros(end_cycle,npos);
emission_corr = zeros(end_cycle,npos);
mean_emission_corr = zeros(end_cycle,1);
COV_ALL = zeros(npos,npos,end_cycle);
COV_bk_ALL = zeros(npos,npos,end_cycle);
COV_ntz_ALL = zeros(npos,npos,end_cycle);
recomb_COV_ALL = zeros(npos,npos,end_cycle);
recomb_COV_bk_ALL = zeros(npos,npos,end_cycle);
recomb_COV_ntz_ALL = zeros(npos,npos,end_cycle);
mut_COV_ALL = zeros(npos,npos,end_cycle);
mut_COV_bk_ALL = zeros(npos,npos,end_cycle);
mut_COV_ntz_ALL = zeros(npos,npos,end_cycle);
cov_COV_ALL = zeros(npos,npos,end_cycle);
cov_COV_bk_ALL = zeros(npos,npos,end_cycle);
cov_COV_ntz_ALL = zeros(npos,npos,end_cycle);
glob_COV_ALL = zeros(npos,npos,end_cycle);
glob_COV_bk_ALL = zeros(npos,npos,end_cycle);
glob_COV_ntz_ALL = zeros(npos,npos,end_cycle);
mutcov_COV_ALL = zeros(npos,npos,end_cycle);
mutcov_COV_bk_ALL = zeros(npos,npos,end_cycle);
mutcov_COV_ntz_ALL = zeros(npos,npos,end_cycle);

fcov_COV = zeros(end_cycle,3);
fcov_recomb_COV = zeros(end_cycle,3);
fcov_cov_COV = zeros(end_cycle,3);
fcov_mut_COV = zeros(end_cycle,3);
fcov_glob_COV = zeros(end_cycle,3);

cov_zscore_COV = zeros(end_cycle,ncov_coord);
covar_vec_zscore_COV = zeros(ncov_coord,3,end_cycle);
cov_zscore_recomb_COV = zeros(end_cycle,ncov_coord);
covar_vec_zscore_recomb_COV = zeros(ncov_coord,3,end_cycle);
cov_zscore_mut_COV = zeros(end_cycle,ncov_coord);
covar_vec_zscore_mut_COV = zeros(ncov_coord,3,end_cycle);
cov_zscore_cov_COV = zeros(end_cycle,ncov_coord);
covar_vec_zscore_cov_COV = zeros(ncov_coord,3,end_cycle);
cov_zscore_glob_COV = zeros(end_cycle,ncov_coord);
covar_vec_zscore_glob_COV = zeros(ncov_coord,3,end_cycle);
 
corr_REF_MSA_rel_entropy = zeros(end_cycle,1);
corr_REF_MSA_s_consensus_score = zeros(end_cycle,1);
neigmod_ALL = zeros(end_cycle,1);

MSA1_nseq = sum(bncopies_1);

MSA1 = zeros(sum(bncopies_1),npos);
MSA2 = zeros(nselect,npos);
MSA1_history = zeros(sum(bncopies_1),npos);
MSA2_history = zeros(nselect,npos);

MSA1_recomb_history = zeros(sum(bncopies_1),npos);
MSA2_recomb_history = zeros(nselect,npos);

MSA1_mut_history = zeros(sum(bncopies_1),npos);
MSA2_mut_history = zeros(nselect,npos);

MSA1_cov_history = zeros(sum(bncopies_1),npos);
MSA2_cov_history = zeros(nselect,npos);

MSA2_ALL = zeros(nselect,npos,end_cycle);
MSA2_history_ALL = zeros(nselect,npos,end_cycle);
MSA2_recomb_history_ALL = zeros(nselect,npos,end_cycle);
MSA2_mut_history_ALL = zeros(nselect,npos,end_cycle);
MSA2_cov_history_ALL = zeros(nselect,npos,end_cycle);

MSA_select_ALL = zeros(nselect,npos,end_cycle);
MSA_select_orig_ALL = zeros(nselect,npos,end_cycle);
MSA_cov_mat_ALL = zeros(nselect,nselect,end_cycle);
ancestor_ALL = zeros(end_cycle,npos);

select_orig_ind = false(nselect,end_cycle);

branch_emission_corr = zeros(nbranches,npos,end_cycle);
mean_branch_emission_corr = zeros(nbranches,end_cycle);
branch_profile_corr = zeros(nbranches,npos,end_cycle);
mean_branch_profile_corr = zeros(nbranches,end_cycle);

if SAVE_ALL_COV  % CAREFUL! Huge memory requirements if more than 1 cycle.
COV_3D_ALL= zeros(npos,npos,nselect,end_cycle);
cov_COV_3D_ALL= zeros(npos,npos,nselect,end_cycle);
mut_COV_3D_ALL= zeros(npos,npos,nselect,end_cycle);
recomb_COV_3D_ALL= zeros(npos,npos,nselect,end_cycle);
end

% End of the declaration of variables cell. If you are running in cell mode
% and you skipped the 'RUN' cell, now you should jump to the 'evolution'
% cell.
%
end % end of NEW_RUN.

%% 'evolution' section: tree build

% Don't touch anything in the following section unless you really know what 
% you are doing!

WRITE_RNG = true;

% Change the start_cycle from 1 to the last saved cycle to restart an 
% interrupted run of multiple cycles, and run just this cell.

% start_cycle = 1;              % Default = assigned by RUNPAR
% end_cycle = no. of cycles;    % Default = assigned by RUNPAR

if start_cycle > 1

% Clear some variables possibly left over from an earlier cycle.    

clear *orig MSA0 MSA1 MSA2
clear REF_tree ncopies_11 ncopies_12 ncopies
clear ncopies_21 ncopies_22 ncopies_31 ncopies_32
clear COV_sum glob_COV_sum mut_COV_sum
clear cov_COV_sum mutcov_COV_sum recomb_COV_sum
clear sumCOVdif*
clear i j h k m A B D V gap name
clear vec1 vec2 a b ind
clear jc_MSA jc_MSA0 data
clear algn algn_rnd neigmod
clear ev lbd ev_rnd lbd_rnd ev_rnd_tmp lbd_rnd_tmp lbd_mat
clear N_ev N_samples N_pos N_seq C_rnd C_tmp perm_seq
clear calc_start_time cov_calc_time

% Restore the correct paths.

% rmpath('SCA_Sharma');
% rmpath(genpath('SCA_4.5'));
addpath('DWINNEL_MI');
addpath('MATLAB_LOCAL');

% Restore the last random number sequence used.

rng(srng(start_cycle));
WRITE_RNG = false;

end

for n = start_cycle:end_cycle

% CRITICAL!!! We restore here the value of covar_vec that is used in the
% calculations of the statistical properties of all covariation matrices to 
% covar_vec_large to match all the arrays.

covar_vec = covar_vec_large;

% Keep track of the cycle number:
fprintf('Evolution cycle %d \n', n);

% We save the random generator state in case we want to reproduce a
% specific run. The WRITE_RNG variable allows a bitwise restart from an 
% interrupted run.

if WRITE_RNG
 srng(n) = rng;
else
 WRITE_RNG = true;    
end

% For example if we wanted to reproduce exactly the 3rd run we would unflag
% the following line.
% rng(srng(3));

% An ancestor is created based on the background probabilities or the 
% emission probabilities of the hmm model.
% The ancestor has 'nbranches' copies of the same sequence.

if ANCESTOR_BG_PROB    
ancestor = make_ancestor(nbranches,npos,bg_prob_PD);
else
ancestor = make_ancestor_2(nbranches,npos,mut_prob_PD);
end

ancestor_ALL(n,:) = ancestor(1,:);

% Covarions can be assigned randomly at each cycle. If we want 10% of
% all the positions (eg 300) to be covarying with another 10%, 
% we will randomly pick two vectors of 30 columns (ncov = 30). In practice, 
% since there can be some overlap in the random choice of the positions,
%  we may end up with less than 60 unique columns. For this reason we
% multiply ncov by ncov_mult. A value of 3 is more than safe, but if  
% there is an error because some covar_vec's have less than ncov elements
% increase 'ncov_mult' to 4 or 5. However 'ncov_mult*ncov' cannot be larger 
% than the total number of positions npos.

if RAND_COVARIONS
    ncov_mult = 3.0;
    covar_vec_rand = ...
        pick_covarions_2(REF_nmsa,ncov_coord,ncov_mult,REF_conserved_ind);    
    run PROCESS_COVARIONS_2    
    covar_vec_ALL(:,:,n) = covar_vec;
end

% Recombination zones can be assigned randomly at each evolution cycle.

if RAND_RECOMB_ZONES
    nedges = nzones - 1;
    hrf = randi(npos,1,nedges);
    hrf = sort(hrf,2);
    hrf = [1 hrf npos];
    hrf_ALL(n,:) = hrf;
end

% Here we start the clock.

evol_start_time = tic;
    
MSA0 = ancestor;
tnseq = size(MSA0,1);

% Here we create a numeric array with all zeros of the same size as
% ancestor for the various types of hystory. It is critical to make sure
% that the history arrays are always numeric in the evolution and
% recombination functions below to avoid unpredictable results.

MSA0_history = zeros(tnseq,npos);
MSA0_recomb_history = zeros(tnseq,npos);
MSA0_mut_history = zeros(tnseq,npos);
MSA0_cov_history = zeros(tnseq,npos);
COV = zeros(npos,npos,tnseq);
recomb_COV = zeros(npos,npos,tnseq);
mut_COV = zeros(npos,npos,tnseq);
cov_COV = zeros(npos,npos,tnseq);
glob_COV = zeros(npos,npos,tnseq);
mutcov_COV = zeros(npos,npos,tnseq);

% 1st node of evolution. 
% At the beginning the tree includes only the ancestor. We also set up the 
% coevolution and history matrices that we will keep updating.

%%%% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
%%%%

% mut_rate = randi(mut_rate_range_1);

for i = 1:nbranches
    mut_rate = randi(mut_rate_range_1);
    mut_prob_PD = REF_branch_mut_prob_PD(:,i);
    mut_freq_PD = REF_branch_mut_freq_PD(:,:,i);
    cov_prob_PD = REF_branch_cov_prob_PD(:,i);
    cov_alphabet = REF_branch_cov_alphabet(:,:,:,i);

    for m = 1:mut_cycles    
    [MSA0(i,:),MSA0_history(i,:),MSA0_mut_history(i,:),...
        MSA0_cov_history(i,:),COV(:,:,i),mut_COV(:,:,i),...
        cov_COV(:,:,i),glob_COV(:,:,i),mutcov_COV(:,:,i)] = ...
        evolve_6s_3D_cov_hist_multiplets(MSA0(i,:),MSA0_history(i,:),...
        MSA0_mut_history(i,:),MSA0_cov_history(i,:),...
        COV(:,:,i),mut_COV(:,:,i),cov_COV(:,:,i),...
        glob_COV(:,:,i),mutcov_COV(:,:,i),...
        mut_rate,covar_vec,mut_prob_PD,...
        cov_prob_PD,cov_alphabet,mut_freq_PD);
    end        
end

%%%% For testing:
% for i = 1:nbranches
% dif_COV_HIST = sum(diag(COV(:,:,i))) - sum(MSA0_history(i,:))
% dif_globCOV_HIST = sum(diag(glob_COV(:,:,i))) - sum(MSA0_history(i,:))
% dif_mutCOV_HIST = sum(diag(mut_COV(:,:,i))) - sum(MSA0_mut_history(i,:))
% dif_covCOV_HIST = sum(diag(cov_COV(:,:,i))) - sum(MSA0_cov_history(i,:))
% end

% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

% Recombination between all the members of the original subgroup (MSA0);

%%%% For testing:
% COV_orig = COV;
% glob_COV_orig = glob_COV;
% recomb_COV_orig = recomb_COV;
%%%%

for r = 1:rec_cycles
nhrt = randi(nhrt_range_1);
 [MSA0,MSA0_history,MSA0_recomb_history,COV,recomb_COV,glob_COV] = ...
    hor_transfer_6s_3D_cov_hist(MSA0,MSA0_history,MSA0_recomb_history,...
    COV,recomb_COV,glob_COV,hrf,nzones,nhrt);
end

%%%% For testing:
% for i = 1:nbranches
% dif_COV_HIST = sum(diag(COV(:,:,i))) - sum(MSA0_history(i,:))
% dif_globCOV_HIST = sum(diag(glob_COV(:,:,i))) - sum(MSA0_history(i,:))
% dif_recombCOV_HIST = ...
%     sum(diag(recomb_COV(:,:,i))) - sum(MSA0_recomb_history(i,:))
% end

% sumCOVdif = testrec(COV,COV_orig,recomb_COV,recomb_COV_orig)
% sumCOVdif2 = testrec_2(glob_COV,glob_COV_orig,recomb_COV,recomb_COV_orig)
%%%%

% Evolution after the recombination.

%%%% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
%%%%

for i = 1:nbranches
    mut_rate = randi(mut_rate_range_2);
    mut_prob_PD = REF_branch_mut_prob_PD(:,i);
    mut_freq_PD = REF_branch_mut_freq_PD(:,:,i);
    cov_prob_PD = REF_branch_cov_prob_PD(:,i);
    cov_alphabet = REF_branch_cov_alphabet(:,:,:,i);

    for m = 1:mut_cycles    
    [MSA0(i,:),MSA0_history(i,:),MSA0_mut_history(i,:),...
        MSA0_cov_history(i,:),COV(:,:,i),mut_COV(:,:,i),...
        cov_COV(:,:,i),glob_COV(:,:,i),mutcov_COV(:,:,i)] = ...
        evolve_6s_3D_cov_hist_multiplets(MSA0(i,:),MSA0_history(i,:),...
        MSA0_mut_history(i,:),MSA0_cov_history(i,:),...
        COV(:,:,i),mut_COV(:,:,i),cov_COV(:,:,i),...
        glob_COV(:,:,i),mutcov_COV(:,:,i),...
        mut_rate,covar_vec,mut_prob_PD,...
        cov_prob_PD,cov_alphabet,mut_freq_PD);
    end    
end

%%%% For testing:
% for i = 1:nbranches
% dif_COV_HIST = sum(diag(COV(:,:,i))) - sum(MSA0_history(i,:))
% dif_globCOV_HIST = sum(diag(glob_COV(:,:,i))) - sum(MSA0_history(i,:))
% dif_mutCOV_HIST = sum(diag(mut_COV(:,:,i))) - sum(MSA0_mut_history(i,:))
% dif_covCOV_HIST = sum(diag(cov_COV(:,:,i))) - sum(MSA0_cov_history(i,:))
% end

% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

% Each of the 'nbranches' different sequences derived from evolving 
% the ancestor gets copied a number of times equal to "ncopies_1". 
% COV is the matrix that holds the true coevolution between 
% positions. Accordingly, we also multiply by "ncopies_1" the COV matrix.

MSA1 = zeros(MSA1_nseq,npos);
ind = zeros(nbranches,2);
for i = 1:nbranches
    ncopies_1 = bncopies_1(i);
    ind(i,:) = [(ind(i,1)+1) (ind(i,2)+ncopies_1)];
    ind(i+1,:) = [ind(i,2) ind(i,2)]; 
end    
b_ind = ind(1:nbranches,:);
clear ind
for i = 1:nbranches
 ncopies_1 = bncopies_1(i);
 ind = b_ind(i,1):b_ind(i,2);
 MSA1(ind,:) = repmat(MSA0(i,:),ncopies_1,1);
 MSA1_history(ind,:) = repmat(MSA0_history(i,:),ncopies_1,1);
 MSA1_recomb_history(ind,:) = repmat(MSA0_recomb_history(i,:),ncopies_1,1);
 MSA1_mut_history(ind,:) = repmat(MSA0_mut_history(i,:),ncopies_1,1);
 MSA1_cov_history(ind,:) = repmat(MSA0_cov_history(i,:),ncopies_1,1);
 clear ind
end
    
% Here we proliferate the COV, mut_COV, cov_COV, and recomb_COV matrix by the 
% number of ncopies.

COV_orig = COV;
mut_COV_orig = mut_COV;
cov_COV_orig = cov_COV;
glob_COV_orig = glob_COV;
mutcov_COV_orig = mutcov_COV;
recomb_COV_orig = recomb_COV;

for i = 1:nbranches
ncopies_1 = bncopies_1(i);
ind = b_ind(i,1):b_ind(i,2);
COV(:,:,ind) = repmat(COV_orig(:,:,i),[1 1 ncopies_1]);
mut_COV(:,:,ind) = repmat(mut_COV_orig(:,:,i),[1 1 ncopies_1]);
cov_COV(:,:,ind) = repmat(cov_COV_orig(:,:,i),[1 1 ncopies_1]);
glob_COV(:,:,ind) = repmat(glob_COV_orig(:,:,i),[1 1 ncopies_1]);
mutcov_COV(:,:,ind) = repmat(mutcov_COV_orig(:,:,i),[1 1 ncopies_1]);
recomb_COV(:,:,ind) = repmat(recomb_COV_orig(:,:,i),[1 1 ncopies_1]);
clear ind
end

% Free up memory:
clear COV_orig mut_COV_orig cov_COV_orig glob_COV_orig 
clear mutcov_COV_orig recomb_COV_orig 

% 2nd node of evolution.

%%%% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
% recomb_COV_orig = recomb_COV;
%%%%

for i = 1:nbranches
ind = b_ind(i,1):b_ind(i,2);
mut_rate = randi(mut_rate_range_3);
mut_prob_PD = REF_branch_mut_prob_PD(:,i);
mut_freq_PD = REF_branch_mut_freq_PD(:,:,i);
cov_prob_PD = REF_branch_cov_prob_PD(:,i);
cov_alphabet = REF_branch_cov_alphabet(:,:,:,i);

for m = 1:mut_cycles
[MSA1(ind,:),MSA1_history(ind,:),MSA1_mut_history(ind,:),...
    MSA1_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
    cov_COV(:,:,ind),glob_COV(:,:,ind),mutcov_COV(:,:,ind)] = ...
    evolve_6s_3D_cov_hist_multiplets(MSA1(ind,:),...
    MSA1_history(ind,:),MSA1_mut_history(ind,:),MSA1_cov_history(ind,:),...
    COV(:,:,ind),mut_COV(:,:,ind),cov_COV(:,:,ind),glob_COV(:,:,ind),...
    mutcov_COV(:,:,ind),mut_rate,covar_vec,mut_prob_PD,...
    cov_prob_PD,cov_alphabet,mut_freq_PD);
end    
end

%%%% For testing:
% for i = 1:nbranches
% ind = b_ind(i,1):b_ind(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_mutCOV_HIST = ...
%     sum(diag(sum(mut_COV(:,:,ind),3))) - sum(sum(MSA1_mut_history(ind,:)))
% dif_covCOV_HIST = ...
%     sum(diag(sum(cov_COV(:,:,ind),3))) - sum(sum(MSA1_cov_history(ind,:)))
% end

% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

% n2nd repeats 'n2nd_cycle' times only the recombination between
% the main branches and the following evolution. This step is
% conceptually similar to 'nfinal_cycle'.

for n2nd = 1:n2nd_cycle

% Recombination between the branches. First we sum all the
% branches into a single array.

%%%% For testing:
% COV_orig = COV;
% glob_COV_orig = glob_COV;
% recomb_COV_orig = recomb_COV;
%%%%

% Here we recombine.

for r = 1:rec_cycles
nhrt = randi(nhrt_range_2);
[MSA1,MSA1_history,MSA1_recomb_history,COV,...
    recomb_COV,glob_COV] = ...
    hor_transfer_6s_3D_cov_hist(MSA1,MSA1_history,...
    MSA1_recomb_history,COV,recomb_COV,glob_COV,...
    hrf,nzones,nhrt);
end

%%%% For testing:
% for i = 1:nbranches
% ind = b_ind(i,1):b_ind(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_recombCOV_HIST = ...
%     sum(diag(sum(recomb_COV(:,:,ind),3))) - ...
%     sum(sum(MSA1_recomb_history(ind,:)))
% end

% sumCOVdif = testrec(COV,COV_orig,recomb_COV,recomb_COV_orig)
% sumCOVdif2 = testrec_2(glob_COV,glob_COV_orig,recomb_COV,recomb_COV_orig)
%%%%

% Evolution after the recombination.

% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
% recomb_COV_orig = recomb_COV;
%%%%

for i = 1:nbranches
ind = b_ind(i,1):b_ind(i,2);
mut_rate = randi(mut_rate_range_4);
mut_prob_PD = REF_branch_mut_prob_PD(:,i);
mut_freq_PD = REF_branch_mut_freq_PD(:,:,i);
cov_prob_PD = REF_branch_cov_prob_PD(:,i);
cov_alphabet = REF_branch_cov_alphabet(:,:,:,i);

for m = 1:mut_cycles
[MSA1(ind,:),MSA1_history(ind,:),MSA1_mut_history(ind,:),...
    MSA1_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
    cov_COV(:,:,ind),glob_COV(:,:,ind),mutcov_COV(:,:,ind)] = ...
    evolve_6s_3D_cov_hist_multiplets(MSA1(ind,:),...
    MSA1_history(ind,:),MSA1_mut_history(ind,:),MSA1_cov_history(ind,:),...
    COV(:,:,ind),mut_COV(:,:,ind),cov_COV(:,:,ind),glob_COV(:,:,ind),...
    mutcov_COV(:,:,ind),mut_rate,covar_vec,mut_prob_PD,...
    cov_prob_PD,cov_alphabet,mut_freq_PD);
end    
end

%%%% For testing:
% for i = 1:nbranches
% ind = b_ind(i,1):b_ind(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_mutCOV_HIST = ...
%     sum(diag(sum(mut_COV(:,:,ind),3))) - sum(sum(MSA1_mut_history(ind,:)))
% dif_covCOV_HIST = ...
%     sum(diag(sum(cov_COV(:,:,ind),3))) - sum(sum(MSA1_cov_history(ind,:)))
% end

% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

end % end of n2nd_cycle: we proceed with the rest of the 2nd node of 
    % evolution.
    
% Recombination inside each branch.

%%%% For testing:
% COV_orig = COV;
% glob_COV_orig = glob_COV;
% recomb_COV_orig = recomb_COV;
%%%%

for i = 1:nbranches
ind = b_ind(i,1):b_ind(i,2);
    
for r = 1:rec_cycles
nhrt = randi(round(nhrt_range_3*REF_branch_sim_score(i)));
[MSA1(ind,:),MSA1_history(ind,:),MSA1_recomb_history(ind,:),...
    COV(:,:,ind),recomb_COV(:,:,ind),glob_COV(:,:,ind)] = ...
    hor_transfer_6s_3D_cov_hist(MSA1(ind,:),MSA1_history(ind,:),...
    MSA1_recomb_history(ind,:),COV(:,:,ind),recomb_COV(:,:,ind),...
    glob_COV(:,:,ind),hrf,nzones,nhrt);
end
end

%%%% For testing:
% for i = 1:nbranches
% ind = b_ind(i,1):b_ind(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_recombCOV_HIST = ...
%     sum(diag(sum(recomb_COV(:,:,ind),3))) - ...
%     sum(sum(MSA1_recomb_history(ind,:)))
% end

% sumCOVdif = testrec(COV,COV_orig,recomb_COV,recomb_COV_orig)
% sumCOVdif2 = testrec_2(glob_COV,glob_COV_orig,recomb_COV,recomb_COV_orig)
%%%%

% Evolution after the recombination.

%%%% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
%%%%

for i = 1:nbranches
ind = b_ind(i,1):b_ind(i,2);
mut_rate = randi(mut_rate_range_5);
mut_prob_PD = REF_branch_mut_prob_PD(:,i);
mut_freq_PD = REF_branch_mut_freq_PD(:,:,i);
cov_prob_PD = REF_branch_cov_prob_PD(:,i);
cov_alphabet = REF_branch_cov_alphabet(:,:,:,i);

for m = 1:mut_cycles
[MSA1(ind,:),MSA1_history(ind,:),MSA1_mut_history(ind,:),...
    MSA1_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
    cov_COV(:,:,ind),glob_COV(:,:,ind),mutcov_COV(:,:,ind)] = ...
    evolve_6s_3D_cov_hist_multiplets(MSA1(ind,:),...
    MSA1_history(ind,:),MSA1_mut_history(ind,:),MSA1_cov_history(ind,:),...
    COV(:,:,ind),mut_COV(:,:,ind),cov_COV(:,:,ind),glob_COV(:,:,ind),...
    mutcov_COV(:,:,ind),mut_rate,covar_vec,mut_prob_PD,...
    cov_prob_PD,cov_alphabet,mut_freq_PD);
end    
end

%%%% For testing:
% for i = 1:nbranches
% ind = b_ind(i,1):b_ind(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA1_history(ind,:)))
% dif_mutCOV_HIST = ...
%     sum(diag(sum(mut_COV(:,:,ind),3))) - sum(sum(MSA1_mut_history(ind,:)))
% dif_covCOV_HIST = ...
%     sum(diag(sum(cov_COV(:,:,ind),3))) - sum(sum(MSA1_cov_history(ind,:)))
% end

% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

% Here each of the ncopies_1 different sequences inside each of the 
% original 'nbranches' branches gets copied a number of times equal to 
% 'ncopies_2'. We will also proliferate the coevolution matrices.
 
ind2 = zeros(nbranches,2);
for i = 1:nbranches
    ncopies_1_2 = bncopies_1(i)*bncopies_2(i);
    ind2(i,:) = [(ind2(i,1)+1) (ind2(i,2)+ncopies_1_2)];
    ind2(i+1,:) = [ind2(i,2) ind2(i,2)]; 
end    
b_ind2 = ind2(1:nbranches,:);
clear ind2

nsubbranches = sum(bncopies_2);
ind3 = zeros(nsubbranches,2);
b_ind3 = [];
for i = 1:nbranches
    for j = 1:bncopies_2(i)
    ind3(j,:) = [(ind3(j,1)+1) (ind3(j,2)+bncopies_1(i))];
    ind3(j+1,:) = [ind3(j,2) ind3(j,2)]; 
    end
%end
b_ind3_part = ind3(1:bncopies_2(i),:);
b_ind3 = [b_ind3 ; b_ind3_part];
ind3 = max(ind3);
end
clear ind3

for k = 1:nbranches
 ncopies_2 = bncopies_2(k);
 ind = b_ind(k,1):b_ind(k,2);
 ind2 = b_ind2(k,1):b_ind2(k,2);
 MSA2(ind2,:) = repmat(MSA1(ind,:),ncopies_2,1);    
 MSA2_history(ind2,:) = repmat(MSA1_history(ind,:),ncopies_2,1);        
 MSA2_recomb_history(ind2,:) = ...
     repmat(MSA1_recomb_history(ind,:),ncopies_2,1);
 MSA2_mut_history(ind2,:) = ...
     repmat(MSA1_mut_history(ind,:),ncopies_2,1);        
 MSA2_cov_history(ind2,:) = ...
         repmat(MSA1_cov_history(ind,:),ncopies_2,1);        
end

% Here we proliferate the COV, mut_COV, cov_COV, and recomb_COV matrix by the 
% number of ncopies_2

COV_orig = COV;
mut_COV_orig = mut_COV;
cov_COV_orig = cov_COV;
glob_COV_orig = glob_COV;
mutcov_COV_orig = mutcov_COV;
recomb_COV_orig = recomb_COV;

for k = 1:nbranches
 ncopies_2 = bncopies_2(k);
 ind = b_ind(k,1):b_ind(k,2);
 ind2 = b_ind2(k,1):b_ind2(k,2);
 COV(:,:,ind2) = repmat(COV_orig(:,:,ind),[1 1 ncopies_2]);
 mut_COV(:,:,ind2) = repmat(mut_COV_orig(:,:,ind),[1 1 ncopies_2]);
 cov_COV(:,:,ind2) = repmat(cov_COV_orig(:,:,ind),[1 1 ncopies_2]);
 glob_COV(:,:,ind2) = repmat(glob_COV_orig(:,:,ind),[1 1 ncopies_2]);
 mutcov_COV(:,:,ind2) = repmat(mutcov_COV_orig(:,:,ind),[1 1 ncopies_2]);
 recomb_COV(:,:,ind2) = repmat(recomb_COV_orig(:,:,ind),[1 1 ncopies_2]);
end

% Free up some memory:
clear COV_orig mut_COV_orig cov_COV_orig 
clear glob_COV_orig mutcov_COV_orig recomb_COV_orig 

% 3rd node of evolution. It can be repeated 'n3rd_cycle' times.

for n3rd = 1:n3rd_cycle

for k = 1:nbranches
    
mut_prob_PD = REF_branch_mut_prob_PD(:,k);
mut_freq_PD = REF_branch_mut_freq_PD(:,:,k);
cov_prob_PD = REF_branch_cov_prob_PD(:,k);
cov_alphabet = REF_branch_cov_alphabet(:,:,:,k);

%%%% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
%%%%

% Here we cycle through each sub-branch.

a = find(b_ind3(:,1) == b_ind2(k,1));
b = find(b_ind3(:,2) == b_ind2(k,2));

    for i = a:b
        ind = b_ind3(i,1):b_ind3(i,2);

        mut_rate = randi(mut_rate_range_6);    
        
        for m = 1:mut_cycles
        mut_rate_1 = mut_rate + randi([nmrm,pmrm]);    
            
        [MSA2(ind,:),MSA2_history(ind,:),MSA2_mut_history(ind,:),...
            MSA2_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
            cov_COV(:,:,ind),glob_COV(:,:,ind),...
            mutcov_COV(:,:,ind)] = ...
            evolve_6s_3D_cov_hist_multiplets(MSA2(ind,:),...
            MSA2_history(ind,:),MSA2_mut_history(ind,:),...
            MSA2_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
            cov_COV(:,:,ind),glob_COV(:,:,ind),...
            mutcov_COV(:,:,ind),...
            mut_rate_1,covar_vec,mut_prob_PD,...
            cov_prob_PD,cov_alphabet,mut_freq_PD);
        end
    end

%%%% For testing:
% for i = a:b
%       ind = b_ind3(i,1):b_ind3(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
% dif_mutCOV_HIST = ...
%     sum(diag(sum(mut_COV(:,:,ind),3))) - sum(sum(MSA2_mut_history(ind,:)))
% dif_covCOV_HIST = ...
%     sum(diag(sum(cov_COV(:,:,ind),3))) - sum(sum(MSA2_cov_history(ind,:)))
% end
    
% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

% Recombination between the sub-branches. 

%%%% For testing:
% COV_orig = COV;
% glob_COV_orig = glob_COV;
% recomb_COV_orig = recomb_COV;
%%%%
    
% Here we recombine.

ind = b_ind2(k,1):b_ind2(k,2);

    for r = 1:rec_cycles
    nhrt = randi(round(nhrt_range_4*REF_branch_sim_score(k)));
    [MSA2(ind,:),MSA2_history(ind,:),MSA2_recomb_history(ind,:),...
    COV(:,:,ind),recomb_COV(:,:,ind),glob_COV(:,:,ind)] = ...
    hor_transfer_6s_3D_cov_hist(MSA2(ind,:),MSA2_history(ind,:),...
    MSA2_recomb_history(ind,:),COV(:,:,ind),recomb_COV(:,:,ind),...
    glob_COV(:,:,ind),hrf,nzones,nhrt);
    end
    

%%%% For testing:
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
% dif_recombCOV_HIST = ...
%     sum(diag(sum(recomb_COV(:,:,ind),3))) - ...
%     sum(sum(MSA2_recomb_history(ind,:)))

% sumCOVdif = testrec(COV,COV_orig,recomb_COV,recomb_COV_orig)
% sumCOVdif2 = testrec_2(glob_COV,glob_COV_orig,recomb_COV,recomb_COV_orig)
%%%%

% Evolution after the recombination.

%%%% For testing:
% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
%%%%

% Here we cycle through each sub-branch.

a = find(b_ind3(:,1) == b_ind2(k,1));
b = find(b_ind3(:,2) == b_ind2(k,2));

    for i = a:b
        ind = b_ind3(i,1):b_ind3(i,2);

        mut_rate = randi(mut_rate_range_7);    
        
        for m = 1:mut_cycles
        mut_rate_1 = mut_rate + randi([nmrm,pmrm]);    
            
        [MSA2(ind,:),MSA2_history(ind,:),MSA2_mut_history(ind,:),...
            MSA2_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
            cov_COV(:,:,ind),glob_COV(:,:,ind),...
            mutcov_COV(:,:,ind)] = ...
            evolve_6s_3D_cov_hist_multiplets(MSA2(ind,:),...
            MSA2_history(ind,:),MSA2_mut_history(ind,:),...
            MSA2_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
            cov_COV(:,:,ind),glob_COV(:,:,ind),...
            mutcov_COV(:,:,ind),...
            mut_rate_1,covar_vec,mut_prob_PD,...
            cov_prob_PD,cov_alphabet,mut_freq_PD);
        end
    end

%%%% For testing:
% for i = a:b
%       ind = b_ind3(i,1):b_ind3(i,2);
% dif_COV_HIST = ...
%     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
% dif_globCOV_HIST = ...
%     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
% dif_mutCOV_HIST = ...
%     sum(diag(sum(mut_COV(:,:,ind),3))) - sum(sum(MSA2_mut_history(ind,:)))
% dif_covCOV_HIST = ...
%     sum(diag(sum(cov_COV(:,:,ind),3))) - sum(sum(MSA2_cov_history(ind,:)))
% end
    
% sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
%     cov_COV,cov_COV_orig)
% sumCOVdif2 = testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
%%%%

% The following two steps can be repeated 'nfinal_cycle' times to fine
% tune the final correlation between the reference and the simulated msa.

for nfinal = 1:nfinal_cycle

% Recombination inside each sub-branch.

    %%%% For testing:
    % COV_orig = COV;
    % glob_COV_orig = glob_COV;
    % recomb_COV_orig = recomb_COV;
    %%%%

% Here we cycle through each sub-branch.

a = find(b_ind3(:,1) == b_ind2(k,1));
b = find(b_ind3(:,2) == b_ind2(k,2));

    for i = a:b
        ind = b_ind3(i,1):b_ind3(i,2);
    
    for r = 1:rec_cycles
    nhrt = randi(round(nhrt_range_5*REF_branch_sim_score(k)));
    [MSA2(ind,:),MSA2_history(ind,:),MSA2_recomb_history(ind,:),...
    COV(:,:,ind),recomb_COV(:,:,ind),glob_COV(:,:,ind)] = ...
    hor_transfer_6s_3D_cov_hist(MSA2(ind,:),MSA2_history(ind,:),...
    MSA2_recomb_history(ind,:),COV(:,:,ind),recomb_COV(:,:,ind),...
    glob_COV(:,:,ind),hrf,nzones,nhrt);
    end
    end

    
    %%%% For testing:
    % for i = a:b
    %     ind = b_ind3(i,1):b_ind3(i,2);
    % dif_COV_HIST = ...
    %     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
    % dif_globCOV_HIST = ...
    %     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
    % dif_recombCOV_HIST = ...
    %     sum(diag(sum(recomb_COV(:,:,ind),3))) - ...
    %     sum(sum(MSA2_recomb_history(ind,:)))
    % end

    % sumCOVdif = testrec(COV,COV_orig,recomb_COV,recomb_COV_orig)
    % sumCOVdif2 = ...
    %    testrec_2(glob_COV,glob_COV_orig,recomb_COV,recomb_COV_orig)
    %%%%

    % Evolution after the recombination.

    %%%% For testing:
    % COV_orig = COV;
    % mut_COV_orig = mut_COV;
    % cov_COV_orig = cov_COV;
    % glob_COV_orig = glob_COV;
    % mutcov_COV_orig = mutcov_COV;
    %%%%

% Here we cycle through each sub-branch.

a = find(b_ind3(:,1) == b_ind2(k,1));
b = find(b_ind3(:,2) == b_ind2(k,2));

    for i = a:b
        ind = b_ind3(i,1):b_ind3(i,2);

        mut_rate = randi(mut_rate_range_8);    
        
        for m = 1:mut_cycles
        mut_rate_1 = mut_rate + randi([nmrm,pmrm]);    
            
        [MSA2(ind,:),MSA2_history(ind,:),MSA2_mut_history(ind,:),...
            MSA2_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
            cov_COV(:,:,ind),glob_COV(:,:,ind),...
            mutcov_COV(:,:,ind)] = ...
            evolve_6s_3D_cov_hist_multiplets(MSA2(ind,:),...
            MSA2_history(ind,:),MSA2_mut_history(ind,:),...
            MSA2_cov_history(ind,:),COV(:,:,ind),mut_COV(:,:,ind),...
            cov_COV(:,:,ind),glob_COV(:,:,ind),...
            mutcov_COV(:,:,ind),...
            mut_rate_1,covar_vec,mut_prob_PD,...
            cov_prob_PD,cov_alphabet,mut_freq_PD);
        end
    end

    %%%% For testing:
    % for i = a:b
    %       ind = b_ind3(i,1):b_ind3(i,2);
    % dif_COV_HIST = ...
    %     sum(diag(sum(COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
    % dif_globCOV_HIST = ...
    %     sum(diag(sum(glob_COV(:,:,ind),3))) - sum(sum(MSA2_history(ind,:)))
    % dif_mutCOV_HIST = ...
    %     sum(diag(sum(mut_COV(:,:,ind),3))) - sum(sum(MSA2_mut_history(ind,:)))
    % dif_covCOV_HIST = ...
    %     sum(diag(sum(cov_COV(:,:,ind),3))) - sum(sum(MSA2_cov_history(ind,:)))
    % end
    
    % sumCOVdif = testmut(COV,COV_orig,mut_COV,mut_COV_orig,...
    % cov_COV,cov_COV_orig)
    % sumCOVdif2 = ...
    %     testmut_2(glob_COV,glob_COV_orig,mutcov_COV,mutcov_COV_orig )
    %%%%

end % end of nfinal_cycle

end % end of nbranches

    %%%% For testing:
    % dif_COV_HIST = ...
    % sum(diag(sum(COV(:,:,:),3))) - sum(sum(MSA2_history(:,:)))
    % dif_globCOV_HIST = ...
    % sum(diag(sum(glob_COV(:,:,:),3))) - sum(sum(MSA2_history(:,:)))
    % dif_mutCOV_HIST = ...
    % sum(diag(sum(mut_COV(:,:,:),3))) - ...
    % sum(sum(MSA2_mut_history(:,:)))
    % dif_covCOV_HIST = ...
    % sum(diag(sum(cov_COV(:,:,:),3))) - ...
    % sum(sum(MSA2_cov_history(:,:)))
    % dif_recombCOV_HIST = ...
    % sum(diag(sum(recomb_COV(:,:,:),3))) - ...
    % sum(sum(MSA2_recomb_history(:,:)))
    %%%%

end % end of n3rd_cycle

% Here we check that the sum of the partial COV's is equal to the global
% COV.

%%%% For testing:
% sum_COV = mut_COV + cov_COV + recomb_COV;
% sum_COV2 = mutcov_COV + recomb_COV;
% COVdif = COV - sum_COV;
% COVdif2 = glob_COV - sum_COV2;
% sumCOVdif = sum(COVdif(1:end))
% sumCOVdif2 = sum(COVdif2(1:end))

% First we back up all the coevolution matrices.

% COV_orig = COV;
% mut_COV_orig = mut_COV;
% cov_COV_orig = cov_COV;
% glob_COV_orig = glob_COV;
% mutcov_COV_orig = mutcov_COV;
% recomb_COV_orig = mutcov_COV;
%%%%

    
MSA_all = [ancestor(1,:);MSA0;MSA1;MSA2];

[rows,~] = size(MSA2);

% Here we first save everything, then we flatten all the coevolution 
% matrices.
MSA2_ALL(:,:,n) = MSA2;
MSA2_history_ALL(:,:,n) = MSA2_history;
MSA2_recomb_history_ALL(:,:,n) = MSA2_recomb_history;
MSA2_mut_history_ALL(:,:,n) = MSA2_mut_history;
MSA2_cov_history_ALL(:,:,n) = MSA2_cov_history;

COV_3D = COV;
mut_COV_3D = mut_COV;
cov_COV_3D = cov_COV;
glob_COV_3D = glob_COV;
mutcov_COV_3D = mutcov_COV;
recomb_COV_3D = recomb_COV;

clear COV mut_COV cov_COV glob_COV mutcov_COV recomb_COV

if SAVE_ALL_COV
COV_3D_ALL(:,:,:,n) = COV_3D;
cov_COV_3D_ALL(:,:,:,n) = cov_COV_3D;
mut_COV_3D_ALL(:,:,:,n) = mut_COV_3D;
recomb_COV_3D_ALL(:,:,:,n) = recomb_COV_3D;
end

COV = sum(COV_3D,3);
mut_COV = sum(mut_COV_3D,3);
cov_COV = sum(cov_COV_3D,3);
glob_COV = sum(glob_COV_3D,3);
mutcov_COV = sum(mutcov_COV_3D,3);
recomb_COV = sum(recomb_COV_3D,3);

%%%% For testing:
% And we check that all the sums are still holding

% sum_COV = mut_COV + cov_COV + recomb_COV;
% sum_COV2 = mutcov_COV + recomb_COV;
% COVdif = COV - sum_COV;
% COVdif2 = glob_COV - sum_COV2;
% sumCOVdif = sum(COVdif(1:end))
% sumCOVdif2 = sum(COVdif2(1:end))

% dif_COV_HIST = ...
%     sum(diag(COV)) - sum(sum(MSA3_history_sum))
% dif_COV_HIST = ...
%     sum(diag(glob_COV)) - sum(sum(MSA3_history_sum))
% dif_COV_HIST = ...
%     sum(diag(mut_COV)) - sum(sum(MSA3_mut_history_sum))
% dif_COV_HIST = ...
%     sum(diag(cov_COV)) - sum(sum(MSA3_cov_history_sum))
% dif_COV_HIST = ...
%     sum(diag(recomb_COV)) - sum(sum(MSA3_recomb_history_sum))
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selective types of flattening can be carried out to obtain phylogenetic
% information. For example, if we are interested only in sequences [20,73,
% 121,230] of the 10th run, we could flatten the 3rd dimension of the array 
% COV_3D_ALL corresponding to those sequences:
% COV_20_73_121_230 = sum(COV_3D_ALL(:,:,[20,73,121,230],10),3)
% Ultimately, we can look at the coevolution of just one sequence (eg.
% 223):
% COV_223 = COV_3D_ALL(:,:,223,10).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear up memory

clear *orig 
clear ncopies
clear sumCOVdif*

%--------------------------------------------------------------------------
% Here we calculate the emission probabilities of
% each branch of the experimental msa and of the evolved msa.

c_MSA2 = int2aa(MSA2);

if GAPS
MSA_branch_t_emissions = zeros(20,npos,nbranches);
MSA_branch_profile = zeros(21,npos,nbranches);

for k = 1:nbranches
    ind = b_ind2(k,1):b_ind2(k,2);

MSA_branch_hmm_model = nmsa_to_hmm_model(MSA2(ind,:));
MSA_branch_t_emissions(:,:,k) = MSA_branch_hmm_model.MatchEmission';
MSA_branch_profile(:,:,k) = seqprofile(c_MSA2(ind,:),'Gaps',gaps_count);
end
    
else
    
MSA_branch_t_emissions = zeros(20,npos,nbranches);
MSA_branch_profile = zeros(20,npos,nbranches);

for k = 1:nbranches
    ind = b_ind2(k,1):b_ind2(k,2);

MSA_branch_hmm_model = nmsa_to_hmm_model(MSA2(ind,:));
MSA_branch_t_emissions(:,:,k) = MSA_branch_hmm_model.MatchEmission';
MSA_branch_profile(:,:,k) = seqprofile(c_MSA2(ind,:),'Gaps','none');
end
end
%--------------------------------------------------------------------------        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if UNEVEN_BRANCHES
    
MSA_select = MSA2;

else    
% Here we select only a fraction of all the sequences in the tree. This is
% not recommended. It is better to simulate a smaller or larger tree as
% needed and then take all the sequences.

if isempty(select_msa)
    select = (1:nselect)';
elseif select_msa<rows
    select = randperm(rows)';
    select = select(1:select_msa,1);
else
    select = (1:nselect)';
end

select_orig_ind(select,n) = 1;
MSA_select_orig_ALL(:,:,n) = MSA2;
MSA_select = MSA2(select,:);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we update the information on the size of MSA_select.

[rows,~] = size(MSA_select);

if size(MSA_select_ALL,1)>rows
MSA_select_ALL(rows+1:nselect,:,:) = [];
end

% Then we store a copy of the current MSA_select in the array for the 
% entire run.

MSA_select_ALL(:,:,n) = MSA_select;

% s_MSA_select is the alignment with each sequence as a separate cell.
% c_MSA_select is the alignment in single letter selex format.

s_MSA_select = nmsa_to_smsa(MSA_select);  
c_MSA_select = int2aa(MSA_select);
  
% If printing the msa and tree as an ascii file.
% multialignwrite('s_MSA_select.aln',s_MSA_select); 
% multialignwrite('s_MSA_select.msf',s_MSA_select);
% phytreewrite('s_MSA_select.tree',MSA_tree);

% Here we calculate the correlations between the emission probabilities of
% the experimental msa and of the evolved msa.

MSA_select_hmm_model = nmsa_to_hmm_model(MSA_select);
emission_corr(n,:) = (diag(corr(REF_hmm_model.MatchEmission',...
    MSA_select_hmm_model.MatchEmission')))';
mean_emission_corr(n) = mean(emission_corr(n,:));

% t_emissions and profiles are already in format [20 x npos x nbranches],
% so we will not take the transpose.
for k = 1:nbranches
branch_emission_corr(k,:,n) = (diag(corr(REF_branch_t_emissions(:,:,k),...
    MSA_branch_t_emissions(:,:,k))))';
branch_profile_corr(k,:,n) = (diag(corr(REF_branch_profile(:,:,k),...
    MSA_branch_profile(:,:,k),'rows','pairwise')))';
mean_branch_emission_corr(k,n) = nanmean(branch_emission_corr(k,:,n),2);
mean_branch_profile_corr(k,n) = nanmean(branch_profile_corr(k,:,n),2);
end
    
% Here we calculate the distance or covariance matrix of the MSA so we can  
% compare the evolved msa with the experimental msa.
 
 if DISTANCE
 MSA_dist = seqpdist(c_MSA_select,'SquareForm',true,'Method','p-distance');
 MSA_cov_mat = 1-MSA_dist;
 else
 MSA_binmsa = nmsa_to_binmsa(MSA_select); 
 MSA_cov_mat = cov(MSA_binmsa',1);    
 end

 MSA_sim_score(n,1) = multidim_MI(MSA_cov_mat);
 if size(MSA_cov_mat_ALL,1)>rows
 MSA_cov_mat_ALL(rows+1:nselect,:,:) = [];
 MSA_cov_mat_ALL(:,rows+1:nselect,:) = []; 
 end
 MSA_cov_mat_ALL(:,:,n) = MSA_cov_mat;

% For example we can identify the clusters in the simulated MSA.

% [MSA_pc,MSA_ev] = spectral(MSA_cov_mat);
% MSA_Clusters = clusterdata(MSA_pc(:,1:DIMENSIONS),'maxclust',MAX_CLUSTERS);
% MSA_Clusters_ALL(:,n) = MSA_Clusters;

% Find the consensus sequence in the MSA. 

[MSA_s_consensus,MSA_s_consensus_score] = seqconsensus(s_MSA_select);
MSA_entropy(n,:) = Entropy(MSA_select);

for i=1:npos
for j=i:npos
    gap=find((MSA_select(:,i)~=25) & (MSA_select(:,j)~=25));
    MSA_joint_entropy(i,j,n) = JointEntropy([MSA_select(gap,i)...
        MSA_select(gap,j)]);
end
end

for j = 1:n
    MSA_joint_entropy(:,:,n) = ...
        MSA_joint_entropy(:,:,n) + MSA_joint_entropy(:,:,n)';
    MSA_joint_entropy_bk = MSA_joint_entropy;
for i = 1:npos;
    MSA_joint_entropy(i,i,n) = MSA_joint_entropy(i,i,n)/2;        
end
end

% We also find the relative entropy.

[MSA_profile,~] = seqprofile(s_MSA_select,'Gaps','none');
[~,MSA_rel_entropy(n,:)] = rel_entropy_nmsa(MSA_select,...
    bg_frequencies,MSA_profile);

% Calculate the corr. coeff. between the two rel. entropy vectors. 
% It should be high if everything went well.

corr_REF_MSA_rel_entropy(n) = corr(REF_rel_entropy',MSA_rel_entropy(n,:)');

% Here we find the consensus sequence and the consensus score.

MSA_consensus(n,:) = aa2int(MSA_s_consensus);
MSA_consensus_score(n,:) = MSA_s_consensus_score;

% "MSA_cons_anc_sim_score" is the percent identity between ancestor and 
% consensus sequence.

MSA_cons_anc_sim_score(n,1) = ...
    sum(MSA_consensus(n,:) == ancestor(1,:))/npos;

% Since we also calculate the consensus score for each position of the 
% consensus sequence we can compare the consensus score vector   
% "s_consensus_score" with that calculated from the consensus sequence 
% derived from the experimental msa "REF_s_consensus_score".
% The following is the corr. coeff. between the two consensus scores. 
% It should be high if everything went well.

corr_REF_MSA_s_consensus_score(n) = ...
    corr(REF_s_consensus_score',MSA_s_consensus_score');
 
clear i j h k m A B D V gap name
clear vec1 vec2 a b ind

% Check the execution time.
tElapsed = toc(evol_start_time);
fprintf('Tree calculation: %d minutes \n', tElapsed/60);

clear evol_start_time evol_cycle_time

%--------------------------------------------------------------------------
% 'Co-evolution' section. 
%--------------------------------------------------------------------------
% Here we calculate various forms of coevolution matrix and the fraction 
% of the 'ncov' top zscores in each coev matrix that corresponds to 
% covarying positions.
% COV is the real coevolution map of the entire MSA derived from the 
% recorded hystory of all mutations and recombinations carried out in the 
% 'evolution' section. Since we NaN the diagonal of COV we also make a
% backup copy of the original matrix.

% Here we bring back the original covar_vec.
calc_start_time = tic;
covar_vec_large = covar_vec;
ncov_large = size(covar_vec_large,1);
% covar_vec = sortrows(covar_vec_2,1);
% ncov = size(covar_vec,1);
covar_vec = cov_coord;
% Reorder covar_vec with the first elements of the pairs in ascending order.    
covar_vec_2 = sortrows(covar_vec,1);
covar_vec_2 = sort(covar_vec_2,2);
covar_vec_2 = sortrows(covar_vec_2,1);
covar_vec_2 = sortrows(covar_vec_2,2);
covar_vec = sortrows(covar_vec_2,1);
ncov = ncov_coord;

clear covar_vec_2

 COV = COV + COV';
 COV_bk = COV;
 recomb_COV = recomb_COV + recomb_COV';
 recomb_COV_bk = recomb_COV;
 mut_COV = mut_COV + mut_COV';
 mut_COV_bk = mut_COV;
 cov_COV = cov_COV + cov_COV';
 cov_COV_bk = cov_COV;
 glob_COV = glob_COV + glob_COV';
 glob_COV_bk = glob_COV;
 mutcov_COV = mutcov_COV + mutcov_COV';
 mutcov_COV_bk = mutcov_COV;
 
    for i = 1:npos;
        COV(i,i) = NaN;
        COV_bk(i,i) = COV_bk(i,i)/2;        
        recomb_COV(i,i) = NaN;
        recomb_COV_bk(i,i) = recomb_COV_bk(i,i)/2;        
        mut_COV(i,i) = NaN;
        mut_COV_bk(i,i) = mut_COV_bk(i,i)/2;        
        cov_COV(i,i) = NaN;
        cov_COV_bk(i,i) = cov_COV_bk(i,i)/2;        
        glob_COV(i,i) = NaN;
        glob_COV_bk(i,i) = glob_COV_bk(i,i)/2;        
        mutcov_COV(i,i) = NaN;
        mutcov_COV_bk(i,i) = mutcov_COV_bk(i,i)/2;        
    end
    
 COV_ntz = nantozero(COV);
 COV_ALL(:,:,n) = COV;
 COV_bk_ALL(:,:,n) = COV_bk;
 COV_ntz_ALL(:,:,n) = COV_ntz;
 recomb_COV_ntz = nantozero(recomb_COV);
 recomb_COV_ALL(:,:,n) = recomb_COV;
 recomb_COV_bk_ALL(:,:,n) = recomb_COV_bk;
 recomb_COV_ntz_ALL(:,:,n) = recomb_COV_ntz;
 mut_COV_ntz = nantozero(mut_COV);
 mut_COV_ALL(:,:,n) = mut_COV;
 mut_COV_bk_ALL(:,:,n) = mut_COV_bk;
 mut_COV_ntz_ALL(:,:,n) = mut_COV_ntz;
 cov_COV_ntz = nantozero(cov_COV);
 cov_COV_ALL(:,:,n) = cov_COV;
 cov_COV_bk_ALL(:,:,n) = cov_COV_bk;
 cov_COV_ntz_ALL(:,:,n) = cov_COV_ntz;
 glob_COV_ntz = nantozero(glob_COV);
 glob_COV_ALL(:,:,n) = glob_COV;
 glob_COV_bk_ALL(:,:,n) = glob_COV_bk;
 glob_COV_ntz_ALL(:,:,n) = glob_COV_ntz;
 mutcov_COV_ntz = nantozero(mutcov_COV);
 mutcov_COV_ALL(:,:,n) = mutcov_COV;
 mutcov_COV_bk_ALL(:,:,n) = mutcov_COV_bk;
 mutcov_COV_ntz_ALL(:,:,n) = mutcov_COV_ntz;

% One last check. The following should all be 0 if everything was OK during 
% the evolution 

sum_COV = mut_COV_bk + cov_COV_bk + recomb_COV_bk;
sum_COV2 = mutcov_COV_bk + recomb_COV_bk;
COVdif = COV_bk - sum_COV;
COVdif2 = glob_COV_bk - sum_COV2;
sumCOVdif = sum(COVdif(1:end));
sumCOVdif2 = sum(COVdif2(1:end));

dif_COV_HIST = ...
 sum(diag(COV_bk_ALL(:,:,n))) - sum(sum(MSA2_history_ALL(:,:,n)));
dif_glob_COV_HIST = ...
 sum(diag(glob_COV_bk_ALL(:,:,n))) - sum(sum(MSA2_history_ALL(:,:,n)));
dif_mut_COV_HIST = ...
 sum(diag(mut_COV_bk_ALL(:,:,n))) - sum(sum(MSA2_mut_history_ALL(:,:,n)));
dif_cov_COV_HIST = ...
 sum(diag(cov_COV_bk_ALL(:,:,n))) - sum(sum(MSA2_cov_history_ALL(:,:,n)));
dif_recomb_COV_HIST = sum(diag(recomb_COV_bk_ALL(:,:,n))) - ...
 sum(sum(MSA2_recomb_history_ALL(:,:,n)));

disp(...
'The following should all be 0 if everything was OK during the evolution:');
fprintf('sumCOVdif = %d \n', sumCOVdif);
fprintf('sumCOVdif2 = %d \n', sumCOVdif2);
fprintf('dif_COV_HIST = %d \n', dif_COV_HIST);
fprintf('dif_glob_COV_HIST = %d \n', dif_glob_COV_HIST);
fprintf('dif_mut_COV_HIST = %d \n', dif_mut_COV_HIST);
fprintf('dif_cov_COV_HIST = %d \n', dif_cov_COV_HIST);
fprintf('dif_recomb_COV_HIST = %d \n', dif_recomb_COV_HIST);

% Free up memory:
clear sum_COV sum_COV2 COVdif COVdif2
clear dif_COV_HIST dif_glob_COV_HIST dif_mut_COV_HIST
clear dif_cov_COV_HIST dif_recomb_COV_HIST

 % Output from "coev_stats:
 % fcov_# 1st field:
 % Percentage of the top 10 zscores in the COV matrix that belong to the
 % covar_vec.
 % fcov_# 2nd field:
 % Percentage of the top N zscores (where N is the length of
 % covar_vec) that belong to the covar_vec.
 % fcov_# 3rd field:
 % Percentage of the top 100 zscores in the COV matrix that belong to the
 % covar_vec.
 % "cov_zscore_#" is the zscore in the COV matrix for each of the covarying
 % positions in covar_vec. Only the score is reported
 % "covar_vec_zscore_#" is the zscore in the COV matrix for each of the 
 % covarying positions in covar_vec. Both the covarying positions (the 
 % covar_vec) and the scores are reported.

 [fcov_COV(n,:),cov_zscore_COV(n,:),...
     covar_vec_zscore_COV(:,:,n)] = ... 
     coev_stats_2(COV,ncov,covar_vec);
 [fcov_recomb_COV(n,:),cov_zscore_recomb_COV(n,:),...
     covar_vec_zscore_recomb_COV(:,:,n)] = ... 
     coev_stats_2(recomb_COV,ncov,covar_vec);
 [fcov_mut_COV(n,:),cov_zscore_mut_COV(n,:),...
     covar_vec_zscore_mut_COV(:,:,n)] = ... 
     coev_stats_2(mut_COV,ncov,covar_vec);
 [fcov_cov_COV(n,:),cov_zscore_cov_COV(n,:),...
     covar_vec_zscore_cov_COV(:,:,n)] = ... 
     coev_stats_2(cov_COV,ncov,covar_vec);
 [fcov_glob_COV(n,:),cov_zscore_glob_COV(n,:),...
     covar_vec_zscore_glob_COV(:,:,n)] = ... 
     coev_stats_2(glob_COV,ncov,covar_vec);
%--------------------------------------------------------------------------
% End of the 'co-evolution' section
%--------------------------------------------------------------------------
 
% Check the execution time.

tElapsed = toc(calc_start_time);
fprintf('Coevolution calculation: %d minutes \n', tElapsed/60);
clear calc_start_time cov_calc_time

end % End of the end_cycle loop

%% 'statistics' section.

% This section runs only if the 'evolution' section is run for multiple
% cycles (end_cycle > 1; recommended end_cycle = 100);
% Here we get the parameters of the two distributions;

if  end_cycle > 1
stat_fcov_COV = [mean(fcov_COV);std(fcov_COV)];
stat_fcov_cov_COV = [mean(fcov_cov_COV);std(fcov_cov_COV)];
stat_fcov_mut_COV = [mean(fcov_mut_COV);std(fcov_mut_COV)];
stat_fcov_recomb_COV = [mean(fcov_recomb_COV);std(fcov_recomb_COV)];
stat_fcov_glob_COV = [mean(fcov_glob_COV);std(fcov_glob_COV)];

% Here we look at the relationship between the consensus sequence in the
% simulated MSA and the ancestor.
    
MSA_cons_anc_gevPD = fitdist(MSA_cons_anc_sim_score,'normal');
MSA_cons_anc_gevPD_params(1,1) = mean(MSA_cons_anc_sim_score);
MSA_cons_anc_gevPD_params(1,2) = std(MSA_cons_anc_sim_score);
MSA_cons_anc_gevPD_params(1,3:4) = MSA_cons_anc_gevPD.Params;
MSA_sim_score_gevPD = fitdist(MSA_sim_score,'normal');
MSA_sim_score_gevPD_params(1,1) = mean(MSA_sim_score);
MSA_sim_score_gevPD_params(1,2) = std(MSA_sim_score);
MSA_sim_score_gevPD_params(1,3:4) = MSA_sim_score_gevPD.Params;

% The following is a logical array which for each position of each sequence
% has 1 when ancestor and consensus are identical and 0 when they are not.

con_anc_int_ind = false(end_cycle,npos);
for k = 1:end_cycle
con_anc_int_ind(k,:) = ancestor_ALL(k,:) == MSA_consensus(k,:);
end

mean_MSA_entropy = mean(MSA_entropy);
mean_MSA_rel_entropy = mean(MSA_rel_entropy);
mean_con_anc_int_ind = mean(con_anc_int_ind);
mean_MSA_cov_mat = mean(MSA_cov_mat_ALL,3);

% The following is the correlation between the mean entropy for each
% position of the msa and the mean identity between consensus and ancestor
% at each position of the msa

corr_AC_MSA_entropy = corr(mean_MSA_entropy',mean_con_anc_int_ind');
corr_AC_MSA_rel_entropy = corr(mean_MSA_rel_entropy',mean_con_anc_int_ind');
corr_REFhmm_MSA_entropy = corr(REF_hmm_rel_entropy',mean_MSA_entropy');
corr_REFhmm_MSA_rel_entropy = corr(REF_hmm_rel_entropy',mean_MSA_rel_entropy');

mean_cov_zscore_COV = mean(cov_zscore_COV,2);
mean_cov_zscore_cov_COV = mean(cov_zscore_cov_COV,2);

end

%% 'save' section
if SAVEFILE
    save(savefilename);
end
%--------------------------------------------------------------------------

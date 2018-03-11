% This script sets a number of general parameters that define how the
% MSAvolve run will be carried out. First we store the program version.
PROGRAM = 'MSAvolve v.3.0a';
%--------------------------------------------------------------------------
SAVEFILE = 1;   % Default is 1: we save the workspace at the end of the run.
savefilename = 'my_protein_family_MSAvolve_run';
%--------------------------------------------------------------------------
% If either of the two following fields is set to [] a new seed is
% initialized for the random number generator.
RANDOM_SEED_FILE = [];  % Default is []. File name must be in quotes ''.
RANDOM_SEED_VAR = [];   % Default is []. Variable name must be in quotes ''.
%--------------------------------------------------------------------------
RUN_NOW = 1;    % If set to 1 MSAvolve will run with the parameters set in 
                % this file. If you want to run MSAvolve in the safer cell
                % mode, step-by-step, set RUN_NOW to 0.
CHECK_BRANCHES = 0; % Default is 0. If 1 control is passed to the keyboard. 
                    % If you are satisfied with the assignment of branches
% based on clusters, of covarions, and if the histogram of the covariance
% matrix is smooth, without peaks at high covariance, you can go ahead.
%
% To continue the program type 'return' and then enter.
%
% If you are not satisfied, without leaving the program, you can check the
% composition of the clusters from the command line. For example you can
% get the msa's corresponding to the clusters: REF_smsa_clust1 =
% REF_smsa(REF_clust1); REF_smsa_clust2 = REF_smsa(REF_clust2); ... and
% print them out in clustal, msf, or fasta format.
% multialignwrite('REF_cluster_1.aln',REF_smsa_clust1);
% multialignwrite('REF_cluster_2.aln',REF_smsa_clust2);
% fastawrite('REF_cluster_1.faln',REF_smsa_clust1);
% fastawrite('REF_cluster_2.faln',REF_smsa_clust2);
%
% If the histogram shows a peak at high covariance it means there is a
% large number of sequences that are too similar to each other. In that
% case, without leaving the program, you can remove sequences that are more
% than 90% identical to any other sequence with the function
% 'trim_nmsa_by_threshold': [REF_trimmed_nmsa,REF_trimmed_smsa] =
% trim_nmsa_by_threshold(REF_nmsa,0.9);
% fastawrite('trimmed90_msa.faln',REF_trimmed_smsa);
%
% To quit the program type 'dbquit' and then enter.
%
% Alternatively, you can trim your original msa with 't_coffee': "t_coffee
% -other_pg seq_reformat -in complete.aln -action +trim
%   _aln_%%90_ +upper > trimmed.aln"
%--------------------------------------------------------------------------
MAKE_MSA = 0;   % If set to 1 a reference MSA is generated randomly using 
                % the dimensions (no_seq,no_pos,no_branches) stated below.
                % Make sure the number of positions is consistent with the
                % crossover points defined in the recombination section.
no_seq = 100;   % However, this is not the best way to run MSAvolve. Much
no_pos = 100;   % better to use an external experimental MSA as reference.
no_branches = 3; 

% If we are going to read an external MSA we declare its name here.

ext_MSA = 'my_protein_family.faln';

fasta = 1;  % If 0 the alignment is Clustal or MSF, and the suffix (.aln 
            % or .msf) defines the type of MSA used as imput. If 1 the
            % alignment is in fasta format.
%--------------------------------------------------------------------------
% We repeat the tree building more than once only if we intend to
% accumulate various kinds of statistics. REMEMBER: the memory requirements
% will increase proportionally to the number of building cycles. Here we
% also decide if it is a completely new run or if it is the restart of a
% run that crushed. If it is a restart, set NEW_RUN to 0 and change
% 'start_cycle' to the last saved cycle of the interrupted run. If you
% restart you cannot increase the original 'end_cycle' value because all
% the arrays are already declared. However, you can change it to a smaller
% number than the original. This section is also very useful if you want to
% repeat one run in particular for the purpose of getting the 4D COV
% matrices. In this case, set NEW_RUN = 0 and set both 'start_cycle' and
% 'end_cycle' to the number of the run you want to repeat. Then set
% SAVE_ALL_COV = 1, which should otherwise always be equal to 0 to avoid
% huge memory loads. PLEASE NOTE: If you don't want to run any evolution
% cycle, but just want to get the random msa's, set 'end_cycle' to 0 and
% 'random_msa_no' to the desired number of random msa's.

NEW_RUN = 1;
start_cycle = 1;
end_cycle = 10;

%--------------------------------------------------------------------------
% Choice of the ancestor: based on background probabilities or on the
% emission probabilities of the hmm model of the experimental msa.
ANCESTOR_BG_PROB = 0;      % Default is 0: ancestor based on hmm emissions.

%--------------------------------------------------------------------------
% Choice of a reference random msa.
RANDOM_MSA = 0;     % Default is 0: if set to 1 random msa's are generated 
                    % with the same criteria adopted to generate the
                    % ancestor. They have the same statistical properties
                    % of the final simulated msa, but they were not
                    % obtained by evolution.
random_msa_no = 1;  % Default is 'end_cycle', the number of evolved 
                            % msa's, but can be set to any number. The
                            % number of sequences in the msa is set
                            % automatically equal to the number of
                            % sequences in the simulated msa.
random_cycles = 2;  % Default is 2: larger number increase the overall 
                    % relative entropy (similarity between sequences). If
                    % the default value produces a higher or lower overall
                    % relative entropy with respect to the reference msa,
                    % usually 1 or 3 random cycles will produce an almost
                    % perfect match. Run time increases linearly with the
                    % number of cycles.
random_cov_cycles = 2;  % Default is 2: larger number increase the chance
                        % that more frequent combinations of A,B at ij
                        % position will be selected. The effect is similar
                        % to that of 'random_cycles', but works only on
                        % positions that covary.
WITH_COVARIONS = 1; % Default is 1: random msa's include covariation based 
                    % on the choice of covarions determined by the
                    % selections in the covariation section. If set to 0
                    % any covariation in the random msa's is purely by
                    % chance, and thus this option is very useful for
                    % detecting the background stochastic coevolution.
%--------------------------------------------------------------------------
% Parameters for tree building.
DISTANCE = 3;   % Default is 0. Phylogenetic trees and clustering are based
                % on the covariance matrix of the binary msa. If set to 1 a
                % distance matrix is calculated with the 'p-distance'
                % method. In this case, if the label 'covariance' appears
                % in any of the plot axes or titles, what is really meant
                % is 'p-distance'. If set to 2 a distance matrix is
                % calculated with the 'alignment-score' method. In this
                % case, if the label 'covariance' appears in any of the
                % plot axes or titles, what is really meant is '1 -
                % distance' by alignment score. If set to 3 a distance
                % matrix is calculated with the 'Jukes-Cantor' method. In
                % this case, if the label 'covariance' appears in any of
                % the plot axes or titles, what is really meant is '1 -
                % distance' by the Jukes-Cantor method. 
SCREE = 0; % Default is 0. If set to 1 a 'scree' plot of the eigenvalues of 
           % the covariance matrix is plotted. It may be useful to decide
           % the number of dimensions to use.
DIMENSIONS = 3; % Default is 3. We use the first 3 principal components of the
                % distance or covariance matrix to cluster the sequences in
                % the msa. Obviously, the maximum that can be visualized is
                % 3, although more dimensions can be used. However, if you
                % use more than three dimensions the cluster assignment may
                % appear odd because you don't see the relationship to the
                % 3+ dimensions. We recommend using 3 dimensions and at
                % most 10 clusters.
ICA = 0; % default is 0. If fastICA is installed, setting ICA to 1 will 
         % trigger the use of Independent Component Analysis to cluster the
         % sequences.
MAX_CLUSTERS = 10; % There is no default, it should be decided by looking at 
                  % the spectral analysis of the distance or covariance
                  % matrix. BE CAREFUL: too many clusters increase the
                  % memory requirements and decrease the validity of the
                  % statistical analysis. REMEMBER: the number of clusters
                  % becomes the number of branches in the simulated msa.
cluster_distance = 'euclidean'; % any distance metric allowed by Matlab 'pdist'.
cluster_linkage = 'ward'; % any distance method allowed by Matlab 'linkage':
                              % useful combinations: [euclidean/ward],
                              % [chebychev/complete],[minkowski/complete].
TREE_ONLY = 0;  % Default is 0. Setting this value to 1 is only useful to 
                % see the phylogenetic tree, but is not recommended to
                % assign the clusters of sequences. The methods used for
                % the tree are defined in the lines below. If changed to 1,
                % remember to set it back to 0 for optimal cluster assignment. 
distance_method = 'alignment-score'; % Options available for Matlab 
                                % seqpdist:'p-distance','Jukes-Cantor',
                                % 'alignment-score', if scoring of gaps is
                                % required. If gap scoring is not required
                                % other options are 'Poisson','Gamma'.
tree_method = 'seqneighjoin'; % Default is 'seqneighjoin'(Neighbor-joining) 
                % Other options are those provided by Matlab 'seqlinkage':
                % 'single' (= Nearest Distance), 'complete' (= Furthest
                % distance), 'average' (= UPGMA), 'weighted' (= WPGMA),
                % 'centroid' (= UPGMC), 'median' (= WPGMC).
neighjoin_method = 'equivar'; % Neighbor-joining method: default is the
                % the 'equivar' option (Studier, J.A., Keppler, K.J.
                % Molecular Biology and Evolution 5(6) 729–731, 1988).
                % Other options are 'firstorder' (Gascuel, O. Molecular
                % Biology and Evolution 14 685–695, 1997) and 'average'
                % (Saitou, N., and Nei, M. Molecular Biology and Evolution
                % 4(4), 406–425, 1987).
criterion = 'silhouette'; % Criterion used for clustering the branches of  
                          % the tree. Any choice available for the Matlab
                          % function cluster(phytree):'maximum','gain',
                          % 'silhouette','ratio','median','average'. Check 
                          % the Matlab documentation for the meaning of 
                          % each choice.
UNEVEN_BRANCHES = 1;    % Default is 1. We select sequences in the 
                        % simulated MSA in order to reproduce the size of
                        % the original main branches. If set to 0 the main
                        % branches of the simulated MSA will have equal
                        % number of sequences and the total number of
                        % sequences in the final simulated msa is given by
                        % the product nbranches*ncopies_1*ncopies_2.
msa_size = [];   % Default is [], in which case the size of the simulated msa
                 % is as close as possible to that of the experimental msa.
                 % If set to a number, the size of the simulated msa is
                 % made as close to that number keeping the same ratio
                 % between branches and between tree expansions. This
                 % selection affects only the UNEVEN choice, and is
                 % recommended if the experimental msa is very large and we
                 % want a much smaller simulated msa.
GAPS = 1;   % Default is 1, in which case mutation probabilities are based 
            % on the profile of the reference msa calculated including
            % gaps. In this case it is necessary to set the flag
            % 'gaps_count' to one of the two allowed values: 'all' or
            % 'noflanks'. If set to 0 gaps are not counted (not
            % recommended) and are not considered symbols that can appear
            % during the evolution.
gaps_count = 'all';  % Default is 'all': all gaps are counted. If set to 
                     % 'noflanks' all gaps are counted except those at the
                     % flanks of every sequence. The 'noflanks' option is
                     % very risky because some branches may have only gaps
                     % at certain positions, and this will crash the
                     % program when attempting to determine the probability
                     % distribution for different aa's at those positions.
                     % The same problem might occur if GAPS is set to 0. In
                     % practice the option 'noflanks' works only if very
                     % few clusters are sought.
%--------------------------------------------------------------------------
% Rules for branch expansion and the selection of the final msa.
% 'ncopies_1' and 'ncopies_2' get rescaled if UNEVEN_BRANCHES is set to 1
% in order to obtain a final MSA of approximately the same size as the
% reference msa.
ncopies_1 = 10;     % 2nd node expansion [10]****.
ncopies_2 = 10;    % 3rd node expansion [15]****.
select_msa = [];   % Default is [], in which case all sequences, which are 
                   % the product nbranches*ncopies_1*ncopies_2, are taken.
                   % If we are using equal branches, we can use
                   % 'select_msa' to select a given number of sequences
                   % randomly. However, this is not recommended: it is
                   % better/faster to simulate a smaller or larger tree as
                   % needed and then take all the sequences. In any case
                   % 'select_msa' cannot be larger than
                   % nbranches*ncopies_1*ncopies_2, where nbranches is
                   % either defined by us in the MAKE_MSA section or is
                   % determined automatically from the number of clusters
                   % in experimental msa.
%--------------------------------------------------------------------------
% Each mutation and recombination cycle internal to one of the three main
% levels of evolution can be repeated 'mut_cycles' and 'rec_cycles' times.
% In this way the balance between point mutations and recombination can be
% changed at will.
mut_cycles = 11;     % [11]
rec_cycles = 1;     % [1]
%----------------NODE 1----------------------------------------------------
mut_rate_range_1 = [6,10]; % Rate for the ancestor 'nbranches' before
                           % recombination [6,10].
nhrt_range_1 = [0,2];      % Recombination between the 'nbranches' of the 
                           % ancestors [0,2].
mut_rate_range_2 = [6,10]; % Rate for the ancestor 'nbranches' after 
                           % recombination [6,10].
%----------------NODE 2----------------------------------------------------
mut_rate_range_3 = [4,8];  % Rate for the 'nbranches' each with 'ncopies_1'
                           % [4,8].
n2nd_cycle = 2;            % Number of times only the following two steps 
                           % are repeated [3]****.
nhrt_range_2 = [1,3];      % Recombination between the 'nbranches' [1,3].
                           % nhrt_range_1 and nhrt_range_2 act at the same
                           % level producing recombination between the main
                           % branches. However 'nhrt_range_2' acts after
                           % the first round of expansion of the data set
                           % and thus affects a larger number of sequences.
mut_rate_range_4 = [4,8];  % Rate for the 'nbranches' after recombination 
                           % [3,6]
nhrt_range_3 = [1,4];      % Recombination inside each main branch [1,4]. 
mut_rate_range_5 = [4,8];  % Rate for the 'nbranches' after recombination
                           % [2,6].
%----------------NODE 3----------------------------------------------------
n3rd_cycle = 1;            % Number of times all the following steps are 
                           % repeated [1]****.
mut_rate_range_6 = [2,6];  % Rate for the 'ncopies_2' expansion of each 
                           % 'ncopies_1' [2,6].
nhrt_range_4 = [1,3];      % Recombination inside each main branch [1,3].
                           % nhrt_range_3 and nhrt_range_4 act at the same
                           % level producing recombination only inside each
                           % of the main branches, but nhrt_range_4 acts on
                           % a much larger number of sequences.
mut_rate_range_7 = [1,4];  % Rate for the 'ncopies_2' expansion of each 
                           %'ncopies_1' after recombination [1,4]
nfinal_cycle = 1;          % Number of times only the following two steps  
                           % are repeated [1]****.
nhrt_range_5 = [2,6];      % Recombination inside each expanded branch [0,2].                             
mut_rate_range_8 = [1,4];  % Rate for the 'ncopies_2' expansion of each 
                           % 'ncopies_1' after recombination [1,4].
nmrm = -1;                 % Negative and positive mutation rate modifiers
pmrm = 2;                  % in the last cycles of mutations [-1,2]. Higher 
                           % absolute values of these numbers produces
                           % larger differences between the leaves of the
                           % tree. However,'nmrm' cannot be larger in
                           % absolute value than min(mut_rate_range_8).
%--------------------------------------------------------------------------
% Rules for recombination.
nzones = 8;       % This is the number of blocks in which the sequence is
                   % divided. One or more blocks can be transferred at one
                   % time from one to one or more sequences. Recombination
                   % blocks are separated by 'crossover points'.
RAND_RECOMB_ZONES = 0; % Default is 0. Recombination zones are assigned 
                       % manually possibly using the SCHEMA algorithm. If
                       % set to 1 the recombination zones are assigned
                       % randomly at each evolution trial. If you want to
                       % assign randomly only once for all the evolution
                       % cycles, leave the default '0' and set the
                       % following to 1.
ONE_RAND_RECOMB_ZONES = 0; % Default is 0. If set to 1 the recombination
                           % zones are assigned randomly only once for all
                           % evolution trials.
RECOMB_SCALING = 0; % Default is 0. If set to 1 the level of recombination 
                    % in each branch is scaled according to the overall
                    % similarity in that branch. CAREFUL: '1' will magnify
                    % the effect of recombination or lack thereoff.
% Crossover points can be obtained from the output of RASSP inside the
% SCHEMA package or from the background regions of the ZPX/ZPX2 matrix of
% the experimental msa. Alternatively you can try the automatic
% identification of crossover points using the function:
% [~,~,~,~,crossover_points] = get_nmsa_recomb_points(REF_nmsa,
% 0.9985,0.01); In either case make sure the number of zones defined above
% corresponds to the number of crossover point. For example if you define
% crossovers_points = [89 193 216 259] there will be 5 zones and 4
% crossover points represented inside MSAvolve by the variable 'hrf'. hrf =
% [1 89 193 216 259 npos]. IMPORTANT: If crossover_points = [] there is no
% recombination.

crossover_points = [30 68 100 118 139 175 205]; % at .9999
% crossover_points = [];
%--------------------------------------------------------------------------
% Rules for covariation. Make sure selections are consistent: for example
% if RAND_COVARIONS = 1 then SET_COVARIONS = 0, and viceversa. Only one of
% ENTROPY, REL_ENTROPY, METHOD, and VECTOR should be set to 1.

fcov =  15;  % If 'fcov = 15' 30% of all positions will be set as covarions,
            % that is 15% covarying with another 15%.
             
RAND_COVARIONS = 0;     % Default is 0. If set to 1 covarions are
                        % assigned randomly at each evolution trial. In
                        % this case overlapping pairs (multiples) are
                        % automatically removed.

SET_COVARIONS = 1;  % If set to  1 covarions are assigned based on 
                    % entropy value or the result of some method to detect
                    % coevolution, or a vector of covarions.
                    
ENTROPY = 0;        % If set to 1 covarions are assigned based on entropy. 
REL_ENTROPY = 1;    % If set to 1 covarions are assigned based on relative 
                    % entropy.
cov_gaps = 1; % Defult is 0, covarions are picked only in columns where there
              % are no gaps. If 1 covarions can be in positions of the msa
              % where there are some gaps. BE CAREFUL: we recommend chosing
              % 1 only if the entropy/rel_entropy range is set to 'mid' or
              % to some mixture including 'mid'. DO NOT use with the
              % combinations (ENTROPY = 1 ; cov_entropy range = high) and
              % (REL_ENTROPY = 1 ; cov_entropy range = low).
cov_entropy_range = 'mid';   % Choices for the entropy or relative entropy  
                             % range are: 'low','mid','high', and random
                             % mixtures 'low-mid','low-high','mid-high',
                             % 'low-mid-high', 'mix'. 'mix' randomizes a
                             % little more than 'low-mid-high'. Remember:
                             % 'low entropy = high relative entropy', so
                             % the meaning of this selection depends on
                             % whether you choose to use entropy or
                             % relative entropy to select the covarions.
                             
METHOD = 0;          % If set to 0 the value of 'cov_method' is ignored. 
                     % If set to 1 covarions are selected using one of the
                     % available methods to detect coevolution in an msa.
cov_method = 'ZPX2';    % The recommended method for this task is 'ZPX2'.  
                        % All methods can be accessed independently using
                        % the function: [coev_mat,covar_vec] =
                        % get_nmsa_covar_vec(nmsa,fcov,method), which will
                        % output the coevolution matrix and a vector
                        % containing a number of covarions equal to the
                        % percentage (e.g. fcov = 15) of the sequence
                        % positions.
                     
VECTOR = 0;     % If set to 0, the value of 'covarions' is ignored. 
                % If set to 1 we have to provide a vector of covarions.
                % Useful if we want to test the effect of triples or
                % multiples of covarions. The syntax of the covarions
                % vector is such that entries can appear in multiple ways
                % to produce multiples. For example they can be introduced
                % as overlapping pairs = [23 45;23 70;70 67] or directly as
                % multiplets = [23 45 70; 67 123 180]. If the multiplets
                % are entered as mixed pairs with higher orders, all the
                % entries of the vector must have the same size by padding
                % the entries with zeros. For example: = [72 120 0 0;120
                % 130 142 0;123 120 72 34]. If there are overlapping
                % entries the combination of those entries are also
                % calculated if MULTIPLETS is turned on. In the end all
                % possible covarying pairs and multiplets are stored in a
                % variable called 'covar_vec_large'. You should check this
                % variable to really know what was covarying during the
                % simulation.
MULTIPLETS = 0;     % Default is 0. If set to 1 a different calculation is  
                    % used in an attempt to give a better approximation of
                    % the covariation within multiplets. The MULTIPLETS
                    % option is still experimental: we recommend to turn it
                    % off (= 0) unless either METHOD or VECTOR are turned
                    % on and there are multiplets in either one.
multiples = 20;     % This is the maximum number of element in a multiple 
                    % to be considered. It should be set to some number
                    % slightly larger than the total number of unique
                    % positions in the 'covarions' vector. You should
                    % increase this number if you get an error about
                    % mismatched arrays.
covarions = [50,72;50,123;50,177;50,181;50,230;72,123;72,177]; 
%--------------------------------------------------------------------------
SAVE_ALL_COV = 0;       % Default is 0. If set to 1 memory requirements  
                        % become huge for data sets with a large number of
                        % cycles. This variable should be set to 1 only if
                        % you are doing or repeating a single run to get
                        % the 4-dimensional COV matrices. These 4D matrices
                        % are useful only if you are interested in a
                        % phylogenetic analysis of the covariation.
%--------------------------------------------------------------------------
% The section below should not require any changes unless you want to
% override or change the effects of some of the decisions made in the
% previous sections.
%--------------------------------------------------------------------------
% The next line safeguards against misuse of the SAVE_ALL_COV variable.
if (end_cycle-start_cycle) > 0
    SAVE_ALL_COV = 0;
end
%--------------------------------------------------------------------------
% We are done! Now we can run MSAvolve.
if RUN_NOW
    run MSAvolve_v3_0
end
%--------------------------------------------------------------------------


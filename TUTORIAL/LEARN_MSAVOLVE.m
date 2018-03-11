% This is a tutorial that introduces the basic use of MSAvolve
% to simulate the evolution of a protein and derive a msa that
% mimics the msa of a given protein family. MSAvolve requires three Matlab
% toolboxes, Bioinformatics, Curve Fitting, and Statistics.
%
%--------------------------------------------------------------------------
% IMPORTANT!!!
% Each time you open one of the original functions or scripts of the
% toolbox to study it, immediately save it in your working directory with a
% local name. If you modify one of the native functions/scripts and you
% made a mistake, the Toolbox will not work properly any longer.
%--------------------------------------------------------------------------
%
%%
% Start with setting the path:

addpath(genpath('../../MSAVOLVE_v3.0a'));
% The packages GREMLIN and gplmDCA are inconsistent with some MSAvolve
% functions. We will reinstall when necessary.
rmpath(genpath('../../MSAVOLVE_v3.0a/COEVOLUTION_METHODS/FUNCTIONS/GREMLIN'));
rmpath(genpath('../../MSAVOLVE_v3.0a/COEVOLUTION_METHODS/FUNCTIONS/gplm_DCA_asymm'));

%%
% A small MSA in fasta format (ArsC_95.faln) is provided as a test set. We
% will start running MSAvolve on this test set. First we open in the editor
% the script RUN_MSAVOLVE_ARSC.m, which drives all the operations carried
% out by MSAvolve. There is a template for this script in the directory
% MSAVOLVE, which contains also the program itself. The template can be
% adapted to analyze other msa's.

open RUN_MSAvolve_ARSC

% The script describes succintly the meaning of the various variables. It's
% a good idea to read the explanation of each variable/flag before running
% the script. The meaning and correct use of each flag will become
% progressively more evident as you use MSAvolve with multiple msa's.
%
% Make sure the following variables are set as follows (use the binocular 
% icon in the action bar to find them):
%
% CHECK_BRANCHES = 1;
% end cycle = 1;
% DIMENSIONS = 3;
% MAX_CLUSTERS = 3;
% TREE_ONLY = 0;
%
% This means we will simulate only 1 msa, and we want to identify 3
% clusters of sequences (these are going to be the branches of the
% evolutionary tree) using 3 dimensions of the covariance matrix and
% ignoring the phylogenetic tree. Even before starting the simulation, the
% execution will stop giving us some information on the properties of the
% msa and its possible partition in branches. Execute the script by
% clicking on the green triangle in the upper bar of the editor, or
% alternatively:

run RUN_MSAvolve_ARSC

% A 4-panel figure opens: the top-left panel shows a histogram plot of the
% sequence covariance matrix. If you see only one symmetric peak, most 
% likely it means that there is only one population of sequences. The
% top-right panel shows a heat map of the sequence covariance matrix that
% allows a better visualization of the various groups of sequences. The
% bottom-left panel shows a cluster fragmentation of the msa based on a
% spectral analysis of the covariance matrix. Sequence clusters are shown
% with different colors. You can rotate the cluster box by clicking on the 
% round arrow button in the figure bar. The bottom-right panel shows the 
% relative entropy ( = conservation degree) of the msa and the positions of
% the covarions picked on the basis of the variable:
%
% cov_entropy_range = 'mid';
%
% See what happens if you change the parameters for cluster assignments.
% At this point the command windows shows a prompt K>>. It means Matlab is
% waiting for you to proceed with the simulation or quit. We will quit.
% In the command window type:

dbquit;

% and at the next prompt:

close all; 

% then change the number of clusters sought in the spectral analysis of the
% covariance matrix:
%
% MAX_CLUSTERS = 8; (then save)
%
% and execute again clicking on the green triangle. It is a good idea to
% clear everything in the workspace to avoid conflict between arrays that
% may be set differently with different values of various flags. You can
% type in the command window:

clear
run RUN_MSAvolve_ARSC

% It looks decent; 8 clusters are recognized. Perhaps we could have used 
% fewer clusters, and, in general it is better not to look for more than 
% 10 clusters.
%
% Before proceeding with the simulation let's set up the recombination
% zones. We can use three methods: if we have a representative X-ray
% structure we can use Arnold's SCHEMA package (freely downloadable
% from http://www.che.caltech.edu/groups/fha/media/schema-tools.zip).
% Alternatively we can use an internal utility of MSAvolve, based on our 
% observation suggesting that recombination points tend to coincide with 
% positions of the msa with high relative entropy. At this point
% we still have the K>> prompt. The external msa has already been read in 
% and is stored as the variable 'REF_nmsa'. Before quitting, let's click on 
% the arrow symbol in the upper bar and then on the bottom right panel of 
% the four-panel figure (until is highlighted) and let's run:

[~,~,~,~,crossover_points] = get_nmsa_recomb_points(REF_nmsa, 0.999,0.01);

% it is a smoothing fit of the relative entropy plot on the right. The
% positions of the peaks can be retrieved by opening the variable:

openvar crossover_points

% see what happens if we change the parameters of the function:

[~,~,~,~,crossover_points] = get_nmsa_recomb_points(REF_nmsa, 0.9999,0.01);

% now you have more crossover points

[~,~,~,~,crossover_points] = get_nmsa_recomb_points(REF_nmsa, 0.9997,0.5);

% and now again fewer. You can learn about the meaning of the two
% parameters by opening the function:

open get_nmsa_recomb_points

% and its subfunction:

open peakdet

% The initial trial looked good (in general we want crossover points to be 
% separated by 10-15 residues). We can rerun it and copy the values of
% the variable 'crossover_points' into the corresponding variable of the 
% RUN_MSAvolve_ARSC.m script.
%
% crossover_points = [12 34 64 93 107 127];
%
% We also want to finalize our choice of covarions. We set the fraction of
% covarions to 15%:
%
% fcov =  15;
%
% We want the positions of the covarions to be based on the positions in
% the msa with medium levels of relative entropy (conservation), and also
% to include positions where there are some gaps in the alignment. We set:
%
% REL_ENTROPY = 1;
% cov_entropy_range = 'mid';
% cov_gaps = 1;
%
% If we are finally satisfied with our selections, we set:
%
% CHECK_BRANCHES = 0; 
% end cycle = 10;
%
% as we don't need to check the branches any more, and we start a trial run
% simulating 10 msa's. 
 
dbquit
close all
% save from the editor window
run RUN_MSAvolve_ARSC

% At the end of the run, the 10 msa's are store in the variable
% MSA_select_ALL in the internal Matlab numeric format.
% If we want to see msa # 7 we can type:

openvar MSA_select_ALL(:,:,7)

% We can convert the msa to a 'human' format:

MSA_no7 = int2aa(MSA_select_ALL(:,:,7));
openvar MSA_no7

%%
% We can now plot a few features of the trial run using some plotting
% functions provided in the PLOT directory of the toolbox: first we plot
% the distribution among all 10 MSAs of the similarity score among all
% sequences in each MSA: we compare this to the similarity score among the
% sequences in the Reference MSA.

Anc_Cons_Rel = figure;
set(Anc_Cons_Rel,'Units','normalized','Position',[0 0.2 0.5 0.8 ],...
    'Name','Ancestor Consensus Relationships'); clf;

subplot1 = subplot(2,1,1,'Parent',figure(gcf));
box(subplot1,'on');
hold(subplot1,'all');

histfit(MSA_sim_score,10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'FaceColor','b')
set(get(get(h(1),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude fit from legend

h = stem(REF_sim_score,max(hist(MSA_sim_score,10)+6),'fill','--');
%h = get(gca,'Children');
set(h,'LineWidth',1,'Color','b');
set(get(h,'BaseLine'),'LineStyle','none')
set(h,'MarkerFaceColor','red')
% set(get(get(h,'Annotation'),'LegendInformation'),...
%      'IconDisplayStyle','off'); % Exclude fit from legend
xlim([0.4 0.6]);
ylim([0 n+1]);
title('Simulated MSA versus Reference MSA Sim. Score ');

legend(subplot1,...
    'MSA Sim. Score','REF Sim. Score','Location','Best');
xlabel('Similarity Score ')
ylabel('Number of sequences ')
grid on

% We also plot the distribution of the mean correlation between the HMM
% emissions calculated from each MSA and the HMM emissions calculated from
% the Reference MSA, and the correlation between the relative entropy of
% each MSA and the relative entropy of the Reference MSA.

subplot2 = subplot(2,1,2,'Parent',figure(gcf));
box(subplot2,'on');
hold(subplot2,'all');

histfit(mean_emission_corr,10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'FaceColor',[1.0 0.5 0])
set(get(get(h(1),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude fit from legend

histfit(corr_REF_MSA_rel_entropy,10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'FaceColor','y')
set(get(get(h(1),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude fit from legend

title('REF-MSA Corr. Rel. Entropy/HMM Emissions');

legend(subplot2,'REF-MSA Corr. HMM Emiss.',...
    'REF-MSA Corr. Rel. Entr.','Location','Best');
xlabel('Correlation ')
ylabel('Number of sequences ')
grid on

%%
close all

%%
% Looking at the top plot we notice that the distribution of the simulated
% MSAs similarity score (MSA_sim_score) is centered at a value slightly
% smaller than the similarity score of the experimental msa (REF_sim_score).
saveas(gcf,'MSAvolve_MSA_statistics','fig');

%%
% Now, let's select the any msa (for example no. 6) by setting set n = 6,
% and execute the function:

n = 6;
run make_REF_MSA_comp_00
saveas(gcf,'REF_to_MSAvolve_MSA_comparison','fig');
close all

% The upper row of panels shows the properties of the experimental msa
% (REF_msa), while the lower row those of the selected simulated msa. In
% both you can rotate the spectral analysis box by clicking on the curved
% arrow button in the figure bar. The outermost panel on the bottom right
% shows the relative entropy of the REF_msa in blue and that of the
% simulated msa in red. We can close the figure and try another msa.

% We can also look at the average of all 10 runs. Execution of this
% command may take a long time if there are significantly more than just 10
% different MSAs.

run make_REF_MSA_comp_000
close all

% The covariance matrix histogram of these MSAs is reasonable. We are also
% happy about the fact the cluster analysis shows a nice distinction
% between the clusters selected originally.

% Altogether, these figures suggest the parameters we chose were very good,
% and we can be satisfied with our simulation: if we really wanted to
% improve further we could decrease a little the level of mutation and/or
% increase the level of recombination in order to increase the similarity
% among the sequences in each MSA to bring it even closer to the level in
% the Reference MSA.
%
% We notice at this point that it is also possible to fine tune the
% progression of the evolution by modifying the parameters that control the
% various evolution nodes (NODES 1,2,3). However these parameters have been
% optimized to yield good results with most proteins and in general there
% is no need to change them. In most cases a very good match with the
% properties of the experimental msa can be obtained by playing only with
% the cycles of mutations (mut_cycles), the cycles of recombination
% (rec_cycles), the number of crossover points (nzones, crossover_points),
% and the scaling of recombination (RECOMB_SCALING). Scaling recombination
% tends to amplify the relative effects of recombination with respect to
% mutations: if the level of recombination exceed that of mutations scaling
% will increase it even more. Conversely, if the level of recombination 
% is smaller than that of mutations scaling will decrease it further.
%
% We can also compare the evolutionary tree of a simulated MSA with the tree
% of the Reference MSA:

n = 6; % n = 8; n = 10; ... and so on.
make_REF_MSA_tree(REF_nmsa,MSA_select_ALL(:,:,n));    
saveas(gcf,'REF_to_MSAvolve_tree_comparison','fig');
close all

% We notice that some of branches of the tree are much longer in the
% simulated msa. If we wanted to eliminate this difference we could change
% the individual parameters in the 'NODES' section of RUN_MSAvolve_ARSC.
% Namely, we could increase the recombination rate inside the 3rd NODE
% (e.g. nhrt_range_5 = [2,8]), and decrease the mutation rate (e.g,
% mut_rate_range_8 = [0,3]). However, since this will bring the overall
% similarity between sequences up, we would also need to increase the
% mutation rate in the 1st and 2nd NODE of the evolution tree. If we are
% really interested in getting a perfect match between experimental and
% simulated msa, we could spend some time on this, but the stage at which
% we are is sufficient to study the coevolution process. Also, in a real
% study we would want to get a larger set of simulated msa's (maybe
% 100 msa's), but for our example the 10 msa's we already have serve us well.
%
%%
save 'ARSC_MSAvolve_run';

%%
close all
%% Internal coevolution matrices of synthetic MSAs produced by MSAvolve
% Here we look at the heat maps of the Reference MSA and the selected
% MSA generated by MSAvolve.

REF_SYNTHETIC_MSA = figure;
set(REF_SYNTHETIC_MSA,'Units','normalized','Position',[0 0.2 0.8 0.5 ],...
    'Name','Reference versus Synthetic MSA'); clf;

subplot1 = subplot(1,2,1,'Parent',figure(gcf));
imagesc(REF_nmsa);set(gca,'YDir','Normal');
box(subplot1,'on');
hold(subplot1,'all');

subplot2 = subplot(1,2,2,'Parent',figure(gcf));
imagesc(MSA_select_ALL(:,:,n));set(gca,'YDir','Normal');
box(subplot2,'on');
hold(subplot2,'all');

saveas(gcf,'Reference_versus_Synthetic_MSA','fig');

%%
% Here we plot various coevolution matrices generated by MSA volve. The
% meaning of these matrices can be found in the 2012 PLOS One paper:
% PLoS One. 2012;7(10):e47108. doi: 10.1371/journal.pone.0047108. 
% Epub 2012 Oct 16.

Coevolution_Total = figure; 
    set(Coevolution_Total,'Units','normalized','Position',[0.2 0.2 0.6 1 ],...
    'Name','Coevolution_Total'); clf;

        axes1 = axes('Parent',Coevolution_Total,...
                'Position',[0.02 0.52 0.46 0.46]);    
        imagesc(mut_COV_ALL(:,:,n));figure(gcf)
        axis equal tight
        set(gca,'Ydir','Normal')
        
        axes2 = axes('Parent',Coevolution_Total,...
                'Position',[0.52 0.52 0.46 0.46]);        
        imagesc(recomb_COV_ALL(:,:,n))
        set(gca,'Ydir','Normal')        
        axis equal tight
        
        axes3 = axes('Parent',Coevolution_Total,...
                'Position',[0.02 0.02 0.46 0.46]);        
        imagesc(cov_COV_ALL(:,:,n))
        set(gca,'Ydir','Normal')        
        axis equal tight

        axes4 = axes('Parent',Coevolution_Total,...
                'Position',[0.52 0.02 0.46 0.46]);
        plot_coev_noaxes(COV_ALL(:,:,n),...
            covar_vec);
        axis equal tight
        
        % The following is the same as COV_ALL
%         plot_coev_noaxes((recomb_COV_ALL(:,:,n)+mut_COV_ALL(:,:,n)+cov_COV_ALL(:,:,n)),...
%             covar_vec,MSA_rel_entropy(n,:),npos,hrf);
        
saveas(gcf,'MSAvolve_coevolution_maps','fig');

%%
close all

%% Coevolution matrices of synthetic MSAs from different methods
% Let's start by calculating the coevolution matrices using different
% methods (3D_MI, hpPCA, plmDCA_asym, GREMLIN). For the last three methods
% we choose a level of correction for gaps by setting the variable
% 'gap_corr_level': higher values (up to 3) provide stronger correction, 0 
% produces no correction. 

gap_corr_level = 3;

gW = zeros(npos,npos,end_cycle);
gapW0 = ones(npos,npos);
for n = 1:end_cycle
    gapW1 = correct_coevmat_forgaps(MSA_select_ALL(:,:,n));
    gapW2 = gapW1.^2;
    gapW3 = gapW1.^3;

    if exist('gap_corr_level','var')
        switch gap_corr_level
            case 0
        gW(:,:,n) = gapW0;
            case 1
        gW(:,:,n) = gapW1;
            case 2
        gW(:,:,n) = gapW2;
            case 3        
        gW(:,:,n) = gapW3;
        end
    else
        gW(:,:,n) = gapW0;        
    end
        
end

% set the similarity cutoff and gap correction levels for 3D_MI
sim_cutoff = 0.9;
gapcorr_1 = 0;
gapcorr_2 = 0;
gapcorr_3 = 3;
run add_md3_ZPX2_method

run add_hpPCA_method

rmpath(genpath('../COEVOLUTION_METHODS/FUNCTIONS/gplm_DCA_asymm'));
run add_plmDCA_asym_method

addpath(genpath('../COEVOLUTION_METHODS/FUNCTIONS/GREMLIN'));
run add_GREMLIN_method
rmpath(genpath('../COEVOLUTION_METHODS/FUNCTIONS/GREMLIN'));

%%
% We save the results:
 
save 'ARSC_MSAvolve_run';
 
%%
% We can plot a few coevolution matrices from different methods:

n = 5;

Coevolution_Methods = figure; 
    set(Coevolution_Methods,'Units','normalized','Position',[0.2 0.2 0.6 1 ],...
    'Name','Coevolution_Methods'); clf;

        axes1 = axes('Parent',Coevolution_Methods,...
                'Position',[0.02 0.52 0.46 0.46]);        
        plot_coev_noaxes(md3_ZPX2_ALL(:,:,n),covar_vec)
        clim_range = get(axes1,'CLim');
        coev_map_temp = squeeze(md3_ZPX2_ALL(:,:,n));
        clim_range(1) = nanmean(coev_map_temp(:));        
        set(gca,'Ydir','Normal')
        set(gca,'Ydir','Normal','CLim',clim_range)
        pause(1);
        axis equal tight
        
        axes2 = axes('Parent',Coevolution_Methods,...
                'Position',[0.52 0.52 0.46 0.46]);
        plot_coev_noaxes(hpPCA_ALL(:,:,n),...
            covar_vec);
        clim_range = get(axes2,'CLim');
        coev_map_temp = squeeze(hpPCA_ALL(:,:,n));
        clim_range(1) = nanmean(coev_map_temp(:));
        set(gca,'Ydir','Normal')
        set(gca,'Ydir','Normal','CLim',clim_range)
        axis equal tight
        
        axes3 = axes('Parent',Coevolution_Methods,...
                'Position',[0.02 0.02 0.46 0.46]);    
        plot_coev_noaxes(plmDCA_ALL(:,:,n),covar_vec)
        clim_range = get(axes3,'CLim');
        coev_map_temp = squeeze(plmDCA_ALL(:,:,n));
        clim_range(1) = nanmean(coev_map_temp(:));        
        set(gca,'Ydir','Normal')
        set(gca,'Ydir','Normal','CLim',clim_range)
        axis equal tight
        
        axes4 = axes('Parent',Coevolution_Methods,...
                'Position',[0.52 0.02 0.46 0.46]);        
        plot_coev_noaxes(GREMLIN_ALL(:,:,n),covar_vec)
        clim_range = get(axes4,'CLim');
        coev_map_temp = squeeze(GREMLIN_ALL(:,:,n));
        clim_range(1) = nanmean(coev_map_temp(:));        
        set(gca,'Ydir','Normal')
        set(gca,'Ydir','Normal','CLim',clim_range)        
        axis equal tight

% title('MSA no. 7: plmDCA coevolution matrix.');
% 
% xlabel('Column no. ')
% ylabel('Row no. ')
%  
saveas(gcf,'Methods_coevolution_maps_simulations','fig');

%%   
close all
  
%%
% After loading a simulation run of MSAvolve, we can  compares the
% performance of several methods on simulated data sets. In one type of
% plot we look at the distribution of true covarions among the top scores
% of the coevolution matrices produced by each method. We consider only a
% number of top scores equal to the total number of true covarying pairs
% generated by the program. The dashed orange line represents the
% fraction recorded by the program. In a different type of plot we use
% directly the COV or cov_COV matrix rather than the original vector of
% covarions 'covar_vec'. Results from all the runs are averaged. We will
% get a plot of the cumulative count of true covariation events recognized
% by the different methods, averaged over all 10 simulations. The dashed
% orange line represents the counts of true coevolution events recorded by
% the program.

% If it was not applied already in the scripts that launch each method, we
% can implement a correction of the coevolution matrices based on the
% presence of gaps at this point.

%%
totpairs = npos*npos;
cum_covarions = zeros(totpairs,end_cycle);
cum_plmDCA = zeros(totpairs,end_cycle);
cum_md3_ZPX2 = zeros(totpairs,end_cycle);
cum_hpPCA = zeros(totpairs,end_cycle);
cum_GREMLIN = zeros(totpairs,end_cycle);


for n = 1:end_cycle
% We can run the comparison on the total COV matrix ...
% covarions = sort_matrix_descend_2(COV_ALL(:,:,n),1);
% ... or just the true covarions matrix.
covarions = sort_matrix_descend_2(cov_COV_ALL(:,:,n),1);

cum_covarions(:,n) = cumsum(covarions(:,1));


% Run plmDCA
plmDCA = plmDCA_ALL(:,:,n);
plmDCA = plmDCA - min(plmDCA(:));
% plmDCA = plmDCA.*gW(:,:,n);
sorted_plmDCA = sort_matrix_descend_2(plmDCA,1);
sorted_plmDCA(:,4) = 0;
[~,ia,ib] = intersect(sorted_plmDCA(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_plmDCA(ia(i),4) = covarions(ib(i),1);
end
cum_plmDCA(:,n) = cumsum(sorted_plmDCA(:,4));


% Run md3_ZPX2
md3_ZPX2 = md3_ZPX2_ALL(:,:,n);
md3_ZPX2 = md3_ZPX2 - min(md3_ZPX2(:));
% md3_ZPX2 = md3_ZPX2.*gW(:,:,n);
sorted_md3_ZPX2 = sort_matrix_descend_2(md3_ZPX2,1);
sorted_md3_ZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_md3_ZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_md3_ZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_md3_ZPX2(:,n) = cumsum(sorted_md3_ZPX2(:,4));


% Run hpPCA
hpPCA = hpPCA_ALL(:,:,n);
hpPCA = hpPCA - min(hpPCA(:));
% hpPCA = hpPCA.*gW(:,:,n);
sorted_hpPCA = sort_matrix_descend_2(hpPCA,1);
sorted_hpPCA(:,4) = 0;
[~,ia,ib] = intersect(sorted_hpPCA(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_hpPCA(ia(i),4) = covarions(ib(i),1);
end
cum_hpPCA(:,n) = cumsum(sorted_hpPCA(:,4));


% Run GREMLIN
GREMLIN = GREMLIN_ALL(:,:,n);
GREMLIN = GREMLIN - min(GREMLIN(:));
% GREMLIN = GREMLIN.*gW(:,:,n);
sorted_GREMLIN = sort_matrix_descend_2(GREMLIN,1);
sorted_GREMLIN(:,4) = 0;
[~,ia,ib] = intersect(sorted_GREMLIN(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_GREMLIN(ia(i),4) = covarions(ib(i),1);
end
cum_GREMLIN(:,n) = cumsum(sorted_GREMLIN(:,4));

end

%%
% Here we get all the mean of all cumulative sums.
cum_covarions = nanmean(cum_covarions(:,1:n),2);
cum_plmDCA = nanmean(cum_plmDCA(:,1:n),2);
cum_md3_ZPX2 = nanmean(cum_md3_ZPX2(:,1:n),2);
cum_hpPCA = nanmean(cum_hpPCA(:,1:n),2);
cum_GREMLIN = nanmean(cum_GREMLIN(:,1:n),2);

%%
COEV_methods = figure; 
    set(COEV_methods,'Units','normalized','Position',[0 0.2 0.5 0.8 ],...
    'Name','COEV Methods'); clf;

subplot1 = subplot(2,1,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

histfit(fcov_cov_COV(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[1.0 0.5 0.0],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_COV(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_md3_ZPX2(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_hpPCA(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[1,0,0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_plmDCA(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_GREMLIN(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...save 'ARSC_MSAvolve_run1_with_methods';
    'IconDisplayStyle','off'); % Exclude histogram from legend

set(gca,'Xlim',[0.0 1],'Ylim',[0 n/5]);
legend('Ideal','Possible','3D\_MI','hpPCA','plmDCA','GREMLIN',...
       'Location','Best');
title('Coevolution Methods: distributions','FontSize',14,'FontWeight','n');
xlabel('Percentage of true covarions in the top 35 zscores');
ylabel('No. of occurrences in 10 trials');


subplot2 = subplot(2,1,2,'Parent',figure(gcf));
set(subplot2, 'xscale', 'log');
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');
       
	semilogx(cum_covarions,'--','Color',[1.0 0.5 0.0],'LineWidth',2);hold on
    semilogx(cum_md3_ZPX2,'-','Color','b','LineWidth',2);
    semilogx(cum_hpPCA,'-r','LineWidth',2);
    semilogx(cum_plmDCA,'-m','LineWidth',2);
    semilogx(cum_GREMLIN,'-','Color',[0,1,0],'LineWidth',2);
    plot(ncov,(-100:10E2:55E5),'ok','LineWidth',2,'MarkerSize',3);
	hold off
    set(gca,'Xlim',[0,1E2],'Ylim',[0,max(cum_covarions(1:1.1E2))]);

% 	legend('Ideal','3D\_MI','hpPCA','plmDCA','Location','NorthWest');

    title('Covariation Methods: cumulative sums','FontSize',14,'FontWeight','n');
	xlabel('Number of top pairs identified by a method');
	ylabel('Cumulative count of true covariation events');
    box('on');
    grid('on');

saveas(gcf,'Methods_coevolution_detection','fig');
    
%%
close all
save 'ARSC_MSAvolve_run';

%%

%-----------------Analysis of the experimental MSA-------------------------
% We can use some of the utilities of MSAvolve to analyze the experimental
% MSA and compare the coevolution matrix derived by different methods with
% the distance matrix of the corresponding reference structure of the ArsC
% family. 
%--------------------------------------------------------------------------


%% Choise of msa for covariation analysis
% For all subsequent analyses we can use the nmsa that we have derived in 
% the previous section. 
nmsa = REF_nmsa;

% Alternatively we can load an alignment prepared with different programs.  

% msafile = 'ARSC_95.faln';
% msafile_type = 'faln';
% 
% switch msafile_type 
%     case 'faln'
%     [REF_smsa,REF_nmsa] = faln_to_nmsa(msafile);
%     case 'aln'        
%     [REF_smsa,REF_nmsa] = aln_to_nmsa(msafile);
% end
% 
% smsa = REF_smsa;
% nmsa = REF_nmsa;

%%
pdbfile = 'ARSC.pdb';
START = 3;
END = 140;
PDB_START = 3;
PDB_END = 140;
nmsa = nmsa(:,PDB_START:PDB_END);
REF_length = numel(nmsa(1,:));
cmsa = int2aa(nmsa);
[nrows,ncols] = size(nmsa);

%% NMSA contraction to decrease noise 
% Find the columns with less than 25% of gaps. Generally this is not
% necessary, but it may help with MSAs with a lot of gaps.
% gap_fraction_nmsa_vec = sum(nmsa == 25)/nrows;
% lowgap_ind = find(gap_fraction_nmsa_vec<0.25);
% small_nmsa = nmsa(:,lowgap_ind);
% [snrows,sncols] = size(small_nmsa);

%% Covariation_distance maps
% Here we set some file names to store figures.
param_comparison_file = 'COMPARE_CONTACT_PREDICTIONS';
contact_map_file = 'COMPARE_CONTACT_MAP';

%% Covariation analysis

[~,md3_ZPX2] = ...
    NMSA_to_mdMI(nmsa,'GAPS','3D','FULL',0.9,1,22,0,0,3,12);

% Our implementation of the PSICOV algorithm: extremely fast and working
% also with less than 500 sequences.
% [slPSICOV_ZPX2] = NMSA_to_slPSICOV(nmsa,'NOGAPS',0.9,1.0,...
%     20,'SUM','NONE','fro',0.0001,0.45,0.015,100,0.0,0.0001,'RHO',0,1);
% For the original PSICOV by Jones, which requires a minimum of 500 
% sequences in the MSA use the following lines:
% [PSICOV] = get_nmsa_covar_vec(nmsa,30,'PSICOV');
% PSICOV already has a MIP correction, so we only add a ZPX2 correction.
% PSICOV_ZPX2 = MIP_to_ZPX2(PSICOV);

[plmDCA] = get_nmsa_covar_vec(nmsa,30,'plmDCA_asym');
% plmDCA already has a MIP correction, so we only add a ZPX2 correction.
plmDCA_ZPX2 = MIP_to_ZPX2(plmDCA);

[hpPCA] = get_nmsa_covar_vec(nmsa,30,'hpPCA');
% plmDCA already has a MIP correction, so we only add a ZPX2 correction.
hpPCA_ZPX2 = MIP_to_ZPX2(hpPCA);

addpath(genpath('../COEVOLUTION_METHODS/FUNCTIONS/GREMLIN'));
[GREMLIN] = get_nmsa_covar_vec(nmsa,30,'GREMLIN');
% GREMLIN already has a MIP correction, so we only add a ZPX2 correction.
GREMLIN_ZPX2 = MIP_to_ZPX2(GREMLIN);
! rm gremlin.matrix
% GREMLIN redefines matlab int2aa and aa2int; therefore we remove it from
% the path.
rmpath(genpath('../COEVOLUTION_METHODS/FUNCTIONS/GREMLIN'));
! rm MSA_select.selex
%% Save everything
save 'ARSC_MSAvolve_run';

%% Gap correction 
% Here we calculate a matrix of weights to correct for the presence of
% gaps. 

gapW0 = ones(ncols,ncols);
gapW1 = correct_coevmat_forgaps(nmsa);
gapW2 = gapW1.^2;
gapW3 = gapW1.^3;

% Here we choose how to corrects the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW1,gapW2,gapW3) can be applied. Setting gW = gapW0 means
% no correction.
gW = gapW3;

% We only correct DCA and GREMLIN type maps as the other methods are already 
% corrected for gaps.

% PSICOV_ZPX2_1 = PSICOV_ZPX2;
% PSICOV_ZPX2_2 = PSICOV_ZPX2_1 - min(PSICOV_ZPX2_1(:));
% PSICOV_ZPX2 = PSICOV_ZPX2_2.*gW;

hpPCA_ZPX2_1 = hpPCA_ZPX2;
hpPCA_ZPX2_2 = hpPCA_ZPX2_1 - min(hpPCA_ZPX2_1(:));
hpPCA_ZPX2 = hpPCA_ZPX2_2.*gW;

plmDCA_ZPX2_1 = plmDCA_ZPX2;
plmDCA_ZPX2_2 = plmDCA_ZPX2_1 - min(plmDCA_ZPX2_1(:));
plmDCA_ZPX2 = plmDCA_ZPX2_2.*gW;

GREMLIN_ZPX2_1 = GREMLIN_ZPX2;
GREMLIN_ZPX2_2 = GREMLIN_ZPX2_1 - min(GREMLIN_ZPX2_1(:));
GREMLIN_ZPX2 = GREMLIN_ZPX2_2.*gW;

%% Display individual coevolution matrices.

Coevolution_Methods_Experimental = figure; 
    set(Coevolution_Methods_Experimental,'Units','normalized','Position',[0.2 0. 0.6 1 ],...
    'Name','Coevolution Methods - Experimental MSA'); clf;

        axes1 = axes('Parent',Coevolution_Methods_Experimental,...
                'Position',[0.02 0.52 0.46 0.46]);    
        imagesc(md3_ZPX2);
        clim_range = get(axes1,'CLim');
        clim_range(1) = nanmean(md3_ZPX2(:));        
        set(gca,'Ydir','Normal');
        % set(gca,'CLim',clim_range);
        pause(1); % necessary to resolve a graphic bug in Matlab
        axis equal tight
        
        axes2 = axes('Parent',Coevolution_Methods_Experimental,...
                'Position',[0.52 0.52 0.46 0.46]);        
        imagesc(hpPCA_ZPX2);
        clim_range = get(axes2,'CLim');
        clim_range(1) = nanmean(hpPCA_ZPX2(:));        
        set(gca,'Ydir','Normal');
        % set(gca,'CLim',clim_range);
        axis equal tight
                
        axes3 = axes('Parent',Coevolution_Methods_Experimental,...
                'Position',[0.02 0.02 0.46 0.46]);        
        imagesc(plmDCA_ZPX2);
        clim_range = get(axes3,'CLim');
        clim_range(1) = nanmean(plmDCA_ZPX2(:));        
        set(gca,'Ydir','Normal');
        % set(gca,'CLim',clim_range);
        axis equal tight

        axes4 = axes('Parent',Coevolution_Methods_Experimental,...
                'Position',[0.52 0.02 0.46 0.46]);
        imagesc(GREMLIN_ZPX2);
        clim_range = get(axes4,'CLim');
        clim_range(1) = nanmean(GREMLIN_ZPX2(:));
        set(gca,'Ydir','Normal');
        % set(gca,'CLim',clim_range);
        axis equal tight
        
% title('MSA coevolution matrices.');
% 
% xlabel('Column no. ')
% ylabel('Row no. ')
%  
saveas(gcf,'Methods_coevolution_maps_experimental','fig');

%%
close all

%% Save everything
save 'ARSC_MSAvolve_run';

%% Superposition of covariation and distance matrices 

%---------------------1st set----------------------------------------------
% Usage: [c_distances,sorted_c_distances,sorted_mat1,...
%           sorted_mat2,sorted_mat3,sorted_mat4,pdbstruct] = ...
%           coev_distance_matrix_3(pdbfile,chain,first_res_no,last_res_no,...
%           mat1,mat2,mat3,mat4,radius,near,npreds,plot_dist)
% chain: 1,2,3,...
% first_res_no: number of the first residue with coordinates in the pdb file
% last-res_no: number of the last residue with coordinates in the pdb file
% mat1, mat2 ...: coevolution matrices
% radius: threshold distance to select entries in the protein distance matrix
% near: minimum separation between residues in sequence to be included in
% the analysis: 1 = consecutive; 2 = separated by 1 intervening residue; 
% 3 = separated by 2 intervening residues; and so on ...
% npreds: this is 3 times the number of top positions that will be plotted
% on top of a distance matrix.
% plot_dist can have different values:
% 0 : no plot
% 1 : sparsity plot of the distance matrix (yellow on white background).
% 2 : heat map of the distance matrix with everything beyond the radius zeroed.
% and with blue background and colors representing the number of atomic
% interactions
% else : like 2 (but cyan on white background), with the covarions identified 
% by different methods overlaid on top.

%--------------------------------------------------------------------------
near = 1;
near1 = near;
ncov = REF_length;
radius = 8;

[c_distances,sorted_c_distances,sorted_md3_ZPX2,sorted_hpPCA_ZPX2,...
    sorted_plmDCA_ZPX2,sorted_GREMLIN_ZPX2] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,md3_ZPX2,...
    hpPCA_ZPX2,plmDCA_ZPX2,GREMLIN_ZPX2,radius,near,ncov,3,...
    [0,0,1],[1 0 0],[1 0 1],[0 1 0],...
    [0.9 0.9 0.9],[0.9 0.9 0.9],[0.9 0.9 0.9],[0.9 0.9 0.9],...
    [1.0 1.0 1.0],[1.0 1.0 1.0],[1.0 1.0 1.0],[1.0 1.0 1.0]...
    );

saveas(gcf,contact_map_file,'fig');

%---------------------2nd set----------------------------------------------
near = 21;
near2 = near;
ncov = REF_length;
radius = 8;

[c_distances,sorted_c_distances,sorted_md3_ZPX2_B,sorted_hpPCA_ZPX2_B,...
    sorted_plmDCA_ZPX2_B,sorted_GREMLIN_ZPX2_B] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,md3_ZPX2,...
    hpPCA_ZPX2,plmDCA_ZPX2,GREMLIN_ZPX2,radius,near,ncov,0,...
    [0,0,1],[1 0 0],[1 0 1],[0 1 0],...
    [0.9 0.9 0.9],[0.9 0.9 0.9],[0.9 0.9 0.9],[0.9 0.9 0.9],...
    [1.0 1.0 1.0],[1.0 1.0 1.0],[1.0 1.0 1.0],[1.0 1.0 1.0]...
    );

%%
close all

%% Save everything
save 'ARSC_MSAvolve_run';

%% Validation by prediction of close contacts in the reference X-ray structure

COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.8 0.4 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;

subplot1 = subplot(1,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

ncov = REF_length;
X = (1:ncov);
semilogy(X,sorted_md3_ZPX2(X,6),'-b','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2(X,6),'-r','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2(X,6),'-m','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2(X,6),'-g','LineWidth',2);

hold off
separation1 = num2str(near1 - 1);
distance = num2str(radius);
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of all the residues in the sequence ';
string3 = ' ';
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,0.12]);
    legend('3D\_MI','hpPCA','plmDCA','GREMLIN',...
           'Location','NorthWest');
% 	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
xlabel('Number of top pairs identified by a method');
ylabel({string1;string2;string3});


subplot2 = subplot(1,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_B(X,6),'-b','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_B(X,6),'-r','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_B(X,6),'-m','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_B(X,6),'-g','LineWidth',2);

hold off
separation2 = num2str(near2 - 1);
distance = num2str(radius);
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation2 ' intervening positions in the sequence '];
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,0.16]);
xlabel('Number of top pairs identified by a method');
ylabel({string1;string2;string3});


%% 
saveas(gcf,param_comparison_file,'fig');

%%
close all

%% Save everything
save 'ARSC_MSAvolve_run';

%% Visualization of coevolution in the reference structure

%--------------------------------------------------------------------------
% Generation of a Pymol script to display the top 50 coevolving pairs in  
% the reference structure of the protein family.
%--------------------------------------------------------------------------

% To apply to a different family change the 'family_root' name. To generate
% the script with a different method change the 'method_root' name
family_root = 'ARSC_';
method_root = 'plmDCA';
best_map = eval([method_root '_ZPX2']);
best_sorted = eval(['sorted_' method_root '_ZPX2']);

imagesc(best_map);set(gca,'YDir','Normal')
%%
% Information about the pdb file:

pdbfile = [family_root 'simple.pdb'];
ARSC_pdb = pdbread(pdbfile);
START = 1;
END = 138;
PDB_START = 3;
PDB_END = 140;
nmsa = nmsa(:,START:END);
REF_length = numel(nmsa(1,:));
cmsa = int2aa(nmsa);
[nrows,ncols] = size(nmsa);

% Here we determine the maximum pair contribution of each residue:
max_res_score = zeros(ncols,1);
for i = 1:ncols
    max_res_score(i) = max(best_map(i,:));
end

% Here we find out the pairs that are less than a given threshold in
% angstroms apart
dist_threshold = 13;
close_pairs_ind = best_sorted(1:50,4)<=dist_threshold;
far_pairs_ind = best_sorted(1:50,4)>=dist_threshold;
close_pairs = best_sorted(close_pairs_ind,:);
far_pairs = best_sorted(far_pairs_ind,:);
    
mod_pdb = change_pdb_beta(ARSC_pdb,[(PDB_START:PDB_END)' max_res_score],1);
pdbwrite([family_root method_root '.pdb'],mod_pdb);

%%
fileID = fopen([family_root method_root '.pml'],'w');
pml_header_1 = 'cmd.delete("all")';
pml_header_2 = ['load ' family_root method_root '.pdb'];
pml_header_3 = ['cmd.spectrum("b",selection=("' family_root method_root '"),quiet=0)'];
close_pairs_length = size(close_pairs,1);
far_pairs_length = size(far_pairs,1);
pml_text_1 = cell(close_pairs_length*3,1);
pml_text_2 = cell(far_pairs_length,1);

j = 0;
for i = 1:close_pairs_length
    j = j+1;
    pml_text_1{j} = ['dist ////' int2str(close_pairs(i,2)) '/CA/,////' ...
                 int2str(close_pairs(i,3)) '/CA/'];
    j = j+1;        
    pml_text_1{j} = ['select sele, resi ' int2str(close_pairs(i,2)) '+' ...
                        int2str(close_pairs(i,3))];
    j = j+1;                    
    pml_text_1{j} = 'cmd.show("sticks","sele")';             
end    

j = 0;
for i = 1:far_pairs_length
    j = j+1;
    pml_text_2{j} = ['bond ////' int2str(far_pairs(i,2)) '/CA/,////' ...
                 int2str(far_pairs(i,3)) '/CA/'];
end

% Use the following lines to highlight some specific residues.
pml_text_3 = [{'select sele, resi 13'};{'cmd.show("dots","sele")'};...
              {'cmd.show("sticks","sele")'};
              {'set dot_density, 2'};{'set dot_radius, 0'};...
              {'set dot_width, 1'}];

pml_tail_1 = {['cmd.hide("lines","ARSC_' method_root '")']};

pml_tail_2 = cell(far_pairs_length*2,1);
j = 0;
for i = 1:far_pairs_length
    j = j+1;         
    pml_tail_2{j,1} = ['set_bond line_color,white, ////' int2str(far_pairs(i,2)) '/CA/,////' ...
                 int2str(far_pairs(i,3)) '/CA/']; 
    j = j+1;         
    pml_tail_2{j,1} = ['show lines, ////' int2str(far_pairs(i,2)) '/CA/ | ////' ...
                 int2str(far_pairs(i,3)) '/CA/'];                         
end    

pml_tail_3 = {['cmd.show("ribbon","' family_root method_root '")']};
pml_tail_4 = {['cmd.hide("((byres (' family_root method_root '))&(n. c,o,h|(n. n&!r. pro)))")']};
pml_tail_5 = {'cmd.color(13,"dist*")'};
pml_tail_6 = {'cmd.hide("labels","all")'};
pml_tail_7 = {'cmd.bg_color(''gray'')'};
pml_tail_8 = {'cmd.show("labels","all")'};

pml_view_1 = {'set ortho,0'};
pml_view_2 = [{'view 1, store'};
{'set depth_cue = 0'};
{'set ray_trace_fog, 0'};
{'set ray_shadows, 0'}
{'ray 1000,1200'};
{['png ' family_root method_root '.png, width=1000, height=1200, dpi=300, ray=1']}];

pml_tail = [pml_tail_1;pml_tail_2;pml_tail_3;pml_tail_4;pml_tail_5;...
            pml_tail_6;pml_tail_7;pml_view_1;pml_view_2];

MERGED = [pml_header_1;pml_header_2;pml_header_3;pml_text_1;...
          pml_text_2;pml_text_3;pml_tail];
merged_length = length(MERGED);

formatSpec = '\r\n%s%\r\n';

for i = 1:merged_length
fprintf(fileID,formatSpec,MERGED{i});
end

fclose(fileID);

%%
close all

%% Save everything
save 'ARSC_MSAvolve_run';

%%
%--------------------------------------------------------------------------
% IMPORTANT!!!
% Each time you open one of the original functions or scripts of the
% toolbox to study it, immediately save it in your working directory with a
% local name. If you modify one of the native functions and save the
% changes, the Toolbox may not work properly any longer if you made a
% mistake.
%--------------------------------------------------------------------------













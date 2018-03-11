% This script contains examples of how to plot various results from the
% MSAvolve run
%% 
% First we analyze the phylogenetic trees.

[REF_temp_tree] = make_tree(REF_nmsa);    
[MSA_temp_tree] = make_tree(MSA_select_ALL(:,:,n));

% Here we compare the phylogenetic trees of the REF and evolved MSA using
% a radial plot that highlights the clusters.
% The figure is made by first creating two temporary figures and then 
% transfering the axes to a single one. This is very ackward; if you know
% how to do it in a more elegant way, please tell me.

[~,MSA_temp_tree] = make_REF_MSA_tree(REF_nmsa,MSA_select_ALL(:,:,n));    
     
% Here we compare the REF and MSA covariance matrices and relative 
% entropies.
    
% Here we cluster the MSA based on the spectral analysis 
% of its covariance or distance matrix. First we find the relative entropy.

[MSA_Clusters] = make_REF_MSA_comp_3(REF_nmsa,MSA_select_ALL(:,:,n),...
        DISTANCE,DIMENSIONS,MAX_CLUSTERS,...
        REF_rel_entropy,MSA_rel_entropy(n,:),covar_vec);
        
% Here we cluster the simulated MSA based on the spectral analysis of its
% covariance or distance matrix.

[MSA_Clusters] = make_REF_MSA_comp_3(REF_nmsa,MSA_select_ALL(:,:,n),...
        DISTANCE,DIMENSIONS,MAX_CLUSTERS,...
        REF_rel_entropy,MSA_rel_entropy(n,:),covar_vec);

% Here we cluster the simulated MSA based on the consecutive index of the
% sequences in the branches.
    
[MSA_Clusters] = make_REF_MSA_comp_4(REF_nmsa,MSA_select_ALL(:,:,n),...
        DISTANCE,DIMENSIONS,MAX_CLUSTERS,...
        REF_rel_entropy,MSA_rel_entropy(n,:),b_ind2,covar_vec);

% Same comparisons carried out on the means

[~] = make_REF_MSA_comp_2(REF_nmsa,REF_cov_mat,mean_MSA_cov_mat,...
        DIMENSIONS,MAX_CLUSTERS,...
        REF_rel_entropy,mean_MSA_rel_entropy,covar_vec);

%%
% Here we plot some relationships between the ancestor and the consensus
% sequence of the evolved msa. Try also 'gev' as a distribution for
% histfit.

Anc_Cons_Rel = figure;
set(Anc_Cons_Rel,'Units','normalized','Position',[0 0.5 0.5 0.5 ],...
    'Name','Ancestor Consensus Relationships'); clf;

subplot1 = subplot(3,1,1,'Parent',figure(gcf));
box(subplot1,'on');
hold(subplot1,'all');

histfit(MSA_cons_anc_sim_score,10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[1 0 0]);
set(h(2),'FaceColor',[.8 .8 1])
set(get(get(h(1),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude fit from legend

histfit(MSA_sim_score,10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'FaceColor','b')
set(get(get(h(1),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude fit from legend

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

h = stem(REF_sim_score,max(hist(MSA_sim_score,10)+6),'fill','--');
%h = get(gca,'Children');
set(h,'LineWidth',1,'Color','b');
set(get(h,'BaseLine'),'LineStyle','none')
set(h,'MarkerFaceColor','red')
% set(get(get(h,'Annotation'),'LegendInformation'),...
%      'IconDisplayStyle','off'); % Exclude fit from legend

title('MSA Sim. Score - Ancestor-Consensus Sim. Score - REF-MSA Corr. Rel. Entropy/HMM Emissions');


subplot2 = subplot(3,1,2,'Parent',figure(gcf),'YDir','reverse',...
    'Layer','top');
xlim(subplot2,[0.5 (npos + .5)]);
ylim(subplot2,[0.5 (end_cycle + .5)]);
box(subplot2,'on');
hold(subplot2,'all');
imagesc(con_anc_int_ind,'Parent',subplot2);
title('Ancestor-Consensus Matching Positions');

subplot3 = subplot(3,1,3,'Parent',figure(gcf),'YDir','reverse',...
    'Layer','top');
xlim(subplot3,[0.5 (npos + .5)]);
ylim(subplot3,[0.5 (end_cycle + .5)]);
box(subplot3,'on');
hold(subplot3,'all');
imagesc(MSA_rel_entropy,'Parent',subplot3);
title('MSA Relative Entropy');

legend(subplot1,'Anc-Cons Sim. Score',...
    'MSA Sim. Score','REF-MSA Corr. HMM Emiss.',...
    'REF-MSA Corr. Rel. Entr.','REF Sim. Score','Location','best');

%%
% Here we plot various coevolution matrices.

        plot_tot_coev(COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);
        
        plot_tot_coev(glob_COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
    
        plot_tot_coev(recomb_COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
    
        plot_tot_coev(mut_COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
    
        plot_tot_coev(cov_COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);
                
%--------------------------------------------------------------------------        
% The following is very useful if we are interested in the effects of 
% recombination on the covarions and cross-covarions count.
        covrecomb_COV_ALL = zeros(npos,npos,end_cycle);
        truerecomb_COV_ALL = zeros(npos,npos,end_cycle);
        truecov_COV_ALL = zeros(npos,npos,end_cycle);
        for n = 1:end_cycle 
            cov_COV_ind = cov_COV_ALL(:,:,n) ~= 0;
            temp_recomb_COV = recomb_COV_ALL(:,:,n);
            temp_recomb_COV(~cov_COV_ind) = 0;
            covrecomb_COV_ALL(:,:,n) = temp_recomb_COV;
            truerecomb_COV_ALL(:,:,n) = ...
            recomb_COV_ALL(:,:,n) - covrecomb_COV_ALL(:,:,n);
            truecov_COV_ALL(:,:,n) = ...
            cov_COV_ALL(:,:,n) + covrecomb_COV_ALL(:,:,n);
            clear temp_recomb_COV
        end        
%         plot_tot_coev(truerecomb_COV_ALL(:,:,n),...
%             covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(truecov_COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);       
        plot_tot_coev(covrecomb_COV_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
%--------------------------------------------------------------------------        
    
        plot_tot_coev(MI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(bayesMI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(SU_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(NMI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(MIP_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(ZPX_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(ZPX2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf); 
        plot_tot_coev(ZRES2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(ZMI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(nZNMI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(sZNMI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(OMES_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(ELSC_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(McBASC_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(fodorMI_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(fodorSCA_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);            
        plot_tot_coev(SSEM_SCA_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(RSEM_SCA_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(RAMA_SCA_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(RAMA_eigSCA_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(DCA_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(nbZPX2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(dbZPX2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(gbZPX2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(fgbZPX2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
        plot_tot_coev(dgbZPX2_ALL(:,:,n),...
            covar_vec,MSA_rel_entropy(n,:),npos,hrf);        
    
%%
% Here we plot all the types of coevolution matrices

    Coevolution_All = figure; 
    set(Coevolution_All,'Units','normalized','Position',[0 0.0 1 0.7 ],...
    'Name','Coevolution_ALL'); clf;

    axes1 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.0 0.55 0.2 0.4]);
        plot_coev_noaxes(COV_ALL,covar_vec);
    title('COV coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));

    axes2 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.2 0.55 0.2 0.4]);
        plot_coev_noaxes(cov_COV_ALL,covar_vec);
    title('cov COV coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
        
    axes3 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.4 0.55 0.2 0.4]);
        plot_coev_noaxes(MI_ALL,covar_vec);
    title('MI coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes4 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.6 0.55 0.2 0.4]);
        plot_coev_noaxes(ZMI_ALL,covar_vec);        
    title('ZMI coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes5 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.8 0.55 0.2 0.4]);
        plot_coev_noaxes(MIP_ALL,covar_vec);
    title('MIP coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes6 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.0 0.05 0.2 0.4]);
        plot_coev_noaxes(ZPX2_ALL,covar_vec);
    title('ZPX2 coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes7 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.2 0.05 0.2 0.4]);
        plot_coev_noaxes(OMES_ALL,covar_vec);
    title('OMES coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes8 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.4 0.05 0.2 0.4]);
        plot_coev_noaxes(ELSC_ALL,covar_vec);
    title('ELSC coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes9 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.6 0.05 0.2 0.4]);
        plot_coev_noaxes(McBASC_ALL,covar_vec);
    title('McBASC coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    axes10 = axes('Parent',Coevolution_All,'YTick',zeros(1,0),...
    'Position',[0.8 0.05 0.2 0.4]);
        plot_coev_noaxes(RAMA_SCA_ALL,covar_vec);
    title('RAMA SCA coevolution matrix','FontSize',14,'FontWeight','n');
    set(gca,'YTick',zeros(1,0));
    
    
%%
% We can try to plot some of the statistics about the capacity of 
% each method to identify coevolving positions. However the type of 
% probability distribution to use cannot be generalized. Some examples are 
% given below.

COEV_methods_1 = figure; 
    set(COEV_methods_1,'Units','normalized','Position',[0 0.2 0.8 0.6 ],...
    'Name','COEV Methods: sensitivity'); clf;

subplot1 = subplot(2,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

histfit(fcov_COV(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_cov_COV(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[1.0 0.5 0.0],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_MI(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_MIP(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','g');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_bayesMI(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.3,0.6,0.2]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_ZPX(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','y');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_ZPX2(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_ZRES2(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_ZMI(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','m');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_NMI(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','c');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_nZNMI(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','k');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_DCA(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_nbZPX2(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_dbZPX2(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_gbZPX2(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','y');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_fgbZPX2(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','k');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_dgbZPX2(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend

set(gca,'Xlim',[0,0.7],'Ylim',[0,40]);
%legend('totCOV','covCOV','MI','bayesMI','ZPX','ZPX2','ZRES','ZMI','MIP','NMI','ZNMI','DCA');
legend('totCOV','covCOV','MI','logR','ZPX2','DCA','nbZPX2',...
    'dbZPX2','dgbZPX2');
title('MI Methods','FontSize',14,'FontWeight','n');
xlabel('Percentage of true covarions in the top 21 zscores');
ylabel('No. of occurrences in 100 trials');

%--------------------------------------------------------------------------
subplot2 = subplot(2,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

histfit(fcov_COV(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.5 0.5 0.5],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_cov_COV(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[1.0 0.5 0.0],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(fcov_fodorMI(:,2),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','r');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_MI(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_OMES(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_McBASC(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_ELSC(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_fodorSCA(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(fcov_RSEM_SCA(:,2),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% scatter(0.0005,99.0,60,'y','filled')
set(gca,'Xlim',[0,0.3],'Ylim',[0,40]);
legend('totCOV','covCOV','MI','OMES','McBASC','ELSC','fodorSCA','ramaSCA');
title('Non-MI Methods','FontSize',14,'FontWeight','n');
xlabel('Percentage of true covarions in the top 21 zscores');
ylabel('No. of occurrences in 100 trials');

%--------------------------------------------------------------------------
subplot3 = subplot(2,2,3,'Parent',figure(gcf));
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'all');

histfit(mean_cov_zscore_COV(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[ .5 .5 .5],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_cov_COV(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[ 1 .5 0],'LineStyle','--');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_MI(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_MIP(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','g');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_bayesMI(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.3,0.6,0.2]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_ZPX(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','y');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_ZPX2(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_ZRES2(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_ZMI(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','m');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_NMI(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','c');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_nZNMI(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','k');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_DCA(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_nbZPX2(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_dbZPX2(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_gbZPX2(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','y');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
% histfit(mean_cov_zscore_fgbZPX2(:,1),10,'logistic');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','k');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_dgbZPX2(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
hold off

set(gca,'Xlim',[0 20]);
legend('totCOV','covCOV','MI','logR','ZPX2','DCA','nbZPX2',...
    'dbZPX2','dgbZPX2');
title('MI Methods','FontSize',14,'FontWeight','n');
xlabel('Mean zscore of all true covarying pairs');
ylabel('No.of occurrences in 100 trials');

%--------------------------------------------------------------------------
subplot4 = subplot(2,2,4,'Parent',figure(gcf));
box(subplot4,'on');
grid(subplot4,'on');
hold(subplot4,'all');

% histfit(mean_cov_zscore_fodorMI(:,1),10,'normal');
% h = get(gca,'Children');
% set(h(1),'LineWidth',2,'Color','r');
% set(h(2),'LineStyle','none','FaceColor','none');
% set(get(get(h(2),'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_MI(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_OMES(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_McBASC(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_ELSC(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_fodorSCA(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(mean_cov_zscore_RSEM_SCA(:,1),10,'logistic');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
hold off

% set(gca,'Xlim',[0 20]);
set(gca,'Ylim',[0 30]);
legend('MI','OMES','McBASC','ELSC','fodorSCA','ramaSCA','Location','best');
title('Non-MI Methods','FontSize',14,'FontWeight','n');
xlabel('Mean zscore of all true covarying pairs');
ylabel('No.of occurrences in 100 trials');

%--------------------------------------------------------------------------
%%
COEV_methods_2 = figure; 
    set(COEV_methods_2,'Units','normalized','Position',[0 0.0 0.8 1.0 ],...
    'Name','COEV Methods: accuracy'); clf;

subplot1 = subplot(3,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

histfit(corr_cov_zscore_MI(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZPX(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZPX2(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZRES2(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZMI(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_MIP(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_NMI(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_nZNMI(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','k');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
hold off
legend('MI','ZPX','ZPX2','ZRES','ZMI','MIP','NMI','ZNMI');
title('MI Methods','FontSize',14,'FontWeight','n');
xlabel('Pearson corr. of method zscores to totCOV zscores for true covarions');
ylabel('No.of occurrences in 100 trials');
set(gca,'Ylim',[0,30]);

subplot2 = subplot(3,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

histfit(corr_cov_zscore_fodorMI(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_OMES(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_McBASC(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ELSC(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_fodorSCA(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_RAMA_SCA(:,1),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend

legend('fodorMI','OMES','McBASC','ELSC','fodorSCA','ramaSCA');
title('Non-MI Methods','FontSize',14,'FontWeight','n');
xlabel('Pearson corr. of method zscores to totCOV zscores for true covarions');
ylabel('No.of occurrences in 100 trials');
set(gca,'Ylim',[0,30]);

subplot3 = subplot(3,2,3,'Parent',figure(gcf));
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'all');

histfit(corr_cov_zscore_MI(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZPX(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZPX2(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZRES2(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZMI(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_MIP(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_NMI(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_nZNMI(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','k');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
hold off
legend('MI','ZPX','ZPX2','ZRES','ZMI','MIP','NMI','ZNMI');
title('MI Methods','FontSize',14,'FontWeight','n');
xlabel('Spearman corr. of method zscores to totCOV zscores for true covarions');
ylabel('No.of occurrences in 100 trials');

subplot4 = subplot(3,2,4,'Parent',figure(gcf));
box(subplot4,'on');
grid(subplot4,'on');
hold(subplot4,'all');

histfit(corr_cov_zscore_fodorMI(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_OMES(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_McBASC(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ELSC(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_fodorSCA(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_RAMA_SCA(:,2),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend

legend('fodorMI','OMES','McBASC','ELSC','fodorSCA','ramaSCA');
title('Non-MI Methods','FontSize',14,'FontWeight','n');
xlabel('Spearman corr. of method zscores to totCOV zscores for true covarions');
ylabel('No.of occurrences in 100 trials');
set(gca,'Ylim',[0,30]);

subplot5 = subplot(3,2,5,'Parent',figure(gcf));
box(subplot5,'on');
grid(subplot5,'on');
hold(subplot5,'all');

histfit(corr_cov_zscore_MI(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZPX(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZPX2(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZRES2(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color',[0.0 0.5 1.0]);
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ZMI(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_MIP(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_NMI(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_nZNMI(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','k');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
hold off
legend('MI','ZPX','ZPX2','ZRES','ZMI','MIP','NMI','ZNMI');
title('MI Methods','FontSize',14,'FontWeight','n');
xlabel('Kendall corr. of method zscores to totCOV zscores for true covarions');
ylabel('No.of occurrences in 100 trials');

subplot6 = subplot(3,2,6,'Parent',figure(gcf));
box(subplot6,'on');
grid(subplot6,'on');
hold(subplot6,'all');

histfit(corr_cov_zscore_fodorMI(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','r');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_OMES(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','g');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_McBASC(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','m');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_ELSC(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','b');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_fodorSCA(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','c','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend
histfit(corr_cov_zscore_RAMA_SCA(:,3),10,'normal');
h = get(gca,'Children');
set(h(1),'LineWidth',2,'Color','y','LineStyle','-');
set(h(2),'LineStyle','none','FaceColor','none');
set(get(get(h(2),'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude histogram from legend

legend('fodorMI','OMES','McBASC','ELSC','fodorSCA','ramaSCA');
title('Non-MI Methods','FontSize',14,'FontWeight','n');
xlabel('Kendall corr. of method zscores to totCOV zscores for true covarions');
ylabel('No.of occurrences in 100 trials');


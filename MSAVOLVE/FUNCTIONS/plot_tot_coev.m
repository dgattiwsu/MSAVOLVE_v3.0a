function [ tot_coev_plot ] = ...
    plot_tot_coev( coev_mat,covar_vec,msa_rel_entropy,npos,hrf )

%   This function plots the coevolution matrix together with the boundaries
%   of the recombination zone. The lower inset shows the relative entropy.

Coevolution_Total = figure; 
    set(Coevolution_Total,'Units','normalized','Position',[0 0.1 0.3937 0.7 ],...
    'Name','Coevolution_Total'); clf;
    axes1 = axes('Parent',Coevolution_Total,...
    'Position',[0.15 0.25 0.7 0.7]);
        plot_coev_noaxes(coev_mat,covar_vec);figure(gcf)
        hold on
        hrfy = [1 npos];
        for i = 1:length(hrf)
        hrfx = [hrf(i) hrf(i)];
        plot(hrfx,hrfy,'--y','LineWidth',1)
        hold on
        end

        hrfx = [1 npos];
        for i = 1:length(hrf)
        hrfy = [hrf(i) hrf(i)];
        plot(hrfx,hrfy,'--r','LineWidth',1)
        hold on
        end    
     hold off
        
    axes2 = axes('Parent',Coevolution_Total,...
    'Position',[0.15 0.1 0.7 0.15]);        
        bar(msa_rel_entropy,'DisplayName','MSA_rel_entropy');figure(gcf)
        xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
        ylabel('Rel. H','FontSize',14,'FontWeight','n');
        set(gca,'Xlim',[1 npos]);

end


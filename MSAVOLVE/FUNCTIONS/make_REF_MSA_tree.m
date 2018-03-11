function [ REF_tree,COMP_tree ] = make_REF_MSA_tree( REF_nmsa, COMP_nmsa )
% This function creates a side by side figure of the phylogenetic trees of
% the REF and evolved MSA in equaldaylight radial form useful to identify
% clusters. The figure is made by creating two temporary figures and then 
% transfering the axes to a single one. This is very ackward; if you know
% how to do it in a more elegant way, please tell me.

 REF_smsa = nmsa_to_smsa(REF_nmsa);
 REF_dist = seqpdist(REF_smsa);
 REF_tree = seqlinkage(REF_dist);    

 COMP_smsa = nmsa_to_smsa(COMP_nmsa);
 COMP_dist = seqpdist(COMP_smsa);
 COMP_tree = seqlinkage(COMP_dist);    
    
REF_MSA_tree_comp = figure; 
set(REF_MSA_tree_comp,'Units','normalized','Position',[0.5 0.5 0.5 0.4],'Name',...
    'REF MSA Phylogenetic Trees');clf;
 [~,j,~] = ...
     cluster(REF_tree,[],'criterion','gain','MaxClust',6);    
 h = plot(REF_tree,'Type','equaldaylight','Orientation','left',...
      'TerminalLabels','false','LeafLabel','false');
 title('REF phylogenetic tree','FontSize',14,'FontWeight','n');figure(gcf);
 figureHandle = gcf;set(figureHandle,'Name','tempplot');
 set(h.axes,'Parent',REF_MSA_tree_comp,...
    'Position',[0.05 0.1 0.425 0.8])
 set(h.BranchLines(j==1),'Color','b')
 set(h.BranchLines(j==2),'Color','r')
 set(h.BranchLines(j==3),'Color','g')
 set(h.BranchLines(j==4),'Color','c')
 set(h.BranchLines(j==5),'Color','y')
 set(h.BranchLines(j==6),'Color','k')
 set(h.BranchLines(j>6),'Color','m')
 close tempplot

 [~,j,~] = ...
     cluster(COMP_tree,[],'criterion','gain','MaxClust',6);    
 h = plot(COMP_tree,'Type','equaldaylight','Orientation','left',...
      'TerminalLabels','false','LeafLabel','false');
 title('Synthetic MSA phylogenetic tree','FontSize',14,'FontWeight','n');figure(gcf)
 figureHandle = gcf;set(figureHandle,'Name','tempplot');
 set(h.axes,'Parent',REF_MSA_tree_comp,...
    'Position',[0.525 0.1 0.425 0.8])
 set(h.BranchLines(j==1),'Color','b')
 set(h.BranchLines(j==2),'Color','r')
 set(h.BranchLines(j==3),'Color','g')
 set(h.BranchLines(j==4),'Color','c')
 set(h.BranchLines(j==5),'Color','y')
 set(h.BranchLines(j==6),'Color','k')
 set(h.BranchLines(j>6),'Color','m')
 close tempplot

end


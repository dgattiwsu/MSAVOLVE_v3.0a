function [ smsa_tree ] = make_tree( nmsa )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

smsa = nmsa_to_smsa(nmsa);
smsa_dist = seqpdist(smsa,'Method','p-distance');
smsa_tree = seqlinkage(smsa_dist);
        
 [~,j,~] = ...
     cluster(smsa_tree,[],'criterion','gain','MaxClust',6);    
 h = plot(smsa_tree,'Orientation','left','TerminalLabels','false'); 
 % h = plot(smsa_tree,'Type','equaldaylight','Orientation','left',...
 %     'TerminalLabels','false','LeafLabel','false');

 set(h.BranchLines(j==1),'Color','b')
 set(h.BranchLines(j==2),'Color','r')
 set(h.BranchLines(j==3),'Color','g')
 set(h.BranchLines(j==4),'Color','c')
 set(h.BranchLines(j==5),'Color','y')
 set(h.BranchLines(j==6),'Color','k')
 set(h.BranchLines(j>6),'Color','m')

end


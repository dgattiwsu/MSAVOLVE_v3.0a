H3 = clustergram(nantozero(COV),'Cluster',1);
H3 = clustergram(nantozero(COV),'Standardize',2,...
'Cluster',1,'Symmetric','false','Colormap','jet','DisplayRange',1545);
H4 = plot(H3);
hold on
hrfx = [1 280];
edges = length(hrf)-1;
for i = 2:edges
    hrfy = [hrf(i) hrf(i)];
    plot(hrfx,hrfy,'--y')
    hold on
    plot(hrfy,hrfx,'--y')
end    
    
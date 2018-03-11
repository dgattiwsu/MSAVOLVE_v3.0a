%%
addpath(genpath('../../../MSAVOLVE_v2.0a'));
addpath(genpath('../../../NETWORKS_CONN_MIT/code'));

%%
matfile = 'TRIADS_TRANSITIVITY_1JZW';
figure_file = 'TRIADS_TRANSITIVITY_1JZW_FIG';

%% TRANSITIVITY BY TRIADS 
% constant minimum sequence distance --------------------------------------

method_cell = {'ZPX2','md3_ZPX2','md4_ZPX2','slPSICOV','plmDCA_ZPX2',...
    'GREMLIN_ZPX2','hpPCA_ZPX2'};

ncov = ncols;
radius = 8;
near = 7;

[~,~,sorted_ZPX2,...
    sorted_md3_ZPX2,sorted_md4_ZPX2,sorted_slPSICOV] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,ZPX2,md3_ZPX2,...
    md4_ZPX2,slPSICOV,radius,near,ncov,0);
    
[~,~,sorted_plmDCA_ZPX2,sorted_GREMLIN_ZPX2,sorted_hpPCA_ZPX2,...
    sorted_MERGE] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,plmDCA_ZPX2,GREMLIN_ZPX2,...
    hpPCA_ZPX2,MERGE,radius,near,ncov,0);


step = floor(ncols/5);
cutoff = 1:step:ncols*3;
ncutoff = length(cutoff) - 1;

nmethods = length(method_cell);
if matlabpool('size') > 0
    matlabpool close
end

matlabpool(12)

% pot_triads_mat = nan(ncutoff,nmethods);
% con_triads_mat = nan(ncutoff,nmethods);
pot_triads_mat2 = nan(ncutoff,nmethods);
con_triads_mat2 = nan(ncutoff,nmethods);

for N = 1:nmethods
    
coev_mat = eval(method_cell{N});
coev_mat_u = triu(coev_mat,near);
coev_mat_l = tril(coev_mat,-near);
coev_mat = coev_mat_u + coev_mat_l;

% Here we retrieve the sorted matrix.
sorted_mat = eval(['sorted_' method_cell{N}]);

parfor n = 1:ncutoff
    
cutoff_low = cutoff(n+1);
    
L_value_low = sorted_mat(cutoff_low,1);
L_coev_mat = coev_mat;
zero_ind = L_coev_mat < L_value_low;
L_coev_mat(zero_ind) = 0;
pos_ind = L_coev_mat ~= 0;
L_coev_mat(pos_ind) = 1;

% % Here we calculate the number of potential triads. First we determine all
% % the unique nodes.
% nodes_vec = unique(sorted_mat(1:cutoff_low,2:3));
% nnodes = numel(nodes_vec);
% temp_mat = sorted_mat(1:cutoff_low,2:3);
% 
% temp = [temp_mat ; temp_mat(:,[2 1])];
% list1 = zeros(nnodes^3,3);
% list2 = zeros(nnodes^3,3);
% count = 0;
% 
% for ii = 1:nnodes
%     i = nodes_vec(ii);
%     for jj = ii+1:nnodes
%         j = nodes_vec(jj);        
%         for kk = ii+1:nnodes
%             k = nodes_vec(kk);            
%             count = count + 1;
%             
%             [r,~,~] = find(temp(:,1) == i & temp(:,2) == j);
%             pair1 = temp(r,:);
%             [r,~,~] = find(temp(:,1) == i & temp(:,2) == k);
%             pair2 = temp(r,:);
%             if ~isempty(pair1) && ~isempty(pair2);
%                 triad = unique([pair1 pair2]);
%                 if numel(triad) > 2
%                 list1(count,:) = triad;
%                 end
%             end
%                         
%             [r,~,~] = find(temp(:,1) == j & temp(:,2) == k);
%             pair2 = temp(r,:);
%             if ~isempty(pair1) && ~isempty(pair2)
%                 triad = unique([pair1 pair2]);
%                 if numel(triad) > 2;
%                 list2(count,:) = triad;
%                 end
%             end
%             
%             
%         end
%     end
% end
% list = [list1; list2];
% unique_list = unique(list,'rows');
% 
% pot_triads_mat(n,N) = size(unique_list,1) - 1;
% con_triads_mat(n,N) = trace(L_coev_mat^3)/6;

% Alternative calculation of potential and connected triads based on MIT 
% Matlab toolbox:
pot_triads_mat2(n,N) = num_conn_triples(L_coev_mat);
con_triads_mat2(n,N) = loops3(L_coev_mat);

end

% clear list1 list2
end

% tr_mat = con_triads_mat./pot_triads_mat;
tr_mat2 = con_triads_mat2./pot_triads_mat2;


matlabpool close

% constant network size|variable sequence distance ------------------------

seqdist = [2:21];
nseqdist = numel(seqdist);
cutoff_low = ncols;

% cl_pot_triads_mat = nan(nseqdist,nmethods);
% cl_con_triads_mat = nan(nseqdist,nmethods);
cl_pot_triads_mat2 = nan(nseqdist,nmethods);
cl_con_triads_mat2 = nan(nseqdist,nmethods);

ncov = ncols;
radius = 8;

for n = 1:nseqdist
    
near = seqdist(n);

[~,~,sorted_ZPX2,...
    sorted_md3_ZPX2,sorted_md4_ZPX2,sorted_slPSICOV] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,ZPX2,md3_ZPX2,...
    md4_ZPX2,slPSICOV,radius,near,ncov,0);
    
[~,~,sorted_plmDCA_ZPX2,sorted_GREMLIN_ZPX2,sorted_hpPCA_ZPX2,...
    sorted_MERGE] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,plmDCA_ZPX2,GREMLIN_ZPX2,...
    hpPCA_ZPX2,MERGE,radius,near,ncov,0);
                 
for N = 1:nmethods
    
coev_mat = eval(method_cell{N});
coev_mat_u = triu(coev_mat,near);
coev_mat_l = tril(coev_mat,-near);
coev_mat = coev_mat_u + coev_mat_l;

% Here we retrieve the sorted matrix.
sorted_mat = eval(['sorted_' method_cell{N}]);
        
L_value_low = sorted_mat(cutoff_low,1);
L_coev_mat = coev_mat;
zero_ind = L_coev_mat < L_value_low;
L_coev_mat(zero_ind) = 0;
pos_ind = L_coev_mat ~= 0;
L_coev_mat(pos_ind) = 1;

% % Here we calculate the number of potential triads. First we determine all
% % the unique nodes.
% nodes_vec = unique(sorted_mat(1:cutoff_low,2:3));
% nnodes = numel(nodes_vec);
% temp_mat = sorted_mat(1:cutoff_low,2:3);
% 
% temp = [temp_mat ; temp_mat(:,[2 1])];
% list1 = zeros(nnodes^3,3);
% list2 = zeros(nnodes^3,3);
% count = 0;
% 
% for ii = 1:nnodes
%     i = nodes_vec(ii);
%     for jj = ii+1:nnodes
%         j = nodes_vec(jj);        
%         for kk = ii+1:nnodes
%             k = nodes_vec(kk);            
%             count = count + 1;
%             
%             [r,~,~] = find(temp(:,1) == i & temp(:,2) == j);
%             pair1 = temp(r,:);
%             [r,~,~] = find(temp(:,1) == i & temp(:,2) == k);
%             pair2 = temp(r,:);
%             if ~isempty(pair1) && ~isempty(pair2);
%                 triad = unique([pair1 pair2]);
%                 if numel(triad) > 2
%                 list1(count,:) = triad;
%                 end
%             end
%                         
%             [r,~,~] = find(temp(:,1) == j & temp(:,2) == k);
%             pair2 = temp(r,:);
%             if ~isempty(pair1) && ~isempty(pair2)
%                 triad = unique([pair1 pair2]);
%                 if numel(triad) > 2;
%                 list2(count,:) = triad;
%                 end
%             end
%             
%             
%         end
%     end
% end
% list = [list1; list2];
% unique_list = unique(list,'rows');
% 
% clear list1 list2
% %--------------------------------------------------------------------------
% 
% cl_pot_triads_mat(n,N) = size(unique_list,1) - 1;
% cl_con_triads_mat(n,N) = trace(L_coev_mat^3)/6;

% Alternative calculation of potential and connected triads based on MIT 
% Matlab toolbox:
cl_pot_triads_mat2(n,N) = num_conn_triples(L_coev_mat);
cl_con_triads_mat2(n,N) = loops3(L_coev_mat);

end
end

% cl_tr_mat = cl_con_triads_mat./cl_pot_triads_mat;
cl_tr_mat2 = cl_con_triads_mat2./cl_pot_triads_mat2;

% matlabpool close

%%

tr_mat0 = tr_mat2;
cl_tr_mat0 = cl_tr_mat2;

TRIADS_TRANSITIVITY = figure; 
    	set(TRIADS_TRANSITIVITY,'Units','normalized','Position',[0 0.2 0.9 0.4 ],...
    	'Name','Triads Transitivity'); clf;
    
subplot1 = subplot(1,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

step = floor(ncols/5);
cutoff = 1:step:ncols*3;

% ZPX2
mean_trans = tr_mat0(:,1);
% mean_trans_std = mean_trans_std_dist(:,1);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h1 = plot(xi(1:end),yi(1:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
mean_trans = tr_mat0(:,2);
% mean_trans_std = mean_trans_std_dist(:,2);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h2 = plot(xi(1:end),yi(1:end),'Color','r','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
mean_trans = tr_mat0(:,3);
% mean_trans_std = mean_trans_std_dist(:,3);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h3 = plot(xi(1:end),yi(1:end),'Color','k','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
mean_trans = tr_mat0(:,4);
% mean_trans_std = mean_trans_std_dist(:,4);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h4 = plot(xi(1:end),yi(1:end),'Color','b','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
mean_trans = tr_mat0(:,5);
% mean_trans_std = mean_trans_std_dist(:,5);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h5 = plot(xi(1:end),yi(1:end),'Color','g','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
mean_trans = tr_mat0(:,6);
% mean_trans_std = mean_trans_std_dist(:,6);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h6 = plot(xi(1:end),yi(1:end),'Color','c','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
mean_trans = tr_mat0(:,7);
% mean_trans_std = mean_trans_std_dist(:,7);

x = cutoff(2:end)';
xi = [cutoff(2):cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h7 = plot(xi(1:end),yi(1:end),'Color','m','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[0 cutoff(end)+cutoff(2)], 'YLim',[-0.005 0.185])
%     legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
%         'plmDCA','GREMLIN','hpPCA',...
%         'Location','NorthWest');
%     legend('boxoff');
    
% Label axes
xlabel( 'Number of top pairs included' );
ylabel( 'Average transitivity' );

vline (ncols,'-y');

%--------------------------------------------------------------------------
subplot2 = subplot(1,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

% ZPX2
mean_trans = cl_tr_mat0(:,1);
% mean_trans_std = mean_trans_std_dist(:,1);

cutoff = [0 (seqdist - 1)];

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h1 = plot(xi(1:end),yi(1:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
mean_trans = cl_tr_mat0(:,2);
% mean_trans_std = mean_trans_std_dist(:,2);

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h2 = plot(xi(1:end),yi(1:end),'Color','r','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
mean_trans = cl_tr_mat0(:,3);
% mean_trans_std = mean_trans_std_dist(:,3);

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h3 = plot(xi(1:end),yi(1:end),'Color','k','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
mean_trans = cl_tr_mat0(:,4);
% mean_trans_std = mean_trans_std_dist(:,4);

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h4 = plot(xi(1:end),yi(1:end),'Color','b','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
mean_trans = cl_tr_mat0(:,5);
% mean_trans_std = mean_trans_std_dist(:,5);

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h5 = plot(xi(1:end),yi(1:end),'Color','g','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
mean_trans = cl_tr_mat0(:,6);
% mean_trans_std = mean_trans_std_dist(:,6);

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h6 = plot(xi(1:end),yi(1:end),'Color','c','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
mean_trans = cl_tr_mat0(:,7);
% mean_trans_std = mean_trans_std_dist(:,7);

x = cutoff(2:end)';
xi = [cutoff(2):0.1:cutoff(end)]';
yi = (interp1(x,mean_trans(1:end),xi,'pchip'))';

h7 = plot(xi(1:end),yi(1:end),'Color','m','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(1:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(1:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[0 cutoff(end)+cutoff(2)],'YLim',[0.01 0.135])
%     legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
%         'plmDCA','GREMLIN','hpPCA',...
%         'Location','NorthEast');
%     legend('boxoff');
    
% Label axes
xlabel( 'Minimum sequence distance' );
ylabel( 'Average transitivity' );

% vline (ncols,'-y');

%% 
saveas(gcf,figure_file,'fig');

%%
save(matfile);

%%
close all
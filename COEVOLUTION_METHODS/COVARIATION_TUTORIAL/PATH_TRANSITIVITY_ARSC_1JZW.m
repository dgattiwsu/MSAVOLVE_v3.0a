%%
addpath(genpath('../../../../MSAVOLVE_v2.0a'));

%%
matfile = 'PATH_TRANSITIVITY_1JZW';
figure_file = 'PATH_TRANSITIVITY_1JZW_FIG';

%% PATH TRANSITIVITY

method_cell = {'ZPX2','md3_ZPX2','md4_ZPX2','slPSICOV','plmDCA_ZPX2',...
    'GREMLIN_ZPX2','hpPCA_ZPX2'};

step = floor(ncols/5);
cutoff = 1:step:ncols*3;
ncutoff = length(cutoff) - 1;

nmethods = length(method_cell);
global_mean_step_length = zeros(ncutoff,nmethods);
global_path_length = zeros(ncutoff,nmethods);
corr_coev_dist = zeros(ncutoff,nmethods);
corr_coev_mean_step_length = zeros(ncutoff,nmethods);
corr_coev_path_length = zeros(ncutoff,nmethods);

if matlabpool('size') > 0
    matlabpool close
end

% if nprocs > 1
matlabpool(12)
% end

for N = 1:nmethods
    
test_mat = eval(method_cell{N});
coev_mat = nantozero(triu(test_mat,0));
coev_vec = coev_mat(:);
[sorted_coev_vec,sorted_coev_vec_ind] = sort(coev_vec,'descend');


parfor k = 1:ncutoff
    
% cutoff_high = cutoff(k);    
cutoff_low = cutoff(k+1);
    
% L_value_high = sorted_coev_vec(cutoff_high);
L_value_low = sorted_coev_vec(cutoff_low);
coev_mat = test_mat;
% select_ind = coev_mat <= L_value_high & coev_mat > L_value_low;
zero_ind = coev_mat < L_value_low;
L_coev_mat = nantozero(coev_mat);
L_coev_mat(zero_ind) = 0;
L_coev_mat_full = L_coev_mat;
% pos_ind selects all the pairs retained for each range
pos_ind = L_coev_mat ~= 0;
L_coev_mat(pos_ind) = 1;
sL_coev_mat = sparse(L_coev_mat);
path_length = zeros(ncols,ncols);
path_nsteps = zeros(ncols,ncols);
mean_step_length = zeros(ncols,ncols);

for i = 1:ncols
    for j = i+1:ncols
    
        if L_coev_mat(i,j) == 0
            path_length(i,j) = 0;
            path_nsteps(i,j) = 0;
            mean_step_length(i,j) = 0;
        else
            L_coev_mat_n = L_coev_mat;
            % Here we remove the pair under analysis from the distance 
            % matrix, otherwise the determination of the shortest path
            % would recognize a single step
            L_coev_mat_n(i,j) = 0;
            L_coev_mat_n(j,i) = 0;
            n_dist_ind = L_coev_mat_n == 0;
            dist_mat_n = c_distances;
            dist_mat_n(n_dist_ind) = 0;
            s_dist_mat_n = sparse(dist_mat_n);
            [dist_n, path_n] = graphshortestpath(sparse(dist_mat_n),i,j);
            path_length(i,j) = dist_n;
            path_nsteps(i,j) = numel(path_n) - 1;
            mean_step_length(i,j) = path_length(i,j)/path_nsteps(i,j);
        end

    end
end

% Here we analyze the mean length of the steps along the shortest path
% between two positions
inf_ind = isinf(mean_step_length);
mean_step_length(zero_ind) = 0;
mean_step_length(inf_ind) = 0;
pos_ind = find(mean_step_length);

% Here we analyze the total length of the steps along the shortest path
% between two positions
path_length(zero_ind) = 0;
path_length(inf_ind) = 0;

% Here we take the mean of all the pairs in the selected range
global_mean_step_length(k,N) = nanmean(mean_step_length(pos_ind));
global_path_length(k,N) = nanmean(path_length(pos_ind));

distance_vec = c_distances(pos_ind);
mean_step_length_vec = mean_step_length(pos_ind); 
path_length_vec = path_length(pos_ind); 
coev_score_vec = L_coev_mat_full(pos_ind);
[coev_score_vec,sort_ind] = sort(coev_score_vec,'descend');
distance_vec = distance_vec(sort_ind);
mean_step_length_vec = mean_step_length_vec(sort_ind);
path_length_vec = path_length_vec(sort_ind);

if ~isempty(pos_ind)
corr_coev_dist(k,N) = corr(distance_vec,coev_score_vec);
corr_coev_mean_step_length(k,N) = corr(mean_step_length_vec,coev_score_vec);
corr_coev_path_length(k,N) = corr(path_length_vec,coev_score_vec);
else
corr_coev_dist(k,N) = NaN;
corr_mean_step_length(k,N) = NaN;    
corr_coev_path_length(k,N) = NaN;    
end

end

end

L_cutoff = cutoff/ncols;
[N,D] = rat(L_cutoff,0.1);

matlabpool close


%%
PATH_TRANSITIVITY = figure; 
    	set(PATH_TRANSITIVITY,'Units','normalized','Position',[0 0.2 0.9 0.8 ],...
    	'Name','Path Transitivity'); clf;

subplot1 = subplot(2,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

% ZPX2
mean_trans = global_mean_step_length(:,1);
% mean_trans_std = mean_trans_std_dist(:,1);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h1 = plot(xi(1:end),yi(1:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
mean_trans = global_mean_step_length(:,2);
% mean_trans_std = mean_trans_std_dist(:,2);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h2 = plot(xi(1:end),yi(1:end),'Color','r','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
mean_trans = global_mean_step_length(:,3);
% mean_trans_std = mean_trans_std_dist(:,3);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h3 = plot(xi(1:end),yi(1:end),'Color','k','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
mean_trans = global_mean_step_length(:,4);
% mean_trans_std = mean_trans_std_dist(:,4);

x = cutoff(4:end)';
xi = [cutoff(4):cutoff(end)]';
yi = (interp1(x,mean_trans(3:end),xi,'pchip'))';

h4 = plot(xi(1:end),yi(1:end),'Color','b','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(3:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
mean_trans = global_mean_step_length(:,5);
% mean_trans_std = mean_trans_std_dist(:,5);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h5 = plot(xi(1:end),yi(1:end),'Color','g','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
mean_trans = global_mean_step_length(:,6);
% mean_trans_std = mean_trans_std_dist(:,6);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h6 = plot(xi(1:end),yi(1:end),'Color','c','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
mean_trans = global_mean_step_length(:,7);
% mean_trans_std = mean_trans_std_dist(:,7);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h7 = plot(xi(1:end),yi(1:end),'Color','m','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[cutoff(2) cutoff(end)+cutoff(2)], 'YLim',[5.5 13.5])
    legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
        'plmDCA','GREMLIN','hpPCA',...
        'Location','NorthEast');
    legend('boxoff');
    
% Label axes
xlabel( 'Number of pairs included' );
ylabel( 'Mean step length' );

vline (ncols,'-y');

%--------------------------------------------------------------------------
subplot2 = subplot(2,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

% ZPX2
mean_trans = global_path_length(:,1);
% mean_trans_std = mean_trans_std_dist(:,1);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h1 = plot(xi(1:end),yi(1:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
mean_trans = global_path_length(:,2);
% mean_trans_std = mean_trans_std_dist(:,2);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h2 = plot(xi(1:end),yi(1:end),'Color','r','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
mean_trans = global_path_length(:,3);
% mean_trans_std = mean_trans_std_dist(:,3);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h3 = plot(xi(1:end),yi(1:end),'Color','k','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
mean_trans = global_path_length(:,4);
% mean_trans_std = mean_trans_std_dist(:,4);

x = cutoff(4:end)';
xi = [cutoff(4):cutoff(end)]';
yi = (interp1(x,mean_trans(3:end),xi,'pchip'))';

h4 = plot(xi(1:end),yi(1:end),'Color','b','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(3:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
mean_trans = global_path_length(:,5);
% mean_trans_std = mean_trans_std_dist(:,5);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h5 = plot(xi(1:end),yi(1:end),'Color','g','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
mean_trans = global_path_length(:,6);
% mean_trans_std = mean_trans_std_dist(:,6);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h6 = plot(xi(1:end),yi(1:end),'Color','c','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
mean_trans = global_path_length(:,7);
% mean_trans_std = mean_trans_std_dist(:,7);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h7 = plot(xi(1:end),yi(1:end),'Color','m','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[cutoff(2) cutoff(end)+cutoff(2)],'YLim',[10 126])
%     legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
%         'plmDCA','GREMLIN','hpPCA',...
%         'Location','NorthEast');
%     legend('boxoff');
    
% Label axes
xlabel( 'Number of pairs included' );
ylabel( 'Path length' );

vline (ncols,'-y');

%--------------------------------------------------------------------------
subplot3 = subplot(2,2,3,'Parent',figure(gcf));
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'all');

% ZPX2
mean_trans = corr_coev_mean_step_length(:,1);
% mean_trans_std = mean_trans_std_dist(:,1);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h1 = plot(xi(1:end),yi(1:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
mean_trans = corr_coev_mean_step_length(:,2);
% mean_trans_std = mean_trans_std_dist(:,2);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h2 = plot(xi(1:end),yi(1:end),'Color','r','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
mean_trans = corr_coev_mean_step_length(:,3);
% mean_trans_std = mean_trans_std_dist(:,3);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h3 = plot(xi(1:end),yi(1:end),'Color','k','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
mean_trans = corr_coev_mean_step_length(:,4);
% mean_trans_std = mean_trans_std_dist(:,4);

x = cutoff(4:end)';
xi = [cutoff(4):cutoff(end)]';
yi = (interp1(x,mean_trans(3:end),xi,'pchip'))';

h4 = plot(xi(1:end),yi(1:end),'Color','b','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(3:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
mean_trans = corr_coev_mean_step_length(:,5);
% mean_trans_std = mean_trans_std_dist(:,5);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h5 = plot(xi(1:end),yi(1:end),'Color','g','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
mean_trans = corr_coev_mean_step_length(:,6);
% mean_trans_std = mean_trans_std_dist(:,6);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h6 = plot(xi(1:end),yi(1:end),'Color','c','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
mean_trans = corr_coev_mean_step_length(:,7);
% mean_trans_std = mean_trans_std_dist(:,7);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h7 = plot(xi(1:end),yi(1:end),'Color','m','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[cutoff(2) cutoff(end)+cutoff(2)],'YLim',[-0.33 0.18])
%     legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
%         'plmDCA','GREMLIN','hpPCA',...
%         'Location','NorthEast');
%     legend('boxoff');
    
% Label axes
xlabel( 'Number of pairs included' );
ylabel( 'Coev. score/mean step length correlation' );

vline (ncols,'-y');

%--------------------------------------------------------------------------
subplot4 = subplot(2,2,4,'Parent',figure(gcf));
box(subplot4,'on');
grid(subplot4,'on');
hold(subplot4,'all');

% ZPX2
mean_trans = corr_coev_path_length(:,1);
% mean_trans_std = mean_trans_std_dist(:,1);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h1 = plot(xi(1:end),yi(1:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
mean_trans = corr_coev_path_length(:,2);
% mean_trans_std = mean_trans_std_dist(:,2);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h2 = plot(xi(1:end),yi(1:end),'Color','r','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
mean_trans = corr_coev_path_length(:,3);
% mean_trans_std = mean_trans_std_dist(:,3);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h3 = plot(xi(1:end),yi(1:end),'Color','k','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
mean_trans = corr_coev_path_length(:,4);
% mean_trans_std = mean_trans_std_dist(:,4);

x = cutoff(4:end)';
xi = [cutoff(4):cutoff(end)]';
yi = (interp1(x,mean_trans(3:end),xi,'pchip'))';

h4 = plot(xi(1:end),yi(1:end),'Color','b','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(3:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
mean_trans = corr_coev_path_length(:,5);
% mean_trans_std = mean_trans_std_dist(:,5);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h5 = plot(xi(1:end),yi(1:end),'Color','g','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
mean_trans = corr_coev_path_length(:,6);
% mean_trans_std = mean_trans_std_dist(:,6);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h6 = plot(xi(1:end),yi(1:end),'Color','c','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
mean_trans = corr_coev_path_length(:,7);
% mean_trans_std = mean_trans_std_dist(:,7);

x = cutoff(3:end)';
xi = [cutoff(3):cutoff(end)]';
yi = (interp1(x,mean_trans(2:end),xi,'pchip'))';

h7 = plot(xi(1:end),yi(1:end),'Color','m','LineStyle','-','LineWidth',1.0);
hold on
% errorbar(x,mean_trans(2:end),mean_trans_std(1:end),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x,mean_trans(2:end),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[cutoff(2) cutoff(end)+cutoff(2)],'YLim',[-0.42 0.2])
%     legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
%         'plmDCA','GREMLIN','hpPCA',...
%         'Location','NorthEast');
%     legend('boxoff');
    
% Label axes
xlabel( 'Number of pairs included' );
ylabel( 'Coev. score/path length correlation' );

vline (ncols,'-y');


%% 
saveas(gcf,figure_file,'fig');

%%
save(matfile);

%%
close all
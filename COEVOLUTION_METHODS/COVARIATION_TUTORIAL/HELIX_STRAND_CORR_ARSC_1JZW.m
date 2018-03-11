%%
addpath(genpath('../../../MSAVOLVE_v2.0a'));

%%
matfile = 'SECONDARY_STRUCTURE_1JZW';
figure_file = 'SECONDARY_STRUCTURE_1JZW_FIG';

%% HELIX CORRELATION

nhelix = size(pdbstruct_1.Helix,2);
nhelix_A = 0;
for i = 1:nhelix
    if strcmp(pdbstruct_1.Helix(1,i).initChainID,'A')
        nhelix_A = nhelix_A + 1;
    end
end
helix_cell = cell(nhelix_A,1);
for i = 1:nhelix_A
    helix_cell{i} = [pdbstruct_1.Helix(1,i).initSeqNum:pdbstruct_1.Helix(1,i).endSeqNum];
    
end
total_vec = cat(2,helix_cell{:});
n_total_vec = size(total_vec,2);

% If there is a difference between the numbering of the coevolution matrix 
% and the numbering of the pdb file find the sequence shift value and the 
% reference sequence identity in the smsa structure

% sequence_shift = START - PDB_START;
% special correction for just ArsC
sequence_shift = -2;

for i = 1:nhelix_A
    helix_cell{i} = helix_cell{i} + sequence_shift;
end

nseqs = size(nmsa,1);

for i = 1:nseqs
    header = REF_smsa(i,1).Header;
    if strcmp('P08692|ARSC1_ECOLX',header)
        seq_ind = i;
    end
end

% Alternatively
for i = 1:nseqs
    sequence = REF_smsa(i,1).Sequence;
    if strcmp(pdbstruct_1.Sequence.Sequence,sequence)
        seq_ind = i;
    end
end

% Check on sequence-coevolution matrix alignment
for i = 1:nhelix_A
    cmsa_seq = cmsa(seq_ind,helix_cell{i});
    pdb_seq_start = pdbstruct_1.Helix(1,i).initSeqNum;
    pdb_seq_end = pdbstruct_1.Helix(1,i).endSeqNum;
    pdb_seq = pdbstruct_1.Sequence.Sequence(pdb_seq_start:pdb_seq_end); 
    if strcmp(cmsa_seq,pdb_seq) 
        display('OK')
    else
        display('ERROR')
    end
end

%%
method_cell = {'ZPX2','md3_ZPX2','md4_ZPX2','slPSICOV','plmDCA_ZPX2',...
    'GREMLIN_ZPX2','hpPCA_ZPX2'};

ac_helix_cell = cell(length(method_cell),1);
ac_helix_std_cell = cell(length(method_cell),1);
lastind_hvec = zeros(length(method_cell),1);

for k = 1: length(method_cell)

coev_mat = eval(method_cell{k});

ac_helix = zeros(50,1);
ac_helix_mat = NaN(50,n_total_vec);
ac_helix_sum = zeros(50,1);
for j = 1:length(helix_cell)
    h = helix_cell{j};
    ac_h = zeros(length(h),1);
    last = h(end);
    vec_length = length(h);
    for i = 1: vec_length
        first = h(i);
        temp_vec = nantozero(coev_mat(first,first:last)');
        % temp_vec = [0;zscore(temp_vec(2:end))];
        % temp_vec = zscore(temp_vec);
        vec_pad = i-1;
        temp_vec_ind = find(temp_vec);
        ac_helix_sum(temp_vec_ind) = ac_helix_sum(temp_vec_ind) +1;
        current_col = ac_helix_sum(2);
        ac_helix_mat(temp_vec_ind,current_col) = temp_vec(temp_vec_ind);
        temp_vec = [temp_vec;zeros(vec_pad,1)];
        ac_h = ac_h + temp_vec;
        % ac_helix_sum(temp_vec_ind) = ac_helix_sum(temp_vec_ind) +1;        
    end
    ac_pad = zeros(50-vec_length,1);
    ac_helix = ac_helix + [ac_h;ac_pad];
end
ac_helix_mean = nantozero(nanmean(ac_helix_mat,2));
ac_helix_std = nantozero(nanstd(ac_helix_mat,1,2));

ac_helix = ac_helix_mean;
ac_helix_ind = find(ac_helix);
lastind_hvec(k) = ac_helix_ind(end);
% ac_helix(ac_helix_ind) = ac_helix(ac_helix_ind)./ac_helix_sum(ac_helix_ind);
ac_helix(ac_helix_ind) = zscore(ac_helix(ac_helix_ind));
ac_helix_std(ac_helix_ind) = abs(zscore(ac_helix_std(ac_helix_ind)));

ac_helix_cell{k} = ac_helix;
ac_helix_std_cell{k} = ac_helix_std;
end


%% STRAND CORRELATION

nstrands = size(pdbstruct_1.Sheet,2);
nstrands_A = 0;
for i = 1:nstrands
    if strcmp(pdbstruct_1.Sheet(1,i).initChainID,'A')
        nstrands_A = nstrands_A + 1;
    end
end
strand_cell = cell(nstrands_A,1);
for i = 1:nstrands_A
    strand_cell{i} = [pdbstruct_1.Sheet(1,i).initSeqNum:pdbstruct_1.Sheet(1,i).endSeqNum];
end
total_vec = cat(2,strand_cell{:});
n_total_vec = size(total_vec,2);

% If there is a difference between the numbering of the coevolution matrix 
% and the numbering of the pdb file find the sequence shift value and the 
% reference sequence identity in the smsa structure

% sequence_shift = START - PDB_START;
% special correction for just ArsC
sequence_shift = -2;

for i = 1:nstrands_A
    strand_cell{i} = strand_cell{i} + sequence_shift;
end

method_cell = {'ZPX2','md3_ZPX2','md4_ZPX2','slPSICOV','plmDCA_ZPX2',...
    'GREMLIN_ZPX2','hpPCA_ZPX2'};

ac_strand_cell = cell(length(method_cell),1);
ac_strand_std_cell = cell(length(method_cell),1);
lastind_svec = zeros(length(method_cell),1);

for k = 1: length(method_cell)

coev_mat = eval(method_cell{k});

ac_strand = zeros(50,1);
ac_strand_mat = NaN(50,n_total_vec);
ac_strand_sum = zeros(50,1);
for j = 1:length(strand_cell)
    h = strand_cell{j};
    ac_h = zeros(length(h),1);
    last = h(end);
    vec_length = length(h);
    for i = 1: vec_length
        first = h(i);
        temp_vec = nantozero(coev_mat(first,first:last)');
        % temp_vec = [0;zscore(temp_vec(2:end))];
        % temp_vec = zscore(temp_vec);
        vec_pad = i-1;
        temp_vec_ind = find(temp_vec);
        ac_strand_sum(temp_vec_ind) = ac_strand_sum(temp_vec_ind) +1;
        current_col = ac_strand_sum(2);
        ac_strand_mat(temp_vec_ind,current_col) = temp_vec(temp_vec_ind);
        temp_vec = [temp_vec;zeros(vec_pad,1)];
        ac_h = ac_h + temp_vec;
        % ac_strand_sum(temp_vec_ind) = ac_strand_sum(temp_vec_ind) +1;        
    end
    ac_pad = zeros(50-vec_length,1);
    ac_strand = ac_strand + [ac_h;ac_pad];
end
ac_strand_mean = nantozero(nanmean(ac_strand_mat,2));
ac_strand_std = nantozero(nanstd(ac_strand_mat,1,2));

ac_strand = ac_strand_mean;
ac_strand_ind = find(ac_strand);
lastind_svec(k) = ac_strand_ind(end);
% ac_strand(ac_strand_ind) = ac_strand(ac_strand_ind)./ac_strand_sum(ac_strand_ind);
ac_strand(ac_strand_ind) = zscore(ac_strand(ac_strand_ind));
ac_strand_std(ac_strand_ind) = abs(zscore(ac_strand_std(ac_strand_ind)));

ac_strand_cell{k} = ac_strand;
ac_strand_std_cell{k} = ac_strand_std;
end


%%
SEC_STRUCT_CORR = figure; 
    	set(SEC_STRUCT_CORR,'Units','normalized','Position',[0 0.2 0.9 0.4 ],...
    	'Name','Secondary structure autocorrelation'); clf;

subplot1 = subplot(1,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

% ZPX2
ac_helix = ac_helix_cell{1};
ac_helix_std = ac_helix_std_cell{1};
lastind = lastind_hvec(1);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h1 = plot(xi(11:end),yi(11:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0);
% hold on
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')

% md3_ZPX2
ac_helix = ac_helix_cell{2};
ac_helix_std = ac_helix_std_cell{2};
lastind = lastind_hvec(2);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h2 = plot(xi(11:end),yi(11:end),'Color','r','LineStyle','-','LineWidth',1.0);
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color','r','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')

% md4_ZPX2
ac_helix = ac_helix_cell{3};
ac_helix_std = ac_helix_std_cell{3};
lastind = lastind_hvec(3);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h3 = plot(xi(11:end),yi(11:end),'Color','k','LineStyle','-','LineWidth',1.0);
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color','k','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')

% slPSICOV
ac_helix = ac_helix_cell{4};
ac_helix_std = ac_helix_std_cell{4};
lastind = lastind_hvec(4);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h4 = plot(xi(11:end),yi(11:end),'Color','b','LineStyle','-','LineWidth',1.0);
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color','b','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')

% plmDCA_ZPX2
ac_helix = ac_helix_cell{5};
ac_helix_std = ac_helix_std_cell{5};
lastind = lastind_hvec(5);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h5 = plot(xi(11:end),yi(11:end),'Color','g','LineStyle','-','LineWidth',1.0);
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color','g','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')

% GREMLIN_ZPX2
ac_helix = ac_helix_cell{6};
ac_helix_std = ac_helix_std_cell{6};
lastind = lastind_hvec(6);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h6 = plot(xi(11:end),yi(11:end),'Color','c','LineStyle','-','LineWidth',1.0);
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color','c','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')

% hpPCA
ac_helix = ac_helix_cell{7};
ac_helix_std = ac_helix_std_cell{7};
lastind = lastind_hvec(7);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_helix(1:lastind),xi,'spline');

h7 = plot(xi(11:end),yi(11:end),'Color','m','LineStyle','-','LineWidth',1.0);
% errorbar(x(2:end),ac_helix(2:lastind),ac_helix_std(2:lastind),...
%     'Color','m','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_helix(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')

set(gca, 'XLim',[0 lastind])
    legend([h1,h2,h3,h4,h5,h6,h7],'2D\_MI','3D\_MI','4D\_MI','PSICOV',...
        'plmDCA','GREMLIN','hpPCA',...
        'Location','NorthEast');
    legend('boxoff');
    
% Label axes
xlabel( 'Sequence space - Helices' );
ylabel( 'Average Normalized Covariation Score' );

vline ([3.6 7.2 10.8 14.4 18],'-y');

%--------------------------------------------------------------------------
subplot2 = subplot(1,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

% ZPX2
ac_strand = ac_strand_cell{1};
ac_strand_std = ac_strand_std_cell{1};
lastind = lastind_svec(1);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color',[0,0.5,1],'LineStyle','-','LineWidth',1.0)
% hold on
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color',[0,0.5,1],'LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor',[0,0.5,1],'MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])

% md3_ZPX2
ac_strand = ac_strand_cell{2};
ac_strand_std = ac_strand_std_cell{2};
lastind = lastind_svec(2);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color','r','LineStyle','-','LineWidth',1.0)
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color','r','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','r','MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])

% md4_ZPX2
ac_strand = ac_strand_cell{3};
ac_strand_std = ac_strand_std_cell{3};
lastind = lastind_svec(3);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color','k','LineStyle','-','LineWidth',1.0)
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color','k','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','k','MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])

% slPSICOV
ac_strand = ac_strand_cell{4};
ac_strand_std = ac_strand_std_cell{4};
lastind = lastind_svec(4);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color','b','LineStyle','-','LineWidth',1.0)
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color','b','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','b','MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])

% plmDCA_ZPX2
ac_strand = ac_strand_cell{5};
ac_strand_std = ac_strand_std_cell{5};
lastind = lastind_svec(5);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color','g','LineStyle','-','LineWidth',1.0)
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color','g','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','g','MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])

% GREMLIN_ZPX2
ac_strand = ac_strand_cell{6};
ac_strand_std = ac_strand_std_cell{6};
lastind = lastind_svec(6);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color','c','LineStyle','-','LineWidth',1.0)
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color','c','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','c','MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])

% hpPCA_ZPX2
ac_strand = ac_strand_cell{7};
ac_strand_std = ac_strand_std_cell{7};
lastind = lastind_svec(7);

x = 0:lastind-1;
xi = 0:0.1:lastind-1;
yi = interp1(x,ac_strand(1:lastind),xi,'spline');

plot(xi(11:end),yi(11:end),'Color','m','LineStyle','-','LineWidth',1.0)
% errorbar(x(2:end),ac_strand(2:lastind),ac_strand_std(2:lastind),...
%     'Color','m','LineStyle','none','LineWidth',1)
plot(x(2:end),ac_strand(2:lastind),'LineStyle','none','LineWidth',1.0,...
    'Marker','o','MarkerSize',4,...
    'MarkerEdgeColor','m','MarkerFaceColor','w')
set(gca, 'XLim',[0 lastind])
set(gca, 'XLim',[0.5 3.5])

% Label axes
xlabel( 'Sequence space - Strands' );
ylabel( 'Average Normalized Covariation Score' );

vline ([2.0],'-y');
% vline ([2.3 4.6 6.8 9.2 11.5],'-y');

%% 
saveas(gcf,figure_file,'fig');

%%
save(matfile);

%%
close all
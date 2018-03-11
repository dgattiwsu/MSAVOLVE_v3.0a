%%
addpath(genpath('../../../MSAVOLVE_v2.0a'));

%%
matfile = 'COMPARE_COEV_DISTMAT_ARSC_1JZW';
figure_file = 'COMPARE_COEV_DISTMAT_ARSC_1JZW_FIG';
cumdist_figure_file = 'COMPARE_COEV_DISTMAT_ARSC_1JZW_CUMDISTFIG';
partdist_figure_file = 'COMPARE_COEV_DISTMAT_ARSC_1JZW_PARTDISTFIG';
contact_map_file_1 = 'COMPARE_CONTACT_MAP_ARSC_1JZW_FIG_1';
contact_map_file_2 = 'COMPARE_CONTACT_MAP_ARSC_1JZW_FIG_2';

msafile = 'ARSC_95.faln';
msafile_type = 'faln';

switch msafile_type 
    case 'faln'
    [REF_smsa,REF_nmsa] = faln_to_nmsa(msafile);
    case 'aln'        
    [REF_smsa,REF_nmsa] = aln_to_nmsa(msafile);
end

smsa = REF_smsa;
nmsa = REF_nmsa;

%%
% Information about the pdb file:
pdbfile = '1JZW.pdb';
START = 3;
END = 140;
PDB_START = 3;
PDB_END = 140;
nmsa = nmsa(:,PDB_START:PDB_END);
REF_length = numel(nmsa(1,:));
cmsa = int2aa(nmsa);
ncols = size(nmsa,2);

%% ------------------------------------------------------------------------
% [ZPX2,md3_ZPX2,md4_ZPX2] = ...
%     NMSA_to_mdMI(nmsa,'GAPS','3D','FULL',0.9,1,22,0,0,3,12);
[ZPX2,md3_ZPX2,md4_ZPX2] = ...
    NMSA_to_mdMI(nmsa,'GAPS','4D','FULL',0.9,1,22,0,0,3,12);
[slPSICOV] = NMSA_to_slPSICOV(nmsa,'NOGAPS',0.9,1.0,...
    20,'SUM','NONE','fro',0.0001,0.35,0.015,100,0.0,0.0001,'BOTH',0,1);
[plmDCA,plmDCA_covar_vec] = get_nmsa_covar_vec(nmsa,30,'plmDCA_asym');
% plmDCA already has a MIP correction, so we only add a ZPX2 correction.
plmDCA_ZPX2 = MIP_to_ZPX2(plmDCA);
[GREMLIN,GREMLIN_covar_vec] = get_nmsa_covar_vec(nmsa,30,'GREMLIN');
% GREMLIN already has a MIP correction, so we only add a ZPX2 correction.
GREMLIN_ZPX2 = MIP_to_ZPX2(GREMLIN);
[hpPCA,hpPCA_covar_vec] = get_nmsa_covar_vec(nmsa,30,'hpPCA');
% hpPCA already has a MIP correction, so we only add a ZPX2 correction.
hpPCA_ZPX2 = MIP_to_ZPX2(hpPCA);

%% Gap correction 
% Here we calculate a matrix of weights to correct for the presence of
% gaps. 

gapW0 = ones(ncols,ncols);
gapW1 = correct_coevmat_forgaps(nmsa);
gapW2 = gapW1.^2;
gapW3 = gapW1.^3;

% Here we choose how to corrects the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW1,gapW2,gapW3) can be applied.
gW = gapW3;

% We only correct DCA type maps as the other methods are already corrected 
% for gaps.
plmDCA_orig_1 = plmDCA;
plmDCA_orig_2 = plmDCA_orig_1 - min(plmDCA_orig_1(:));
plmDCA = plmDCA_orig_2.*gW;

plmDCA_ZPX2_orig_1 = plmDCA_ZPX2;
plmDCA_ZPX2_orig_2 = plmDCA_ZPX2_orig_1 - min(plmDCA_ZPX2_orig_1(:));
plmDCA_ZPX2 = plmDCA_ZPX2_orig_2.*gW;

hpPCA_orig_1 = hpPCA;
hpPCA_orig_2 = hpPCA_orig_1 - min(hpPCA_orig_1(:));
hpPCA = hpPCA_orig_2.*gW;

hpPCA_ZPX2_orig_1 = hpPCA_ZPX2;
hpPCA_ZPX2_orig_2 = hpPCA_ZPX2_orig_1 - min(hpPCA_ZPX2_orig_1(:));
hpPCA_ZPX2 = hpPCA_ZPX2_orig_2.*gW;

GREMLIN_orig_1 = GREMLIN;
GREMLIN_orig_2 = GREMLIN_orig_1 - min(GREMLIN_orig_1(:));
GREMLIN = GREMLIN_orig_2.*gW;

GREMLIN_ZPX2_orig_1 = GREMLIN_ZPX2;
GREMLIN_ZPX2_orig_2 = GREMLIN_ZPX2_orig_1 - min(GREMLIN_ZPX2_orig_1(:));
GREMLIN_ZPX2 = GREMLIN_ZPX2_orig_2.*gW;

%% 
% Here we merge 3D_MI with plmDCA, hpPCA, and GREMLIN
nanind = isnan(md3_ZPX2);
plmDCA_ZPX2_nan = plmDCA_ZPX2;
plmDCA_ZPX2_nan(nanind) = NaN;
GREMLIN_ZPX2_nan = GREMLIN_ZPX2;
GREMLIN_ZPX2_nan(nanind) = NaN;

scale1 = scale_matrices(nantozero(plmDCA_ZPX2),nantozero(md3_ZPX2));
scale2 = scale_matrices(nantozero(GREMLIN_ZPX2),nantozero(md3_ZPX2));

stack = cat(3,md3_ZPX2,plmDCA_ZPX2_nan*scale1,GREMLIN_ZPX2_nan*scale2);
MERGE = nanmean(stack,3);

%%
near = 1;
near1 = near;
ncov = REF_length;
radius = 8;
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

[c_distances,sorted_c_distances,sorted_MERGE_1,...
    sorted_plmDCA_ZPX2_1,sorted_slPSICOV_1,sorted_hpPCA_ZPX2_1,pdbstruct_1] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,MERGE,plmDCA_ZPX2,...
    slPSICOV,hpPCA_ZPX2,radius,near,ncov,3,[1 1 0],[0 1 0],[0 0 1],[1 0 1]);

saveas(gcf,contact_map_file_1,'fig');
% close all

%---------------------2nd set----------------------------------------------

[~,~,sorted_ZPX2_1,sorted_md3_ZPX2_1,sorted_md4_ZPX2_1,sorted_GREMLIN_ZPX2_1] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,ZPX2,md3_ZPX2,...
    md4_ZPX2,GREMLIN_ZPX2,radius,near,ncov,3,...
    [0,0.5,1],[1 0 0],[0 0 0],[0 1 1]);

saveas(gcf,contact_map_file_2,'fig');
% close all

%%
close all

%%
near = 7;
near2 = near;
ncov = REF_length;
radius = 8;

%---------------------1st set----------------------------------------------

[c_distances,sorted_c_distances,sorted_MERGE_2,...
    sorted_plmDCA_ZPX2_2,sorted_slPSICOV_2,sorted_hpPCA_ZPX2_2,pdbstruct_1] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,MERGE,plmDCA_ZPX2,...
    slPSICOV,hpPCA_ZPX2,radius,near,ncov,0,[1 1 0],[0 1 0],[0 0 1],[1 0 1]);

%---------------------2nd set----------------------------------------------

[~,~,sorted_ZPX2_2,sorted_md3_ZPX2_2,sorted_md4_ZPX2_2,sorted_GREMLIN_ZPX2_2] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,ZPX2,md3_ZPX2,...
    md4_ZPX2,GREMLIN_ZPX2,radius,near,ncov,0,...
    [0,0.5,1],[1 0 0],[0 0 0],[0 1 1]);

%%
near = 13;
near3 = near;
ncov = REF_length;
radius = 8;

%---------------------1st set----------------------------------------------

[c_distances,sorted_c_distances,sorted_MERGE_3,...
    sorted_plmDCA_ZPX2_3,sorted_slPSICOV_3,sorted_hpPCA_ZPX2_3,pdbstruct_1] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,MERGE,plmDCA_ZPX2,...
    slPSICOV,hpPCA_ZPX2,radius,near,ncov,0,[1 1 0],[0 1 0],[0 0 1],[1 0 1]);

%---------------------2nd set----------------------------------------------

[~,~,sorted_ZPX2_3,sorted_md3_ZPX2_3,sorted_md4_ZPX2_3,sorted_GREMLIN_ZPX2_3] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,ZPX2,md3_ZPX2,...
    md4_ZPX2,GREMLIN_ZPX2,radius,near,ncov,0,...
    [0,0.5,1],[1 0 0],[0 0 0],[0 1 1]);

%%
near = 21;
near4 = near;
ncov = REF_length;
radius = 8;

%---------------------1st set----------------------------------------------

[c_distances,sorted_c_distances,sorted_MERGE_4,...
    sorted_plmDCA_ZPX2_4,sorted_slPSICOV_4,sorted_hpPCA_ZPX2_4,pdbstruct_1] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,MERGE,plmDCA_ZPX2,...
    slPSICOV,hpPCA_ZPX2,radius,near,ncov,0,[1 1 0],[0 1 0],[0 0 1],[1 0 1]);

%---------------------2nd set----------------------------------------------

[~,~,sorted_ZPX2_4,sorted_md3_ZPX2_4,sorted_md4_ZPX2_4,sorted_GREMLIN_ZPX2_4] = ...
    coev_distance_matrix_3(pdbfile,1,PDB_START,PDB_END,ZPX2,md3_ZPX2,...
    md4_ZPX2,GREMLIN_ZPX2,radius,near,ncov,0,...
    [0,0.5,1],[1 0 0],[0 0 0],[0 1 1]);


%%
save(matfile);

%%
COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.8 0.8 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;

subplot1 = subplot(2,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

ncov = REF_length;
X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_1(X,6),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_1(X,6),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_1(X,6),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_1(X,6),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_1(X,6),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_1(X,6),'-m','LineWidth',2);
semilogy(X,sorted_MERGE_1(X,6),'-y','LineWidth',2);
hold off
separation1 = num2str(near1 - 1);
distance = num2str(radius);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,0.13]);
    legend('3D\_MI','4D\_MI','PSICOV',...
        'plmDCA','GREMLIN','hpPCA','MERGED',...
        'Location','NorthWest');
% 	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
string1 = ['Number of top pairs identified by a method'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation1 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation1 ' intervening positions in the sequence '];
ylabel({string1;string2;string3});


subplot2 = subplot(2,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_2(X,6),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_2(X,6),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_2(X,6),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_2(X,6),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_2(X,6),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_2(X,6),'-m','LineWidth',2);
semilogy(X,sorted_MERGE_2(X,6),'-y','LineWidth',2);
hold off
separation2 = num2str(near2 - 1);
distance = num2str(radius);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,.16]);
string1 = ['Number of top pairs identified by a method'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation2 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation2 ' intervening positions in the sequence '];
ylabel({string1;string2;string3});


subplot3 = subplot(2,2,3,'Parent',figure(gcf));
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_3(X,6),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_3(X,6),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_3(X,6),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_3(X,6),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_3(X,6),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_3(X,6),'-m','LineWidth',2);
semilogy(X,sorted_MERGE_3(X,6),'-y','LineWidth',2);hold off
separation3 = num2str(near3 - 1);
distance = num2str(radius);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,.16]);
string1 = ['Number of top pairs identified by a method'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation3 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation3 ' intervening positions in the sequence '];
ylabel({string1;string2;string3});


subplot4 = subplot(2,2,4,'Parent',figure(gcf));
box(subplot4,'on');
grid(subplot4,'on');
hold(subplot4,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_4(X,6),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_4(X,6),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_4(X,6),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_4(X,6),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_4(X,6),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_4(X,6),'-m','LineWidth',2);
semilogy(X,sorted_MERGE_4(X,6),'-y','LineWidth',2);hold off
separation4 = num2str(near4 - 1);
distance = num2str(radius);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,.18]);
string1 = ['Number of top pairs identified by a method'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation4 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation4 ' intervening positions in the sequence '];
ylabel({string1;string2;string3});

%% 
saveas(gcf,figure_file,'fig');

%%
close all
%% Cumulative distances

sorted_md3_ZPX2_1(:,7) = cumsum(sorted_md3_ZPX2_1(:,4));
sorted_md4_ZPX2_1(:,7) = cumsum(sorted_md4_ZPX2_1(:,4));
sorted_slPSICOV_1(:,7) = cumsum(sorted_slPSICOV_1(:,4));
sorted_plmDCA_ZPX2_1(:,7) = cumsum(sorted_plmDCA_ZPX2_1(:,4));
sorted_GREMLIN_ZPX2_1(:,7) = cumsum(sorted_GREMLIN_ZPX2_1(:,4));
sorted_hpPCA_ZPX2_1(:,7) = cumsum(sorted_hpPCA_ZPX2_1(:,4));

sorted_md3_ZPX2_2(:,7) = cumsum(sorted_md3_ZPX2_2(:,4));
sorted_md4_ZPX2_2(:,7) = cumsum(sorted_md4_ZPX2_2(:,4));
sorted_slPSICOV_2(:,7) = cumsum(sorted_slPSICOV_2(:,4));
sorted_plmDCA_ZPX2_2(:,7) = cumsum(sorted_plmDCA_ZPX2_2(:,4));
sorted_GREMLIN_ZPX2_2(:,7) = cumsum(sorted_GREMLIN_ZPX2_2(:,4));
sorted_hpPCA_ZPX2_2(:,7) = cumsum(sorted_hpPCA_ZPX2_2(:,4));

sorted_md3_ZPX2_3(:,7) = cumsum(sorted_md3_ZPX2_3(:,4));
sorted_md4_ZPX2_3(:,7) = cumsum(sorted_md4_ZPX2_3(:,4));
sorted_slPSICOV_3(:,7) = cumsum(sorted_slPSICOV_3(:,4));
sorted_plmDCA_ZPX2_3(:,7) = cumsum(sorted_plmDCA_ZPX2_3(:,4));
sorted_GREMLIN_ZPX2_3(:,7) = cumsum(sorted_GREMLIN_ZPX2_3(:,4));
sorted_hpPCA_ZPX2_3(:,7) = cumsum(sorted_hpPCA_ZPX2_3(:,4));

sorted_md3_ZPX2_4(:,7) = cumsum(sorted_md3_ZPX2_4(:,4));
sorted_md4_ZPX2_4(:,7) = cumsum(sorted_md4_ZPX2_4(:,4));
sorted_slPSICOV_4(:,7) = cumsum(sorted_slPSICOV_4(:,4));
sorted_plmDCA_ZPX2_4(:,7) = cumsum(sorted_plmDCA_ZPX2_4(:,4));
sorted_GREMLIN_ZPX2_4(:,7) = cumsum(sorted_GREMLIN_ZPX2_4(:,4));
sorted_hpPCA_ZPX2_4(:,7) = cumsum(sorted_hpPCA_ZPX2_4(:,4));


%%
COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.8 0.8 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;

subplot1 = subplot(2,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

ncov = REF_length;
X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_1(X,7),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_1(X,7),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_1(X,7),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_1(X,7),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_1(X,7),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_1(X,7),'-m','LineWidth',2);
hold off
separation1 = num2str(near1 - 1);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,1600]);
    legend('3D\_MI','4D\_MI','PSICOV',...
        'plmDCA','GREMLIN','hpPCA',...
        'Location','NorthWest');
% 	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
string1 = ['Number of top pairs identified by a method within '];
string2 = 'the set of all the residues in the sequence ';
xlabel({string1;string2});
string1 = ['Cumulative sum of pairs distances'];
ylabel({string1});


subplot2 = subplot(2,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_2(X,7),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_2(X,7),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_2(X,7),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_2(X,7),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_2(X,7),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_2(X,7),'-m','LineWidth',2);
hold off
separation2 = num2str(near2 - 1);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,2200]);
string1 = ['Number of top pairs identified by a method within '];
string2 = 'the set of residues separated by at least ';
string3 = [separation2 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Cumulative sum of pairs distances '];
ylabel({string1});


subplot3 = subplot(2,2,3,'Parent',figure(gcf));
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_3(X,7),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_3(X,7),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_3(X,7),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_3(X,7),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_3(X,7),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_3(X,7),'-m','LineWidth',2);
separation3 = num2str(near3 - 1);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,2400]);
string1 = ['Number of top pairs identified by a method within '];
string2 = 'the set of residues separated by at least ';
string3 = [separation3 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Cumulative sum of pairs distances '];
ylabel({string1});


subplot4 = subplot(2,2,4,'Parent',figure(gcf));
box(subplot4,'on');
grid(subplot4,'on');
hold(subplot4,'all');

X = (1:ncov);
semilogy(X,sorted_md3_ZPX2_4(X,7),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_4(X,7),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_4(X,7),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_4(X,7),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_4(X,7),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_4(X,7),'-m','LineWidth',2);
separation4 = num2str(near4 - 1);
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,2500]);
string1 = ['Number of top pairs identified by a method within '];
string2 = 'the set of residues separated by at least ';
string3 = [separation4 ' intervening positions in the sequence '];
xlabel({string1;string2;string3});
string1 = ['Cumulative sum of pairs distances '];
ylabel({string1});

%% 
saveas(gcf,cumdist_figure_file,'fig');

%% Partition in <20 and >20 intervening positions

sorted_ZPX2_1(:,8) = abs(sorted_ZPX2_1(:,2)-sorted_ZPX2_1(:,3))<20;
sorted_ZPX2_1(:,9) = abs(sorted_ZPX2_1(:,2)-sorted_ZPX2_1(:,3))>=20;
sorted_ZPX2_1(:,10) = sorted_ZPX2_1(:,5).*sorted_ZPX2_1(:,8);
sorted_ZPX2_1(:,11) = sorted_ZPX2_1(:,5).*sorted_ZPX2_1(:,9);
sorted_ZPX2_1(:,12) = cumsum(sorted_ZPX2_1(:,10)/sum(sorted_ZPX2_1(:,5)));
sorted_ZPX2_1(:,13) = cumsum(sorted_ZPX2_1(:,11)/sum(sorted_ZPX2_1(:,5)));

sorted_md3_ZPX2_1(:,8) = abs(sorted_md3_ZPX2_1(:,2)-sorted_md3_ZPX2_1(:,3))<20;
sorted_md3_ZPX2_1(:,9) = abs(sorted_md3_ZPX2_1(:,2)-sorted_md3_ZPX2_1(:,3))>=20;
sorted_md3_ZPX2_1(:,10) = sorted_md3_ZPX2_1(:,5).*sorted_md3_ZPX2_1(:,8);
sorted_md3_ZPX2_1(:,11) = sorted_md3_ZPX2_1(:,5).*sorted_md3_ZPX2_1(:,9);
sorted_md3_ZPX2_1(:,12) = cumsum(sorted_md3_ZPX2_1(:,10)/sum(sorted_md3_ZPX2_1(:,5)));
sorted_md3_ZPX2_1(:,13) = cumsum(sorted_md3_ZPX2_1(:,11)/sum(sorted_md3_ZPX2_1(:,5)));

sorted_md4_ZPX2_1(:,8) = abs(sorted_md4_ZPX2_1(:,2)-sorted_md4_ZPX2_1(:,3))<20;
sorted_md4_ZPX2_1(:,9) = abs(sorted_md4_ZPX2_1(:,2)-sorted_md4_ZPX2_1(:,3))>=20;
sorted_md4_ZPX2_1(:,10) = sorted_md4_ZPX2_1(:,5).*sorted_md4_ZPX2_1(:,8);
sorted_md4_ZPX2_1(:,11) = sorted_md4_ZPX2_1(:,5).*sorted_md4_ZPX2_1(:,9);
sorted_md4_ZPX2_1(:,12) = cumsum(sorted_md4_ZPX2_1(:,10)/sum(sorted_md4_ZPX2_1(:,5)));
sorted_md4_ZPX2_1(:,13) = cumsum(sorted_md4_ZPX2_1(:,11)/sum(sorted_md4_ZPX2_1(:,5)));

sorted_slPSICOV_1(:,8) = abs(sorted_slPSICOV_1(:,2)-sorted_slPSICOV_1(:,3))<20;
sorted_slPSICOV_1(:,9) = abs(sorted_slPSICOV_1(:,2)-sorted_slPSICOV_1(:,3))>=20;
sorted_slPSICOV_1(:,10) = sorted_slPSICOV_1(:,5).*sorted_slPSICOV_1(:,8);
sorted_slPSICOV_1(:,11) = sorted_slPSICOV_1(:,5).*sorted_slPSICOV_1(:,9);
sorted_slPSICOV_1(:,12) = cumsum(sorted_slPSICOV_1(:,10)/sum(sorted_slPSICOV_1(:,5)));
sorted_slPSICOV_1(:,13) = cumsum(sorted_slPSICOV_1(:,11)/sum(sorted_slPSICOV_1(:,5)));

sorted_plmDCA_ZPX2_1(:,8) = abs(sorted_plmDCA_ZPX2_1(:,2)-sorted_plmDCA_ZPX2_1(:,3))<20;
sorted_plmDCA_ZPX2_1(:,9) = abs(sorted_plmDCA_ZPX2_1(:,2)-sorted_plmDCA_ZPX2_1(:,3))>=20;
sorted_plmDCA_ZPX2_1(:,10) = sorted_plmDCA_ZPX2_1(:,5).*sorted_plmDCA_ZPX2_1(:,8);
sorted_plmDCA_ZPX2_1(:,11) = sorted_plmDCA_ZPX2_1(:,5).*sorted_plmDCA_ZPX2_1(:,9);
sorted_plmDCA_ZPX2_1(:,12) = cumsum(sorted_plmDCA_ZPX2_1(:,10)/sum(sorted_plmDCA_ZPX2_1(:,5)));
sorted_plmDCA_ZPX2_1(:,13) = cumsum(sorted_plmDCA_ZPX2_1(:,11)/sum(sorted_plmDCA_ZPX2_1(:,5)));

sorted_GREMLIN_ZPX2_1(:,8) = abs(sorted_GREMLIN_ZPX2_1(:,2)-sorted_GREMLIN_ZPX2_1(:,3))<20;
sorted_GREMLIN_ZPX2_1(:,9) = abs(sorted_GREMLIN_ZPX2_1(:,2)-sorted_GREMLIN_ZPX2_1(:,3))>=20;
sorted_GREMLIN_ZPX2_1(:,10) = sorted_GREMLIN_ZPX2_1(:,5).*sorted_GREMLIN_ZPX2_1(:,8);
sorted_GREMLIN_ZPX2_1(:,11) = sorted_GREMLIN_ZPX2_1(:,5).*sorted_GREMLIN_ZPX2_1(:,9);
sorted_GREMLIN_ZPX2_1(:,12) = cumsum(sorted_GREMLIN_ZPX2_1(:,10)/sum(sorted_GREMLIN_ZPX2_1(:,5)));
sorted_GREMLIN_ZPX2_1(:,13) = cumsum(sorted_GREMLIN_ZPX2_1(:,11)/sum(sorted_GREMLIN_ZPX2_1(:,5)));

sorted_hpPCA_ZPX2_1(:,8) = abs(sorted_hpPCA_ZPX2_1(:,2)-sorted_hpPCA_ZPX2_1(:,3))<20;
sorted_hpPCA_ZPX2_1(:,9) = abs(sorted_hpPCA_ZPX2_1(:,2)-sorted_hpPCA_ZPX2_1(:,3))>=20;
sorted_hpPCA_ZPX2_1(:,10) = sorted_hpPCA_ZPX2_1(:,5).*sorted_hpPCA_ZPX2_1(:,8);
sorted_hpPCA_ZPX2_1(:,11) = sorted_hpPCA_ZPX2_1(:,5).*sorted_hpPCA_ZPX2_1(:,9);
sorted_hpPCA_ZPX2_1(:,12) = cumsum(sorted_hpPCA_ZPX2_1(:,10)/sum(sorted_hpPCA_ZPX2_1(:,5)));
sorted_hpPCA_ZPX2_1(:,13) = cumsum(sorted_hpPCA_ZPX2_1(:,11)/sum(sorted_hpPCA_ZPX2_1(:,5)));

X01 = find(sorted_ZPX2_1(X,8));
X11 = find(sorted_md3_ZPX2_1(X,8));
X21 = find(sorted_md4_ZPX2_1(X,8));
X31 = find(sorted_slPSICOV_1(X,8));
X41 = find(sorted_plmDCA_ZPX2_1(X,8));
X51 = find(sorted_GREMLIN_ZPX2_1(X,8));
X61 = find(sorted_hpPCA_ZPX2_1(X,8));

sX01 = 1:size(X01,1);
sX11 = 1:size(X11,1);
sX21 = 1:size(X21,1);
sX31 = 1:size(X31,1);
sX41 = 1:size(X41,1);
sX51 = 1:size(X51,1);
sX61 = 1:size(X61,1);

X02 = find(sorted_ZPX2_1(X,9));
X12 = find(sorted_md3_ZPX2_1(X,9));
X22 = find(sorted_md4_ZPX2_1(X,9));
X32 = find(sorted_slPSICOV_1(X,9));
X42 = find(sorted_plmDCA_ZPX2_1(X,9));
X52 = find(sorted_GREMLIN_ZPX2_1(X,9));
X62 = find(sorted_hpPCA_ZPX2_1(X,9));

sX02 = 1:size(X02,1);
sX12 = 1:size(X12,1);
sX22 = 1:size(X22,1);
sX32 = 1:size(X32,1);
sX42 = 1:size(X42,1);
sX52 = 1:size(X52,1);
sX62 = 1:size(X62,1);


%%
COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.9 0.4 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;

subplot1 = subplot(1,3,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

ncov = REF_length;
X = (1:ncov);
semilogy(X,sorted_ZPX2_1(X,6),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,sorted_md3_ZPX2_1(X,6),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_1(X,6),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_1(X,6),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_1(X,6),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_1(X,6),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_1(X,6),'-m','LineWidth',2);
hold off
separation1 = num2str(near1 - 1);
distance = num2str(radius);
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of all the residues in the sequence ';
string3 = ' ';
set(gca,'Ylim',[-0.00 0.127]);
    legend('2D\_MI','3D\_MI','4D\_MI','PSICOV',...
        'plmDCA','GREMLIN','hpPCA',...
        'Location','NorthWest');
% 	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
xlabel('Number of top scores among all pairs ');
ylabel({string1});

    
subplot2 = subplot(1,3,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

ncov = REF_length;
X = (1:ncov);

semilogy(sX01,sorted_ZPX2_1(X01,12),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(sX11,sorted_md3_ZPX2_1(X11,12),'-r','LineWidth',2);
semilogy(sX21,sorted_md4_ZPX2_1(X21,12),'-k','LineWidth',2);
semilogy(sX31,sorted_slPSICOV_1(X31,12),'-b','LineWidth',2);
semilogy(sX41,sorted_plmDCA_ZPX2_1(X41,12),'-g','LineWidth',2);
semilogy(sX51,sorted_GREMLIN_ZPX2_1(X51,12),'-c','LineWidth',2);
semilogy(sX61,sorted_hpPCA_ZPX2_1(X61,12),'-m','LineWidth',2);
hold off
separation1 = num2str(near1 - 1);
distance = num2str(radius);
set(gca,'Ylim',[-0.00 0.107]);
string1 = ['Number of top scores among the pairs separated by  '];
string2 = ['less than ' separation2 ' intervening positions in sequence '];
string3 = 'within the top L scores ';
xlabel({string1;string2});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by less than ';
string3 = [separation2 ' intervening positions in the sequence '];
ylabel({string1});


subplot3 = subplot(1,3,3,'Parent',figure(gcf));
box(subplot3,'on');
grid(subplot3,'on');
hold(subplot3,'all');

X = (1:ncov);

semilogy(sX02,sorted_ZPX2_1(X02,13),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(sX12,sorted_md3_ZPX2_1(X12,13),'-r','LineWidth',2);
semilogy(sX22,sorted_md4_ZPX2_1(X22,13),'-k','LineWidth',2);
semilogy(sX32,sorted_slPSICOV_1(X32,13),'-b','LineWidth',2);
semilogy(sX42,sorted_plmDCA_ZPX2_1(X42,13),'-g','LineWidth',2);
semilogy(sX52,sorted_GREMLIN_ZPX2_1(X52,13),'-c','LineWidth',2);
semilogy(sX62,sorted_hpPCA_ZPX2_1(X62,13),'-m','LineWidth',2);
separation2 = num2str(near4 - 1);
distance = num2str(radius);
set(gca,'Ylim',[-0.00 0.027]);
xlabel('Number of pairs  ');
string1 = ['Number of top scores among the pairs separated by  '];
string2 = ['at least ' separation2 ' intervening positions in sequence '];
string3 = 'within the top L scores ';
xlabel({string1;string2});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation2 ' intervening positions in the sequence '];
ylabel({string1});

%%
COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.9 0.4 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;

subplot1 = subplot(1,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

ncov = REF_length;
X = (1:ncov);
semilogy(X,sorted_ZPX2_1(X,6),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,sorted_md3_ZPX2_1(X,6),'-r','LineWidth',2);
semilogy(X,sorted_md4_ZPX2_1(X,6),'-k','LineWidth',2);
semilogy(X,sorted_slPSICOV_1(X,6),'-b','LineWidth',2);
semilogy(X,sorted_plmDCA_ZPX2_1(X,6),'-g','LineWidth',2);
semilogy(X,sorted_GREMLIN_ZPX2_1(X,6),'-c','LineWidth',2);
semilogy(X,sorted_hpPCA_ZPX2_1(X,6),'-m','LineWidth',2);
hold off
separation1 = num2str(near1 - 1);
distance = num2str(radius);
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of all the residues in the sequence ';
string3 = ' ';
set(gca,'Ylim',[-0.00 0.127]);
    legend('2D\_MI','3D\_MI','4D\_MI','PSICOV',...
        'plmDCA','GREMLIN','hpPCA',...
        'Location','NorthWest');
xlabel('Number of top scores among all pairs ');
ylabel({string1});

    
subplot2 = subplot(1,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

X = (1:ncov);

semilogy(sX02,sorted_ZPX2_1(X02,13),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(sX12,sorted_md3_ZPX2_1(X12,13),'-r','LineWidth',2);
semilogy(sX22,sorted_md4_ZPX2_1(X22,13),'-k','LineWidth',2);
semilogy(sX32,sorted_slPSICOV_1(X32,13),'-b','LineWidth',2);
semilogy(sX42,sorted_plmDCA_ZPX2_1(X42,13),'-g','LineWidth',2);
semilogy(sX52,sorted_GREMLIN_ZPX2_1(X52,13),'-c','LineWidth',2);
semilogy(sX62,sorted_hpPCA_ZPX2_1(X62,13),'-m','LineWidth',2);
separation2 = num2str(near4 - 1);
distance = num2str(radius);
set(gca,'Ylim',[-0.00 0.027]);
xlabel('Number of pairs  ');
string1 = ['Number of top scores among the pairs separated by  '];
string2 = ['at least ' separation2 ' intervening positions in sequence '];
string3 = 'within the top L scores ';
xlabel({string1;string2});
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = 'within the set of residues separated by at least ';
string3 = [separation2 ' intervening positions in the sequence '];
ylabel({string1});

%% 
saveas(gcf,partdist_figure_file,'fig');

%%
save(matfile);

%%
close all
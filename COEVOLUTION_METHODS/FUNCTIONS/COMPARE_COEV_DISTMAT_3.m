% This script shows an example of how to compare the pairs identified by a
% coevolution methods with the pairs that are found within a
% threshold distance in the corresponding X-ray structure. It differs from
% 'COMPARE_COEV_DISTMAT' as here we are looking at the distance matrix
% directly and we determine what percentage of the pairs separated by 8 
% angstrom or less is present in the top covarions identified by a method.
% It differs further from 'COMPARE_COEV_DISTMAT_2' as it weighs the
% distances in column 7 of the sorted matrices such that the cumulative sum
% in column 8 reflects also the closeness of the components in a pair.
%%
% Read in the msa

[sMDH_95_smsa,sMDH_95_nmsa] = faln_to_nmsa('sMDH_95.faln');

% [sMDH_95_smsa,sMDH_95_nmsa] = aln_to_nmsa('sMDH_95.aln');

% IMPORTANT! Set the range of each sequence that will be used to calculate the 
% coevolution matrices: make sure REF_length is consistent with the actual
% range of residues for which there are coordinates in the pdb file.
% The pdb file must have a header with the correct sequence for this script
% to work. This is a little confusing because different depositors follow 
% different rules to fill the SEQRES field of the header. We suggest that 
% the 'SEQRES' value should be equal to the residue number of the last aa 
% for which there are coordinates. It doesn't matter if there  
% are initial or central gaps in the actual coordinates, because these will 
% correspond to all zero values in the distance matrix. Therefore if you
% have leading or trailing residues (e.g. a HIS-tag) that were not included
% in the msa and therefore in the coevolution matrix, it is crucial to
% specify correctly what are the 1st and last residue numbers in the
% coordinates, and these must correspond to the msa for which there is a 
% coevolution matrix. If the coordinates start at a residue that comes after 
% the beginning of the msa it is better to trim the msa in such a way that
% there is full correspondence between 1st and last residues in the msa and
% in the pdb.

REF_length = numel(sMDH_95_nmsa(1,1:353));

% or just

% REF_length = 235;

% The following is the sequence in the msa corresponding to the X-ray 
% structure. It is not used, but it is convenient to have for checking
% various things.

ref_seq_no = 391;

nmsa = sMDH_95_nmsa(:,1:REF_length);
smsa = int2aa(nmsa);
ncols = size(nmsa,2);
pdbfile = '1HUV_mod.pdb';

%%

[MI,MI_covar_vec] = get_nmsa_covar_vec(nmsa,30,'MI');
[logR,logR_covar_vec] = get_nmsa_covar_vec(nmsa,30,'logR');
[ZRES,ZRES_covar_vec] = get_nmsa_covar_vec(nmsa,30,'ZRES2');
[SCA,SCA_covar_vec] = get_nmsa_covar_vec(nmsa,30,'rsemSCA');
[MIP,MIP_covar_vec] = get_nmsa_covar_vec(nmsa,30,'MIP');
[ZPX2,ZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'ZPX2');
[nbZPX2,nbZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'nbZPX2');
% For higher contrast and all positive values use 'nbsZPX2'.
% [nbZPX2,nbZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'nbsZPX2');
[gbZPX2,gbZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'gbZPX2');
[fgbZPX2,fgbZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'fgbZPX2');
[dbZPX2,dbZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'dbZPX2');
[dgbZPX2,dgbZPX2_covar_vec] = get_nmsa_covar_vec(nmsa,30,'dgbZPX2');
[DCA,DCA_covar_vec] = get_nmsa_covar_vec(nmsa,30,'DCA');
[OMES,OMES_covar_vec] = get_nmsa_covar_vec(nmsa,30,'OMES');
[McBASC,McBASC_covar_vec] = get_nmsa_covar_vec(nmsa,30,'McBASC');
[ELSC,ELSC_covar_vec] = get_nmsa_covar_vec(nmsa,30,'ELSC');
[fodorSCA,fodorSCA_covar_vec] = get_nmsa_covar_vec(nmsa,30,'fodorSCA');

%%
% Here we calculate a matrix of weights to correct for the presence of
% gaps. ramaSCA, rsemSCA, and logR do not need a correction. Not clear
% whether OMES, McBASC, ELSC, and fodorSCA need a correction; probably they
% do.

    gapW0 = ones(ncols,ncols);
    gapW1 = correct_coevmat_forgaps(nmsa);
    gapW2 = gapW1.^2;
    gapW3 = gapW1.^3;
    
 MI_orig_1 = MI;
 logR_orig_1 = logR;
 MIP_orig_1 = MIP;
 ZRES_orig_1 = ZRES;
 ZPX2_orig_1 = ZPX2;
 nbZPX2_orig_1 = nbZPX2;
 gbZPX2_orig_1 = gbZPX2;
 dbZPX2_orig_1 = dbZPX2;
 fgbZPX2_orig_1 = fgbZPX2;
 dgbZPX2_orig_1 = dgbZPX2;
 DCA_orig_1 = DCA;
 OMES_orig_1 = OMES;
 McBASC_orig_1 = McBASC;
 ELSC_orig_1 = ELSC;
 SCA_orig_1 = SCA;
 fodorSCA_orig_1 = fodorSCA;
 
%%
% Here we make sure all the matrices are positive to prevent odd phenomena
% in the sorting.

 MI_orig_2 = MI_orig_1 - min(MI_orig_1(:));
 logR_orig_2 = logR_orig_1 - min(logR_orig_1(:));
 MIP_orig_2 = MIP_orig_1 - min(MIP_orig_1(:));
 ZRES_orig_2 = ZRES_orig_1 - min(ZRES_orig_1(:));
 ZPX2_orig_2 = ZPX2_orig_1 -min(ZPX2_orig_1(:));
 nbZPX2_orig_2 = nbZPX2_orig_1 - min(nbZPX2_orig_1(:));
 gbZPX2_orig_2 = gbZPX2_orig_1 - min(gbZPX2_orig_1(:));
 dbZPX2_orig_2 = dbZPX2_orig_1 - min(dbZPX2_orig_1(:));
 fgbZPX2_orig_2 = fgbZPX2_orig_1 - min(fgbZPX2_orig_1(:));
 dgbZPX2_orig_2 = dgbZPX2_orig_1 - min(dgbZPX2_orig_1(:));
 DCA_orig_2 = DCA_orig_1 - min(DCA_orig_1(:));
 OMES_orig_2 = OMES_orig_1 - min(OMES_orig_1(:));
 McBASC_orig_2 = McBASC_orig_1 - min(McBASC_orig_1(:));
 ELSC_orig_2 = ELSC_orig_1 - min(ELSC_orig_1(:));
 SCA_orig_2 = SCA_orig_1 - min(SCA_orig_1(:));
 fodorSCA_orig_2 = fodorSCA_orig_1 - min(fodorSCA_orig_1(:));
 
%%
% Here we corrects the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW1,gapW2,gapW3) can be applied. Use gapW0 if you don't 
% want any correction. logR, and SCA do not need a correction.

 gW = gapW2;

 MI = MI_orig_2.*gW;
%  logR = logR_orig_2.*gW;
 logR = logR_orig_2;
 MIP = MIP_orig_2.*gW;
 ZRES = ZRES_orig_2.*gW;
 ZPX2 = ZPX2_orig_2.*gW;
 nbZPX2 = nbZPX2_orig_2.*gW;
 gbZPX2 = gbZPX2_orig_2.*gW;
 dbZPX2 = dbZPX2_orig_2.*gW;
 fgbZPX2 = fgbZPX2_orig_2.*gW;
 dgbZPX2 = dgbZPX2_orig_2.*gW;
 DCA = DCA_orig_2.*gW;
%  OMES = OMES_orig_2.*gW;
%  McBASC = McBASC_orig_2.*gW;
%  ELSC = ELSC_orig_2.*gW;
 OMES = OMES_orig_2;
 McBASC = McBASC_orig_2;
 ELSC = ELSC_orig_2; 
% SCA = SCA_orig_2.*gW;
 SCA = SCA_orig_2;
 fodorSCA = fodorSCA_orig_2.*gW;

%%
near = 1;
near1 = near;
ncov = round(REF_length/2);
radius = 8;
%---------------------1st set----------------------------------------------
% Usage: [c_distances,sorted_c_distances,sorted_mat1,...
%           sorted_mat2,sorted_mat3,sorted_mat4,pdbstruct] = ...
%           coev_distance_matrix_4(pdbfile,chain,first_res_no,last_res_no,...
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

[c_distances,sorted_c_distances,sorted_ZPX2_1,...
    sorted_dbZPX2_1,sorted_fgbZPX2_1,sorted_DCA_1,pdbstruct_1] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,ZPX2,dbZPX2,...
    fgbZPX2,DCA,radius,near,ncov,0);

%---------------------2nd set----------------------------------------------

[~,~,sorted_MIP_1,sorted_nbZPX2_1,sorted_gbZPX2_1,sorted_dgbZPX2_1] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,MIP,nbZPX2,...
    gbZPX2,dgbZPX2,radius,near,ncov,0);

%---------------------3rd set----------------------------------------------

[~,~,sorted_MI_1,sorted_logR_1,sorted_ZRES_1,sorted_SCA_1] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,MI,logR,...
    ZRES,SCA,radius,near,ncov,0);

%---------------------4th set----------------------------------------------

[~,~,sorted_OMES_1,sorted_McBASC_1,sorted_ELSC_1,sorted_fodorSCA_1] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,OMES,McBASC,...
    ELSC,fodorSCA,radius,near,ncov,0);

%%
near = 4;
near2 = near;
ncov = round(REF_length/2);
radius = 8;

%---------------------1st set----------------------------------------------

[~,~,sorted_ZPX2_2,sorted_dbZPX2_2,sorted_fgbZPX2_2,...
    sorted_DCA_2,pdbstruct_2] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,ZPX2,dbZPX2,...
    fgbZPX2,DCA,radius,near,ncov,0);

%---------------------2nd set----------------------------------------------

[~,~,sorted_MIP_2,sorted_nbZPX2_2,sorted_gbZPX2_2,sorted_dgbZPX2_2] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,MIP,nbZPX2,...
    gbZPX2,dgbZPX2,radius,near,ncov,0);

%---------------------3rd set----------------------------------------------

[~,~,sorted_MI_2,sorted_logR_2,sorted_ZRES_2,sorted_SCA_2] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,MI,logR,...
    ZRES,SCA,radius,near,ncov,0);

%---------------------4th set----------------------------------------------

[~,~,sorted_OMES_2,sorted_McBASC_2,sorted_ELSC_2,sorted_fodorSCA_2] = ...
    coev_distance_matrix_4(pdbfile,1,4,356,OMES,McBASC,...
    ELSC,fodorSCA,radius,near,ncov,0);


%%
COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.8 0.4 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;

subplot1 = subplot(1,2,1,'Parent',figure(gcf));
box(subplot1,'on');
grid(subplot1,'on');
hold(subplot1,'all');

X = (1:ncov);
semilogy(X,sorted_MI_1(X,8),'-r','LineWidth',2);
semilogy(X,sorted_logR_1(X,8),'-','Color',[0.3,0.6,0.2],'LineWidth',2);
% semilogy(X,MIP_1(X,8),'-','Color',[0.3,0.6,0.2],'LineWidth',2);
semilogy(X,sorted_ZPX2_1(X,8),'-m','LineWidth',2);
semilogy(X,sorted_DCA_1(X,8),'-g','LineWidth',2);
semilogy(X,sorted_nbZPX2_1(X,8),'-c','LineWidth',2);
semilogy(X,sorted_dbZPX2_1(X,8),'-y','LineWidth',2);
% semilogy(X,sorted_fgbZPX2_1(X,8),'-k','LineWidth',2);
semilogy(X,sorted_dgbZPX2_1(X,8),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,sorted_SCA_1(X,8),'--','Color','y','LineWidth',2);
semilogy(X,sorted_OMES_1(X,8),'--','Color','g','LineWidth',2);
semilogy(X,sorted_McBASC_1(X,8),'--','Color','m','LineWidth',2);
semilogy(X,sorted_ELSC_1(X,8),'--','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,sorted_fodorSCA_1(X,8),'--','Color','c','LineWidth',2);
hold off
separation1 = num2str(near1 - 1);
distance = num2str(radius);
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = ['within the set of residues separated by at least '];
string3 = [separation1 ' intervening positions in the sequence '];
% string2 = ['within the set of consecutive residues in the sequence '];
% string3 = [' '];
set(gca,'Xlim',[0,max(X)],'Ylim',[0.001,.06]);
% set(gca,'Xlim',[0,max(X)],'Ylim',[0,1.0E3]);
% set(gca,'Xlim',[0,max(X)],'Ylim',[0,2.0E3]);
legend('MI','logR','ZPX2','DCA','nbZPX2',...
    'dbZPX2','dgbZPX2','SCA',...
    'OMES','McBASC','ELSC','fodorSCA','Location','NorthWest');
% 	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
	xlabel('Number of top pairs identified by a method');
	ylabel({string1;string2;string3});


subplot2 = subplot(1,2,2,'Parent',figure(gcf));
box(subplot2,'on');
grid(subplot2,'on');
hold(subplot2,'all');

X = (1:ncov);
semilogy(X,sorted_MI_2(X,8),'-r','LineWidth',2);
semilogy(X,sorted_logR_2(X,8),'-','Color',[0.3,0.6,0.2],'LineWidth',2);
% semilogy(X,MIP_2(X,8),'-','Color',[0.3,0.6,0.2],'LineWidth',2);
semilogy(X,sorted_ZPX2_2(X,8),'-m','LineWidth',2);
semilogy(X,sorted_DCA_2(X,8),'-g','LineWidth',2);
semilogy(X,sorted_nbZPX2_2(X,8),'-c','LineWidth',2);
semilogy(X,sorted_dbZPX2_2(X,8),'-y','LineWidth',2);
% semilogy(X,sorted_fgbZPX2_2(X,8),'-k','LineWidth',2);
semilogy(X,sorted_dgbZPX2_2(X,8),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,sorted_SCA_2(X,8),'--','Color','y','LineWidth',2);
semilogy(X,sorted_OMES_2(X,8),'--','Color','g','LineWidth',2);
semilogy(X,sorted_McBASC_2(X,8),'--','Color','m','LineWidth',2);
semilogy(X,sorted_ELSC_2(X,8),'--','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,sorted_fodorSCA_2(X,8),'--','Color','c','LineWidth',2);
hold off
separation2 = num2str(near2 - 1);
distance = num2str(radius);
string1 = ['Percentage of all pairs separated by less than ' ...
    distance ' angstroms'];
string2 = ['within the set of residues separated by at least '];
string3 = [separation2 ' intervening positions in the sequence '];
set(gca,'Xlim',[0,max(X)],'Ylim',[-0.001,.08]);
% set(gca,'Xlim',[0,max(X)],'Ylim',[0,1.0E3]);
% set(gca,'Xlim',[0,max(X)],'Ylim',[0,2.0E3]);
% legend('Ideal','MI','logR','ZPX2','DCA','nbZPX2',...
%     'dbZPX2','dgbZPX2','SCA',...
%     'OMES','McBASC','ELSC','fodorSCA','Location','NorthWest');
% 	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
	xlabel('Number of top pairs identified by a method');
	ylabel({string1;string2;string3});


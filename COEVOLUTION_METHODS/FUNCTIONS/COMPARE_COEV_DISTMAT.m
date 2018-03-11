% This script shows an example of how to compare the pairs identified by a
% coevolution methods with the pairs that are found within a
% threshold distance in the corresponding X-ray structure.
%%
% Read in the msa

% [ATP12_95_smsa,ATP12_95_nmsa] = faln_to_nmsa('atp12_95.faln');

[ATP12_95_smsa,ATP12_95_nmsa] = aln_to_nmsa('atp12_95.aln');

% IMPORTANT! Set the range of each sequence that will be used to calculate the 
% coevolution matrices: make sure REF_length is consistent with the 'SEQRES' 
% field in the pdb file by modifying accordingly the following line.

REF_length = numel(ATP12_95_nmsa(1,1:235));

% or just

% REF_length = 235;

%%
nmsa = ATP12_95_nmsa(:,1:REF_length);
smsa = int2aa(nmsa);
ncols = size(nmsa,2);
pdbfile = '2R31_renumbered.pdb';

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
 RSEM_SCA_orig_1 = SCA;
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
 RSEM_SCA_orig_2 = RSEM_SCA_orig_1 - min(RSEM_SCA_orig_1(:));
 fodorSCA_orig_2 = fodorSCA_orig_1 - min(fodorSCA_orig_1(:));
 
%%
% Here we corrects the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW1,gapW2,gapW3) can be applied. Use gapW0 if you don't 
% want any correction. logR, and SCA do not need a correction.

 gW = gapW0;

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
 OMES = OMES_orig_2.*gW;
 McBASC = McBASC_orig_2.*gW;
 ELSC = ELSC_orig_2.*gW;
% RSEM_SCA = fodorSCA_orig_2.*gW;
 RSEM_SCA = fodorSCA_orig_2;
 fodorSCA = fodorSCA_orig_2.*gW;

%%
near = 1;
ncov = round(REF_length/2);
radius = 8;
%---------------------1st set----------------------------------------------
% The last input for 'coev_distance_matrix' can be:
% 0 = no distance matrix output
% 1 = distance matrix output with yellow background
% 2 = distance matrix output with blue background and colors representing
% the number of atomic interactions
% >2 = distance matrix output with blue background and colors representing
% the number of atomic interactions, with overlaid covarions identified by
% a method.

[score1,entropy1,relentropy1,interactions,sfu_interactions_1,...
    sfu_interactions_2,sfu_interactions_3,...
    sfu_interactions_4,PDB_structure] = ...
    coev_distance_matrix(nmsa,smsa,pdbfile,1,...
    ZPX2,dbZPX2,fgbZPX2,DCA,radius,near,ncov,2);

[ZPX2_s_coev_mat,ZPX2_s_inter_mat] = ...
    cum_dist_coev(ZPX2,sfu_interactions_1,near);
[dbZPX2_s_coev_mat,dbZPX2_s_inter_mat] = ...
    cum_dist_coev(dbZPX2,sfu_interactions_2,near);
[fgbZPX2_s_coev_mat,fgbZPX2_s_inter_mat] = ...
    cum_dist_coev(fgbZPX2,sfu_interactions_3,near);
[DCA_s_coev_mat,DCA_s_inter_mat] = ...
    cum_dist_coev(DCA,sfu_interactions_4,near);

%---------------------2nd set----------------------------------------------
[score2,entropy2,relentropy2,~,sfu_interactions_1,...
    sfu_interactions_2,sfu_interactions_3,sfu_interactions_4] = ...
    coev_distance_matrix(nmsa,smsa,pdbfile,1,...
    MIP,nbZPX2,gbZPX2,dgbZPX2,radius,near,ncov,0);

[MIP_s_coev_mat,MIP_s_inter_mat] = ...
    cum_dist_coev(MIP,sfu_interactions_1,near);
[nbZPX2_s_coev_mat,nbZPX2_s_inter_mat] = ...
    cum_dist_coev(nbZPX2,sfu_interactions_2,near);
[gbZPX2_s_coev_mat,gbZPX2_s_inter_mat] = ...
    cum_dist_coev(gbZPX2,sfu_interactions_3,near);
[dgbZPX2_s_coev_mat,dgbZPX2_s_inter_mat] = ...
    cum_dist_coev(dgbZPX2,sfu_interactions_4,near);

%---------------------3rd set----------------------------------------------
[score3,entropy3,relentropy3,~,sfu_interactions_1,...
    sfu_interactions_2,sfu_interactions_3,sfu_interactions_4] = ...
    coev_distance_matrix(nmsa,smsa,pdbfile,1,...
    MI,logR,ZRES,SCA,radius,near,ncov,0);

[MI_s_coev_mat,MI_s_inter_mat] = ...
    cum_dist_coev(MI,sfu_interactions_1,near);
[logR_s_coev_mat,logR_s_inter_mat] = ...
    cum_dist_coev(logR,sfu_interactions_2,near);
[ZRES_s_coev_mat,ZRES_s_inter_mat] = ...
    cum_dist_coev(ZRES,sfu_interactions_3,near);
[SCA_s_coev_mat,SCA_s_inter_mat] = ...
    cum_dist_coev(SCA,sfu_interactions_4,near);

%---------------------4th set----------------------------------------------
[score4,entropy4,relentropy4,~,sfu_interactions_1,...
    sfu_interactions_2,sfu_interactions_3,sfu_interactions_4] = ...
    coev_distance_matrix(nmsa,smsa,pdbfile,1,...
    OMES,McBASC,ELSC,fodorSCA,radius,near,ncov,0);

[OMES_s_coev_mat,OMES_s_inter_mat] = ...
    cum_dist_coev(OMES,sfu_interactions_1,near);
[McBASC_s_coev_mat,McBASC_s_inter_mat] = ...
    cum_dist_coev(McBASC,sfu_interactions_2,near);
[ELSC_s_coev_mat,ELSC_s_inter_mat] = ...
    cum_dist_coev(ELSC,sfu_interactions_3,near);
[fodorSCA_s_coev_mat,fodorSCA_s_inter_mat] = ...
    cum_dist_coev(fodorSCA,sfu_interactions_4,near);

%--------------------------------------------------------------------------

s_inter_mat = sort_matrix_descend_2(full(interactions),near);
s_inter_mat(:,4) = 0;
s_inter_mat(:,5) = cumsum(s_inter_mat(:,1));

%%
COEV_DIST_MAT = figure; 
    	set(COEV_DIST_MAT,'Units','normalized','Position',[0 0.2 0.5 0.5 ],...
    	'Name','COEV versus DISTANCE MATRIX'); clf;
box('on');
grid('on');
hold('all');

X = (1:ncov);
semilogy(s_inter_mat(X,5),'-b','LineWidth',2);
hold on
semilogy(X,MI_s_coev_mat(X,5),'-r','LineWidth',2);
semilogy(X,logR_s_coev_mat(X,5),'-','Color',[0.3,0.6,0.2],'LineWidth',2);
% semilogy(X,MIP_s_coev_mat(X,5),'-','Color',[0.3,0.6,0.2],'LineWidth',2);
semilogy(X,ZPX2_s_coev_mat(X,5),'-m','LineWidth',2);
semilogy(X,DCA_s_coev_mat(X,5),'-g','LineWidth',2);
semilogy(X,nbZPX2_s_coev_mat(X,5),'-c','LineWidth',2);
semilogy(X,dbZPX2_s_coev_mat(X,5),'-y','LineWidth',2);
% semilogy(X,fgbZPX2_s_coev_mat(X,5),'-k','LineWidth',2);
% semilogy(X,gbZPX2_s_coev_mat(X,5),'-','Color',[1 0.5 0],'LineWidth',2);
semilogy(X,dgbZPX2_s_coev_mat(X,5),'-','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,SCA_s_coev_mat(X,5),'--','Color','y','LineWidth',2);
semilogy(X,OMES_s_coev_mat(X,5),'--','Color','g','LineWidth',2);
semilogy(X,McBASC_s_coev_mat(X,5),'--','Color','m','LineWidth',2);
semilogy(X,ELSC_s_coev_mat(X,5),'--','Color',[0,0.5,1],'LineWidth',2);
semilogy(X,fodorSCA_s_coev_mat(X,5),'--','Color','c','LineWidth',2);
hold off
separation = num2str(near - 1);
distance = num2str(radius);
string1 = ['Cumulative count of atomic interactions within ' ...
    distance ' angstroms'];
string2 = ['between residues separated by ' ...
    separation ' positions in the sequence '];
set(gca,'Xlim',[0,max(X)],'Ylim',[0,2E3]);
% set(gca,'Xlim',[0,max(X)],'Ylim',[0,1.0E3]);
% set(gca,'Xlim',[0,max(X)],'Ylim',[0,2.0E3]);
legend('Ideal','MI','logR','ZPX2','DCA','nbZPX2',...
    'dbZPX2','dgbZPX2','SCA',...
    'OMES','McBASC','ELSC','fodorSCA','Location','NorthWest');
	title('COEV DISTANCE MATRIX COMPARISON','FontSize',14,'FontWeight','n');
	xlabel('Number of top pairs identified by a method');
	ylabel({string1;string2});


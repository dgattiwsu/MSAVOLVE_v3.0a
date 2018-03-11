% This script can be called after loading a simulation run of MSAvolve. 
% Here we  compares the performance of several MI methods on simulated 
% data sets. In this case we use directly the COV or cov_COV matrix rather  
% than the original vector of covarions 'covar_vec'. Results from all the 
% runs are averaged. This script differs from COMPARE_COEV_METHODS_2 as it
% can implement a correction of the coevolution matrices based on the presence
% of gaps. Three progressively stronger corrections (gapW,gapW2,gapW3) can 
% be applied. logR, and SCA do not need a correction.

%%
[nrows,ncols] = size(MSA_select_ALL(:,:,1));

gW = zeros(ncols,ncols,end_cycle); 
gapW0 = ones(ncols,ncols);

for n = 1:end_cycle
    gapW = correct_coevmat_forgaps(MSA_select_ALL(:,:,n));
    gapW2 = gapW.^2;
    gapW3 = gapW.^3;
    gW(:,:,n) = gapW3;
end

%%
totpairs = npos*npos;
cum_covarions = zeros(totpairs,end_cycle);
cum_MI = zeros(totpairs,end_cycle);
cum_MIP = zeros(totpairs,end_cycle);
cum_logR = zeros(totpairs,end_cycle);
cum_DCA = zeros(totpairs,end_cycle);
cum_ZPX2 = zeros(totpairs,end_cycle);
cum_nbZPX2 = zeros(totpairs,end_cycle);
cum_dbZPX2 = zeros(totpairs,end_cycle);
cum_gbZPX2 = zeros(totpairs,end_cycle);
cum_fgbZPX2 = zeros(totpairs,end_cycle);
cum_dgbZPX2 = zeros(totpairs,end_cycle);
cum_RSEM_SCA = zeros(totpairs,end_cycle);
% cum_RAMA_SCA = zeros(totpairs,end_cycle);
cum_fodorSCA = zeros(totpairs,end_cycle);
cum_OMES = zeros(totpairs,end_cycle);
cum_McBASC = zeros(totpairs,end_cycle);
cum_ELSC = zeros(totpairs,end_cycle);


for n = 1:end_cycle
% We can run the comparison on the total COV matrix ...
% covarions = sort_matrix_descend_2(COV_ALL(:,:,n),1);
% ... or just the true covarions matrix.
covarions = sort_matrix_descend_2(cov_COV_ALL(:,:,n),1);

cum_covarions(:,n) = cumsum(covarions(:,1));

% Run MI
MI = MI_ALL(:,:,n);
MI = MI - min(MI(:));
% MI = MI.*gW(:,:,n);
sorted_MI = sort_matrix_descend_2(MI,1);
sorted_MI(:,4) = 0;
[~,ia,ib] = intersect(sorted_MI(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_MI(ia(i),4) = covarions(ib(i),1);
end
cum_MI(:,n) = cumsum(sorted_MI(:,4));


% Run MIP
MIP = MIP_ALL(:,:,n);
MIP = MIP - min(MIP(:));
% MIP = MIP.*gW(:,:,n);
sorted_MIP = sort_matrix_descend_2(MIP,1);
sorted_MIP(:,4) = 0;
[~,ia,ib] = intersect(sorted_MIP(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_MIP(ia(i),4) = covarions(ib(i),1);
end
cum_MIP(:,n) = cumsum(sorted_MIP(:,4));


% Run logR
logR = logR_ALL(:,:,n);
logR = logR - min(logR(:));
% logR = logR.*gW(:,:,n);
sorted_logR = sort_matrix_descend_2(logR,1);
sorted_logR(:,4) = 0;
[~,ia,ib] = intersect(sorted_logR(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_logR(ia(i),4) = covarions(ib(i),1);
end
cum_logR(:,n) = cumsum(sorted_logR(:,4));


% Run ZPX2
ZPX2 = ZPX2_ALL(:,:,n);
ZPX2 = ZPX2 - min(ZPX2(:));
% ZPX2 = ZPX2.*gW(:,:,n);
sorted_ZPX2 = sort_matrix_descend_2(ZPX2,1);sorted_ZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_ZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_ZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_ZPX2(:,n) = cumsum(sorted_ZPX2(:,4));


% Run DCA
DCA = DCA_ALL(:,:,n);
DCA = DCA - min(DCA(:));
% DCA = DCA.*gW(:,:,n);
sorted_DCA = sort_matrix_descend_2(DCA,1);
sorted_DCA(:,4) = 0;
[~,ia,ib] = intersect(sorted_DCA(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_DCA(ia(i),4) = covarions(ib(i),1);
end
cum_DCA(:,n) = cumsum(sorted_DCA(:,4));


% Run nbZPX2
nbZPX2 = nbZPX2_ALL(:,:,n);
nbZPX2 = nbZPX2 - min(nbZPX2(:));
% nbZPX2 = nbZPX2.*gW(:,:,n);
sorted_nbZPX2 = sort_matrix_descend_2(nbZPX2,1);
sorted_nbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_nbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_nbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_nbZPX2(:,n) = cumsum(sorted_nbZPX2(:,4));


% Run dbZPX2
dbZPX2 = dbZPX2_ALL(:,:,n);
dbZPX2 = dbZPX2 - min(dbZPX2(:));
% dbZPX2 = dbZPX2.*gW(:,:,n);
sorted_dbZPX2 = sort_matrix_descend_2(dbZPX2,1);
sorted_dbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_dbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_dbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_dbZPX2(:,n) = cumsum(sorted_dbZPX2(:,4));


% Run dgbZPX2
dgbZPX2 = dgbZPX2_ALL(:,:,n);
dgbZPX2 = dgbZPX2 - min(dgbZPX2(:));
% dgbZPX2 = dgbZPX2.*gW(:,:,n);
sorted_dgbZPX2 = sort_matrix_descend_2(dgbZPX2,1);
sorted_dgbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_dgbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_dgbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_dgbZPX2(:,n) = cumsum(sorted_dgbZPX2(:,4));


% Run gbZPX2
gbZPX2 = gbZPX2_ALL(:,:,n);
gbZPX2 = gbZPX2 - min(gbZPX2(:));
gbZPX2 = gbZPX2.*gW(:,:,n);
sorted_gbZPX2 = sort_matrix_descend_2(gbZPX2,1);
sorted_gbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_gbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_gbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_gbZPX2(:,n) = cumsum(sorted_gbZPX2(:,4));


% Run fgbZPX2
fgbZPX2 = fgbZPX2_ALL(:,:,n);
fgbZPX2 = fgbZPX2 - min(fgbZPX2(:));
% fgbZPX2 = fgbZPX2.*gW(:,:,n);
sorted_fgbZPX2 = sort_matrix_descend_2(fgbZPX2,1);
sorted_fgbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_fgbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_fgbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_fgbZPX2(:,n) = cumsum(sorted_fgbZPX2(:,4));


% Run RSEM_SCA
RSEM_SCA = RSEM_SCA_ALL(:,:,n);
RSEM_SCA = RSEM_SCA - min(RSEM_SCA(:));
% RSEM_SCA = RSEM_SCA.*gW(:,:,n);
sorted_RSEM_SCA = sort_matrix_descend_2(RSEM_SCA,1);
sorted_RSEM_SCA(:,4) = 0;
[~,ia,ib] = intersect(sorted_RSEM_SCA(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_RSEM_SCA(ia(i),4) = covarions(ib(i),1);
end
cum_RSEM_SCA(:,n) = cumsum(sorted_RSEM_SCA(:,4));


% % Run RAMA_SCA
% RAMA_SCA = RAMA_SCA_ALL(:,:,n);
% RAMA_SCA = RAMA_SCA - min(RAMA_SCA(:));
% RAMA_SCA = RAMA_SCA.*gW(:,:,n);
% sorted_RAMA_SCA = sort_matrix_descend_2(RAMA_SCA,1);
% sorted_RAMA_SCA(:,4) = 0;
% [~,ia,ib] = intersect(sorted_RAMA_SCA(:,2:3),covarions(:,2:3),'rows');    
% npairs=size(ib);
% for i = 1:npairs
%     sorted_RAMA_SCA(ia(i),4) = covarions(ib(i),1);
% end
% cum_RAMA_SCA(:,n) = cumsum(sorted_RAMA_SCA(:,4));


% Run fodorSCA
fodorSCA = fodorSCA_ALL(:,:,n);
fodorSCA = fodorSCA - min(fodorSCA(:));
% fodorSCA = fodorSCA.*gW(:,:,n);
sorted_fodorSCA = sort_matrix_descend_2(fodorSCA,1);
sorted_fodorSCA(:,4) = 0;
[~,ia,ib] = intersect(sorted_fodorSCA(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_fodorSCA(ia(i),4) = covarions(ib(i),1);
end
cum_fodorSCA(:,n) = cumsum(sorted_fodorSCA(:,4));


% Run OMES
OMES = OMES_ALL(:,:,n);
OMES = OMES - min(OMES(:));
% OMES = OMES.*gW(:,:,n);
sorted_OMES = sort_matrix_descend_2(OMES,1);
sorted_OMES(:,4) = 0;
[~,ia,ib] = intersect(sorted_OMES(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_OMES(ia(i),4) = covarions(ib(i),1);
end
cum_OMES(:,n) = cumsum(sorted_OMES(:,4));


% Run McBASC
McBASC = McBASC_ALL(:,:,n);
McBASC = McBASC - min(McBASC(:));
% McBASC = McBASC.*gW(:,:,n);
sorted_McBASC = sort_matrix_descend_2(McBASC,1);
sorted_McBASC(:,4) = 0;
[~,ia,ib] = intersect(sorted_McBASC(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_McBASC(ia(i),4) = covarions(ib(i),1);
end
cum_McBASC(:,n) = cumsum(sorted_McBASC(:,4));


% Run ELSC
ELSC = ELSC_ALL(:,:,n);
ELSC = ELSC - min(ELSC(:));
% ELSC = ELSC.*gW(:,:,n);
sorted_ELSC = sort_matrix_descend_2(ELSC,1);
sorted_ELSC(:,4) = 0;
[~,ia,ib] = intersect(sorted_ELSC(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_ELSC(ia(i),4) = covarions(ib(i),1);
end
cum_ELSC(:,n) = cumsum(sorted_ELSC(:,4));

end

% Here we get all the mean of all cumulative sums.
cum_covarions = nanmean(cum_covarions(:,1:n),2);
cum_DCA = nanmean(cum_DCA(:,1:n),2);
cum_MI = nanmean(cum_MI(:,1:n),2);
cum_MIP = nanmean(cum_MIP(:,1:n),2);
cum_logR = nanmean(cum_logR(:,1:n),2);
cum_ZPX2 = nanmean(cum_ZPX2(:,1:n),2);
cum_nbZPX2 = nanmean(cum_nbZPX2(:,1:n),2);
cum_dbZPX2 = nanmean(cum_dbZPX2(:,1:n),2);
cum_dgbZPX2 = nanmean(cum_dgbZPX2(:,1:n),2);
cum_gbZPX2 = nanmean(cum_gbZPX2(:,1:n),2);
cum_fgbZPX2 = nanmean(cum_fgbZPX2(:,1:n),2);
cum_RSEM_SCA = nanmean(cum_RSEM_SCA(:,1:n),2);
% cum_RAMA_SCA = nanmean(cum_RAMA_SCA(:,1:n),2);
cum_fodorSCA = nanmean(cum_fodorSCA(:,1:n),2);
cum_OMES = nanmean(cum_OMES(:,1:n),2);
cum_McBASC = nanmean(cum_McBASC(:,1:n),2);
cum_ELSC = nanmean(cum_ELSC(:,1:n),2);

%%
COEV_methods = figure; 
    	set(COEV_methods,'Units','normalized','Position',[0 0.2 0.8 0.6 ],...
    	'Name','COEV Methods: sensitivity'); clf;
	semilogx(cum_covarions,'-b','LineWidth',2);hold on
    semilogx(cum_MI,'-r','LineWidth',2);
    % semilogx(cum_MIP,'-g','LineWidth',2);
    semilogx(cum_logR,'-','Color',[0.3,0.6,0.2],'LineWidth',2);
    semilogx(cum_ZPX2,'-m','LineWidth',2);
    semilogx(cum_DCA,'-g','LineWidth',2);
    semilogx(cum_nbZPX2,'-','Color','c','LineWidth',2);
	semilogx(cum_dbZPX2,'-y','LineWidth',2);
% 	semilogx(cum_fgbZPX2,'-k','LineWidth',2);
	semilogx(cum_dgbZPX2,'-','Color',[0.0,0.5,1.0],'LineWidth',2);
%    semilogx(cum_gbZPX2,'-','Color',[1,0.5,0],'LineWidth',2);
	semilogx(cum_RSEM_SCA,'--','Color','y','LineWidth',2);
% 	semilogx(cum_RAMA_SCA,'--','Color','y','LineWidth',2);
	semilogx(cum_fodorSCA,'--','Color','c','LineWidth',2);
	semilogx(cum_OMES,'--','Color','g','LineWidth',2);
	semilogx(cum_McBASC,'--','Color','m','LineWidth',2);
	semilogx(cum_ELSC,'--','Color',[0.0,0.5,1.0],'LineWidth',2);
    plot(ncov,(-100:5E2:40E3),'ok','LineWidth',2,'MarkerSize',3);
	hold off
    % set(gca,'Xlim',[0,1.0],'Ylim',[0,60]);
    % set(gca,'Xlim',[0,0.6E2],'Ylim',[-1E2,4E3]);
    set(gca,'Xlim',[0,1E2],'Ylim',[-1E2,0.45E4]);
    % set(gca,'Xlim',[0,1E2],'Ylim',[-1E2,4.0E4]);

	legend('Ideal','MI','logR','ZPX2',...
        'DCA','nbZPX2','dbZPX2','dgbZPX2',...
        'ramaSCA','fodorSCA','OMES','McBASC','ELSC','Location','Best');
	title('MI Methods','FontSize',14,'FontWeight','n');
	xlabel('Number of top pairs identified by a method');
	ylabel('Cumulative count of true covariation events');
    box('on');
    grid('on');

    


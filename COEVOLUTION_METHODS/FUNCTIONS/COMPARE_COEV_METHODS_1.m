% This script compares the performance of several MI methods (MI,MIP,ZPX2,
% dbZPX2,DCA on simulated data sets. In this case we import directly a
% covar_vec from MSAvolve or Simprot.

%[~,nmsa] = faln_to_nmsa('testcorr320.aln');
%[~,nmsa] = faln_to_nmsa('atp12_phy_sp.aln');
[~,nmsa] = faln_to_nmsa('bt160_2.aln');
tree = 'atp12_95.phy.tre'; %treefile 
totsim = 100; %number of simulations;

import_covar_vec('covar_vec.txt');
ncov = size(covar_vec,1);

% Here we calculate the expected cumulative score
covar_vec = sortrows(covar_vec, -3);
covar_vec(:,4) = cumsum(covar_vec(:,3));

% Run DCA
tStart = tic;
[ DCA ] = NMSA_to_DCA( nmsa,20 );
tElapsed = toc(tStart);
fprintf('DCA execution time: %7.2f s \n',tElapsed);

sorted_DCA = sort_matrix_descend(DCA);
sorted_DCA(:,4) = 0;
[c,ia,ib] = intersect(sorted_DCA(:,2:3),covar_vec(:,1:2),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_DCA(ia(i),4) = covar_vec(ib(i),3);
end
sorted_DCA(:,5)= cumsum(sorted_DCA(:,4));
end

% Run ZPX2
tStart = tic;
[MI] = NMSA_to_MI(nmsa);
[MIP] = MI_to_MIP(MI);
[~,ZPX2] = MI_to_ZPX(MIP);
tElapsed = toc(tStart);
fprintf('ZPX2 execution time: %7.2f s \n',tElapsed);

sorted_ZPX2 = sort_matrix_descend(ZPX2);
sorted_ZPX2(:,4) = 0;
[c,ia,ib] = intersect(sorted_ZPX2(:,2:3),covar_vec(:,1:2),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_ZPX2(ia(i),4) = covar_vec(ib(i),3);
end
sorted_ZPX2(:,5)= cumsum(sorted_ZPX2(:,4));


% Run MI
sorted_MI = sort_matrix_descend(MI);
sorted_MI(:,4) = 0;
[c,ia,ib] = intersect(sorted_MI(:,2:3),covar_vec(:,1:2),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_MI(ia(i),4) = covar_vec(ib(i),3);
end
sorted_MI(:,5)= cumsum(sorted_MI(:,4));


% Run MIP
sorted_MIP = sort_matrix_descend(MIP);
sorted_MIP(:,4) = 0;
[c,ia,ib] = intersect(sorted_MIP(:,2:3),covar_vec(:,1:2),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_MIP(ia(i),4) = covar_vec(ib(i),3);
end
sorted_MIP(:,5)= cumsum(sorted_MIP(:,4));


% Run dbZPX2
tStart = tic;
[dbZPX2] = NMSA_to_dbZPX2( nmsa );
tElapsed = toc(tStart);
fprintf('dbZPX2 execution time: %7.2f s \n',tElapsed);

sorted_dbZPX2 = sort_matrix_descend(dbZPX2);
sorted_dbZPX2(:,4) = 0;
[c,ia,ib] = intersect(sorted_dbZPX2(:,2:3),covar_vec(:,1:2),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_dbZPX2(ia(i),4) = covar_vec(ib(i),3);
end
npairs = size(sorted_dbZPX2,1);
sorted_dbZPX2(:,5)= cumsum(sorted_dbZPX2(:,4));


% Here we plot 
	COEV_methods = figure; 
    	set(COEV_methods,'Units','normalized','Position',[0 0.2 0.8 0.6 ],...
    	'Name','COEV Methods: sensitivity'); clf;
	semilogx(covar_vec(:,4),'b','LineWidth',2);hold on
	semilogx(sorted_MI(:,5),'m','LineWidth',2);
	semilogx(sorted_MIP(:,5),'c','LineWidth',2);
	semilogx(sorted_ZPX2(:,5),':y','LineWidth',2);
	semilogx(sorted_DCA(:,5),':r','LineWidth',2);
	semilogx(sorted_dbZPX2(:,5),':k','LineWidth',1.5);
	hold off
	%set(gca,'Xlim',[0,1.0],'Ylim',[0,60]);
	legend('Ideal','MI','MIP','ZPX2','DCA','dbZPX2');
	title('MI Methods','FontSize',14,'FontWeight','n');
	xlabel('Number of pairs');
	ylabel('Cumulative identification of true covarions');



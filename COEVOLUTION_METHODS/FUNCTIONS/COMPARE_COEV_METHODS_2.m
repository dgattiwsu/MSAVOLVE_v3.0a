% This script can be called after loading a simulation run of MSAvolve 
% to compare the performance of several MI methods on simulated 
% data sets. In this case we use directly the COV or cov_COV matrix rather  
% than the original vector of covarions 'covar_vec'. Results from all the 
% runs are averaged.

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


for n = 1:end_cycle
% We can run the comparison on the total COV matrix ...
% covarions = sort_matrix_descend_2(COV_ALL(:,:,n),1);
% ... or just the true covarions matrix.
covarions = sort_matrix_descend_2(cov_COV_ALL(:,:,n),1);

cum_covarions(:,n) = cumsum(covarions(:,1));

% Run MI
sorted_MI = sort_matrix_descend_2(MI_ALL(:,:,n),1);
sorted_MI(:,4) = 0;
[~,ia,ib] = intersect(sorted_MI(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_MI(ia(i),4) = covarions(ib(i),1);
end
cum_MI(:,n) = cumsum(sorted_MI(:,4));


% Run MIP
sorted_MIP = sort_matrix_descend_2(abs(MIP_ALL(:,:,n)),1);
sorted_MIP(:,4) = 0;
[~,ia,ib] = intersect(sorted_MIP(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_MIP(ia(i),4) = covarions(ib(i),1);
end
cum_MIP(:,n) = cumsum(sorted_MIP(:,4));


% Run logR
sorted_logR = sort_matrix_descend_2(abs(logR_ALL(:,:,n)),1);
sorted_logR(:,4) = 0;
[~,ia,ib] = intersect(sorted_logR(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_logR(ia(i),4) = covarions(ib(i),1);
end
cum_logR(:,n) = cumsum(sorted_logR(:,4));


% Run ZPX2
sorted_ZPX2 = sort_matrix_descend_2(abs(ZPX2_ALL(:,:,n)),1);
sorted_ZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_ZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_ZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_ZPX2(:,n) = cumsum(sorted_ZPX2(:,4));


% Run DCA
sorted_DCA = sort_matrix_descend_2(DCA_ALL(:,:,n),1);
sorted_DCA(:,4) = 0;
[~,ia,ib] = intersect(sorted_DCA(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_DCA(ia(i),4) = covarions(ib(i),1);
end
cum_DCA(:,n) = cumsum(sorted_DCA(:,4));


% Run nbZPX2
sorted_nbZPX2 = sort_matrix_descend_2(abs(nbZPX2_ALL(:,:,n)),1);
sorted_nbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_nbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_nbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_nbZPX2(:,n) = cumsum(sorted_nbZPX2(:,4));


% Run dbZPX2
sorted_dbZPX2 = sort_matrix_descend_2(abs(dbZPX2_ALL(:,:,n)),1);
sorted_dbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_dbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_dbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_dbZPX2(:,n) = cumsum(sorted_dbZPX2(:,4));


% Run dgbZPX2
sorted_dgbZPX2 = sort_matrix_descend_2(abs(dgbZPX2_ALL(:,:,n)),1);
sorted_dgbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_dgbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_dgbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_dgbZPX2(:,n) = cumsum(sorted_dgbZPX2(:,4));


% Run gbZPX2
sorted_gbZPX2 = sort_matrix_descend_2(abs(gbZPX2_ALL(:,:,n)),1);
sorted_gbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_gbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_gbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_gbZPX2(:,n) = cumsum(sorted_gbZPX2(:,4));


% Run fgbZPX2
sorted_fgbZPX2 = sort_matrix_descend_2(abs(fgbZPX2_ALL(:,:,n)),1);
sorted_fgbZPX2(:,4) = 0;
[~,ia,ib] = intersect(sorted_fgbZPX2(:,2:3),covarions(:,2:3),'rows');    
npairs=size(ib);
for i = 1:npairs
    sorted_fgbZPX2(ia(i),4) = covarions(ib(i),1);
end
cum_fgbZPX2(:,n) = cumsum(sorted_fgbZPX2(:,4));


end

% Here we get all the mean of all cumulative sums.
cum_covarions = mean(cum_covarions(:,1:n),2);
cum_DCA = mean(cum_DCA(:,1:n),2);
cum_MI = mean(cum_MI(:,1:n),2);
cum_MIP = mean(cum_MIP(:,1:n),2);
cum_logR = mean(cum_logR(:,1:n),2);
cum_ZPX2 = mean(cum_ZPX2(:,1:n),2);
cum_nbZPX2 = mean(cum_nbZPX2(:,1:n),2);
cum_dbZPX2 = mean(cum_dbZPX2(:,1:n),2);
cum_dgbZPX2 = mean(cum_dgbZPX2(:,1:n),2);
cum_gbZPX2 = mean(cum_gbZPX2(:,1:n),2);
cum_fgbZPX2 = mean(cum_fgbZPX2(:,1:n),2);


COEV_methods = figure; 
    	set(COEV_methods,'Units','normalized','Position',[0 0.2 0.6 0.6 ],...
    	'Name','COEV Methods: sensitivity'); clf;
	semilogx(cum_covarions,'-b','LineWidth',2);hold on
    semilogx(cum_MI,'-r','LineWidth',2);
    % semilogx(cum_MIP,'-g','LineWidth',2);
    semilogx(cum_logR,'-','Color',[0.3,0.6,0.2],'LineWidth',2);
    semilogx(cum_ZPX2,'-m','LineWidth',2);
    semilogx(cum_DCA,'-g','LineWidth',2);
    semilogx(cum_nbZPX2,'-','Color','c','LineWidth',2);
	semilogx(cum_dbZPX2,'-y','LineWidth',2);
	semilogx(cum_fgbZPX2,'-k','LineWidth',2);
	semilogx(cum_dgbZPX2,'-','Color',[0.0,0.5,1.0],'LineWidth',2);
    semilogx(cum_gbZPX2,'-','Color',[1,0.5,0],'LineWidth',2);
    plot(ncov,(-100:15E2:40E3),'ok','LineWidth',2,'MarkerSize',3);
	hold off
    %set(gca,'Xlim',[0,1.0],'Ylim',[0,60]);
    set(gca,'Xlim',[0,0.6E2],'Ylim',[-1E2,4E3]);
    set(gca,'Xlim',[0,1E2],'Ylim',[-1E2,2.5E4]);
    set(gca,'Xlim',[0,1E2],'Ylim',[-1E2,4.0E4]);

	legend('Ideal','MI','logR','ZPX2',...
        'DCA','nbZPX2','dbZPX2','fgbZPX2','dgbZPX2','gbZPX2',...
        'Location','Best');
	title('MI Methods','FontSize',14,'FontWeight','n');
	xlabel('Number of top pairs identified by a method');
	ylabel('Cumulative count of true covariations');
    


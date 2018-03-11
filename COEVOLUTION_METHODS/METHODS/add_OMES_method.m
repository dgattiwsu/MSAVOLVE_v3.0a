% This script adds on the fly the methods included in Fodor's OMES package.
% We calculate also MI also using the OMES package as a test that the package is 
% working correctly. This package is written in Java, therefore ...
% IMPORTANT: remember to set the CLASSPATH env variable to the directory
% where the OMES package was installed. For example, if the package OMES is
% in /usr/local/OMES add the following to your .tcshrc or .cshrc 
% setenv CLASSPATH /usr/local/OMES

OMES_ALL = zeros(npos,npos,end_cycle);
ELSC_ALL = zeros(npos,npos,end_cycle);
McBASC_ALL = zeros(npos,npos,end_cycle);
fodorMI_ALL = zeros(npos,npos,end_cycle);
fodorSCA_ALL = zeros(npos,npos,end_cycle);
fcov_OMES = zeros(end_cycle,3);
fcov_ELSC = zeros(end_cycle,3);
fcov_McBASC = zeros(end_cycle,3);
fcov_fodorMI = zeros(end_cycle,3);
fcov_fodorSCA = zeros(end_cycle,3);
cov_zscore_OMES = zeros(end_cycle,ncov_coord);
covar_vec_zscore_OMES = zeros(ncov_coord,3,end_cycle);
cov_zscore_ELSC = zeros(end_cycle,ncov_coord);
covar_vec_zscore_ELSC = zeros(ncov_coord,3,end_cycle);
cov_zscore_McBASC = zeros(end_cycle,ncov_coord);
covar_vec_zscore_McBASC = zeros(ncov_coord,3,end_cycle);
cov_zscore_fodorMI = zeros(end_cycle,ncov_coord);
covar_vec_zscore_fodorMI = zeros(ncov_coord,3,end_cycle);
cov_zscore_fodorSCA = zeros(end_cycle,ncov_coord);
covar_vec_zscore_fodorSCA = zeros(ncov_coord,3,end_cycle);
corr_cov_zscore_OMES = zeros(end_cycle,3);
corr_cov_zscore_McBASC = zeros(end_cycle,3);
corr_cov_zscore_ELSC = zeros(end_cycle,3);
corr_cov_zscore_fodorMI = zeros(end_cycle,3);
corr_cov_zscore_fodorSCA = zeros(end_cycle,3);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('OMES package cycle %d \n', n);

% We delete some old files and we print the alignment in selex format. 

    delete('*.matrix');
    delete('*.selex');
    
    [rows,~] = size(MSA_select_ALL(:,:,n));
    c_MSA_select = int2aa(MSA_select_ALL(:,:,n));
    
    jc_MSA0 = ['Sequence_' int2str(0) '   '];
    for i = 1:rows
    jc_MSA = ['Sequence_' int2str(i)];
    jc_MSA0 = char(jc_MSA0,jc_MSA);
    end
    jc_MSA0 = jc_MSA0(2:end,:);  
    jc_MSA_select = [jc_MSA0 c_MSA_select];
    
    fid = fopen('MSA_select.selex','w');
    for i = 1:rows
        fprintf(fid, '%s\n', jc_MSA_select(i,:));
    end
    fclose(fid);

    clear jc_MSA jc_MSA0
    
 !java covariance.algorithms.MICovariance MSA_select.selex fodorMI.matrix
 import_omes_matrix('fodorMI.matrix');
 data(:,1:2) = data(:,1:2) + 1;
 fodorMI = zeros(npos,npos);
 for k = 1:size(data,1)
    i = data(k,1);
    j = data(k,2);
        fodorMI(i,j) = data(k,3);
 end
 fodorMI = fodorMI + fodorMI';
 clear data     
     
 !java covariance.algorithms.OmesCovariance MSA_select.selex OMES.matrix
 import_omes_matrix('OMES.matrix');
 data(:,1:2) = data(:,1:2) + 1;
 OMES = zeros(npos,npos);
 for k = 1:size(data,1)
    i = data(k,1);
    j = data(k,2);
        OMES(i,j) = data(k,3);
 end
 OMES = OMES + OMES';
 clear data

 !java covariance.algorithms.ELSCCovariance MSA_select.selex ELSC.matrix
 import_omes_matrix('ELSC.matrix');
 data(:,1:2) = data(:,1:2) + 1;
 ELSC = zeros(npos,npos);
 for k = 1:size(data,1)
    i = data(k,1);
    j = data(k,2);
        ELSC(i,j) = data(k,3);
 end
 ELSC = ELSC + ELSC';
 clear data
 
 !java covariance.algorithms.McBASCCovariance MSA_select.selex McBASC.matrix
 import_omes_matrix('McBASC.matrix');
 data(:,1:2) = data(:,1:2) + 1;
 McBASC = zeros(npos,npos);
 for k = 1:size(data,1)
    i = data(k,1);
    j = data(k,2);
        McBASC(i,j) = data(k,3);
 end
 McBASC = McBASC + McBASC';
 clear data 

 !java covariance.algorithms.JavaSCA MSA_select.selex fodorSCA.matrix
 import_omes_matrix('fodorSCA.matrix');
 data(:,1:2) = data(:,1:2) + 1;
 fodorSCA = zeros(npos,npos);
 for k = 1:size(data,1)
    i = data(k,1);
    j = data(k,2);
        fodorSCA(i,j) = data(k,3);
 end
 fodorSCA = fodorSCA + fodorSCA';
 clear data

 for i = 1:npos
     fodorMI(i,i) = NaN;
     fodorSCA(i,i) = NaN;
     OMES(i,i) = NaN;
     ELSC(i,i) = NaN;
     McBASC(i,i) = NaN;
 end

 OMES_ALL(:,:,n) = OMES;
 ELSC_ALL(:,:,n) = ELSC;
 McBASC_ALL(:,:,n) = McBASC;
 fodorMI_ALL(:,:,n) = fodorMI;
 fodorSCA_ALL(:,:,n) = fodorSCA;

 [fcov_OMES(n,:),cov_zscore_OMES(n,:),...
     covar_vec_zscore_OMES(:,:,n)] = ... 
     coev_stats_2(OMES,ncov,covar_vec);
 [fcov_ELSC(n,:),cov_zscore_ELSC(n,:),...
     covar_vec_zscore_ELSC(:,:,n)] = ... 
     coev_stats_2(ELSC,ncov,covar_vec);
 [fcov_McBASC(n,:),cov_zscore_McBASC(n,:),...
     covar_vec_zscore_McBASC(:,:,n)] = ... 
     coev_stats_2(McBASC,ncov,covar_vec);
 [fcov_fodorMI(n,:),cov_zscore_fodorMI(n,:),...
     covar_vec_zscore_fodorMI(:,:,n)] = ... 
     coev_stats_2(fodorMI,ncov,covar_vec);
 [fcov_fodorSCA(n,:),cov_zscore_fodorSCA(n,:),...
     covar_vec_zscore_fodorSCA(:,:,n)] = ... 
     coev_stats_2(fodorSCA,ncov,covar_vec);  
   
end

stat_fcov_ELSC = [mean(fcov_ELSC);std(fcov_ELSC)];
stat_fcov_McBASC = [mean(fcov_McBASC);std(fcov_McBASC)];
stat_fcov_OMES = [mean(fcov_OMES);std(fcov_OMES)];
stat_fcov_fodorMI = [mean(fcov_fodorMI);std(fcov_fodorMI)];
stat_fcov_fodorSCA = [mean(fcov_fodorSCA);std(fcov_fodorSCA)];

mean_cov_zscore_OMES = nanmean(cov_zscore_OMES,2);
mean_cov_zscore_McBASC = nanmean(cov_zscore_McBASC,2);
mean_cov_zscore_ELSC = nanmean(cov_zscore_ELSC,2);
mean_cov_zscore_fodorMI = nanmean(cov_zscore_fodorMI,2);
mean_cov_zscore_fodorSCA = nanmean(cov_zscore_fodorSCA,2);

for i = 1:end_cycle
corr_cov_zscore_OMES(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_OMES(i,:)','type','Pearson');
corr_cov_zscore_McBASC(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_McBASC(i,:)','type','Pearson');
corr_cov_zscore_ELSC(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_ELSC(i,:)','type','Pearson');
corr_cov_zscore_fodorMI(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorMI(i,:)','type','Pearson');
corr_cov_zscore_fodorSCA(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorSCA(i,:)','type','Pearson');

corr_cov_zscore_OMES(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_OMES(i,:)','type','Spearman');
corr_cov_zscore_McBASC(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_McBASC(i,:)','type','Spearman');
corr_cov_zscore_ELSC(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_ELSC(i,:)','type','Spearman');
corr_cov_zscore_fodorMI(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorMI(i,:)','type','Spearman');
corr_cov_zscore_fodorSCA(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorSCA(i,:)','type','Spearman');

corr_cov_zscore_OMES(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_OMES(i,:)','type','Kendall');
corr_cov_zscore_McBASC(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_McBASC(i,:)','type','Kendall');
corr_cov_zscore_ELSC(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_ELSC(i,:)','type','Kendall');
corr_cov_zscore_fodorMI(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorMI(i,:)','type','Kendall');
corr_cov_zscore_fodorSCA(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_fodorSCA(i,:)','type','Kendall');
end




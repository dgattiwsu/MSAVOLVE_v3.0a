% This script adds on the fly the PSICOV methods.
% IMPORTANT: remember to bring the program psicov in your working directory. 
% The program is in MSAVOLVE_v1.0a/COEVOLUTION_METHODS/FUNCTIONS/PSICOV/ and must be
% compiled according to the directions listed in the README file. We
% provide a linux executable in which the limitation for a minimum number
% of 500 sequences has been removed.

PSICOV_ALL = zeros(npos,npos,end_cycle);
fcov_PSICOV = zeros(end_cycle,3);
cov_zscore_PSICOV = zeros(end_cycle,ncov_coord);
covar_vec_zscore_PSICOV = zeros(ncov_coord,3,end_cycle);
corr_cov_zscore_PSICOV = zeros(end_cycle,3);

[nrows,ncols] = size(MSA_select);

for n = 1:end_cycle
fprintf('PSICOV method cycle %d \n', n);

% We delete some old files and we print the alignment in selex format. 

    delete('*.matrix');
    delete('*.selex');
    
    [rows,~] = size(MSA_select_ALL(:,:,n));
    c_MSA_select = int2aa(MSA_select_ALL(:,:,n));
            
    fid = fopen('MSA_select.selex','w');
    for i = 1:rows
        fprintf(fid, '%s\n', c_MSA_select(i,:));
    end
    fclose(fid);

    clear jc_MSA jc_MSA0
    
% !./psicov -p -d 0.03 MSA_select.selex > PSICOV.matrix
% !./psicov -d 0.03 -r 0.005 -i 62 MSA_select.selex > PSICOV.matrix
% We run here in some relaxed condition to favor convergence
 !./psicov -a -r 0.005 -i 62 MSA_select.selex > psicov.matrix
 import_psicov_matrix('psicov.matrix');
 PSICOV = zeros(npos,npos);
 for k = 1:size(psicov.data,1)
    i = psicov.data(k,1);
    j = psicov.data(k,2);
        PSICOV(i,j) = psicov.data(k,5);
 end
 PSICOV = PSICOV + PSICOV';
 clear psicov

 for i = 1:npos
     PSICOV(i,i) = NaN;
 end

 PSICOV_ALL(:,:,n) = PSICOV;

 [fcov_PSICOV(n,:),cov_zscore_PSICOV(n,:),...
     covar_vec_zscore_PSICOV(:,:,n)] = ... 
     coev_stats_2(PSICOV,ncov,covar_vec);
   
end

stat_fcov_PSICOV = [mean(fcov_PSICOV);std(fcov_PSICOV)];
mean_cov_zscore_PSICOV = nanmean(cov_zscore_PSICOV,2);

for i = 1:end_cycle
corr_cov_zscore_PSICOV(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_PSICOV(i,:)','type','Pearson');
corr_cov_zscore_PSICOV(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_PSICOV(i,:)','type','Spearman');
corr_cov_zscore_PSICOV(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_PSICOV(i,:)','type','Kendall');
end




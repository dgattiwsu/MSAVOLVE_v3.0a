% This script is used to add on the fly the GREMLIN method. 

corr_cov_zscore_GREMLIN = zeros(end_cycle,3);
GREMLIN_ALL = zeros(npos,npos,end_cycle);
fcov_GREMLIN = zeros(end_cycle,3);
cov_zscore_GREMLIN = zeros(end_cycle,ncov);
covar_vec_zscore_GREMLIN = zeros(ncov,3,end_cycle);

[nrows,ncols] = size(MSA_select_ALL(:,:,1));

for n = 1:end_cycle
fprintf('GREMLIN cycle %d \n', n);

 % Here we write out the msa in selex format.   
 
        c_MSA_select = int2aa(MSA_select_ALL(:,:,n));          
        fid = fopen('MSA_select.selex','w');
        for i = 1:nrows
            fprintf(fid, '%s\n', c_MSA_select(i,:));
        end
        fclose(fid);        
        gremlin('MSA_select.selex','GREMLIN.matrix')
        
        GREMLIN = importdata('GREMLIN.matrix');

        for i = 1:ncols
            GREMLIN(i,i) = NaN;
        end
 
 GREMLIN_ZPX2 = MIP_to_ZPX2(GREMLIN);
 
% Gap correction
GREMLIN_ZPX2_1 = GREMLIN_ZPX2;
GREMLIN_ZPX2_2 = GREMLIN_ZPX2_1 - min(GREMLIN_ZPX2_1(:));
GREMLIN_ZPX2 = GREMLIN_ZPX2_2.*gW(:,:,n);
 
 % Here we calculate the usual statistics for GREMLIN.
 GREMLIN_ALL(:,:,n) = GREMLIN_ZPX2;
 
 [fcov_GREMLIN(n,:),cov_zscore_GREMLIN(n,:),...
     covar_vec_zscore_GREMLIN(:,:,n)] = ... 
     coev_stats_2(GREMLIN_ZPX2,ncov,covar_vec); 
 
 % Here we remove the internal and external files generated at each cycle.
 clear c_MSA_select GREMLIN
 !rm MSA_select.selex GREMLIN.matrix
  
end

 stat_fcov_GREMLIN = [mean(fcov_GREMLIN);std(fcov_GREMLIN)];
 mean_cov_zscore_GREMLIN = nanmean(cov_zscore_GREMLIN,2);
 
for i = 1:end_cycle
 corr_cov_zscore_GREMLIN(i,1) = corr(cov_zscore_COV(i,:)',cov_zscore_GREMLIN(i,:)','type','Pearson');
 corr_cov_zscore_GREMLIN(i,2) = corr(cov_zscore_COV(i,:)',cov_zscore_GREMLIN(i,:)','type','Spearman');
 corr_cov_zscore_GREMLIN(i,3) = corr(cov_zscore_COV(i,:)',cov_zscore_GREMLIN(i,:)','type','Kendall');
end



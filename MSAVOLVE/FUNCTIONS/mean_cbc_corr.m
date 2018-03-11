function [ mean_cbc_corr ] = mean_cbc_corr( mat1,mat2 )
% This function calculate the mean column-by-column correlation between two
% MI matrices

mean_cbc_corr=nanmean(diag(corr(mat1,mat2,'rows','pairwise')));

    end
    

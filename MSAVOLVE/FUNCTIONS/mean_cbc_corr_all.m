function [ mean_cbc_corr ] = mean_cbc_corr_all( mat1,mat2 )
% This function calculate the mean column-by-column correlation between two
% MI matrices
mean_cbc_corr = zeros(1,3);
[~,cols] = size(mat1);
pearson = zeros(cols,1);
spearman = zeros(cols,1);
kendall = zeros(cols,1);
for i = 1:cols
%    A = mat1(i+1:rows,i);
%    B = mat2(i+1:rows,i);
    A = mat1(:,i);
    B = mat2(:,i);
    pearson(i,1) = corr(A,B,'rows','pairwise','type','Pearson');
    spearman(i,1) = corr(A,B,'rows','pairwise','type','Spearman');
    kendall(i,1) = corr(A,B,'rows','pairwise','type','Kendall');
end    
mean_cbc_corr(1,1)=nanmean(pearson);
mean_cbc_corr(1,2)=nanmean(spearman);
mean_cbc_corr(1,3)=nanmean(kendall);

    end
    

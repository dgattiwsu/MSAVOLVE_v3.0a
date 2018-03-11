function [ fcov,cov_zscore,covar_vec_zscore ] = coev_stats_3(coev,ncov,covar_vec)
% This function calculates the fraction of the ncov top zscores in each 
% coev matrix that corresponds to covarying positions.
 mean_coev = nanmean(coev(1:end));
 std_coev = nanstd(coev(1:end));
 z_coev = (coev - mean_coev)/std_coev;
 oz_coev = sort_matrix_descend(z_coev);
 oz_coev(:,4) = 0;
    for i = 1:ncov
        ind = oz_coev(:,2) == covar_vec(i,1) ... 
            & oz_coev(:,3) == covar_vec(i,2) ;
        oz_coev(ind,4) = 1;
    end
 % Which percentage of the top 10 zscores in the COV matrix belong to the
 % covar_vec.
 fcov(1) = sum(oz_coev(1:10,4))/10;
 % Which percentage of the top N zscores (where N is the length of
 % covar_vec) belong to the covar_vec.
 fcov(2) = sum(oz_coev(1:ncov,4))/ncov;
 % Which percentage of the top 100 zscores in the COV matrix belong to the
 % covar_vec.
 fcov(3) = sum(oz_coev(1:100,4))/100; 
 
 % "cov_zscore" is the zscore in the COV matrix for each of the covarying
 % positions in covar_vec. Only the score is reported.
 cov_zscore = zeros(1,ncov);
    for i = 1:ncov
        ind = covar_vec(i,:);
        cov_zscore(1,i) = z_coev(ind(1),ind(2));
    end
 % "covar_vec_zscore" is the zscore in the COV matrix for each of the 
 % covarying positions in covar_vec. Both the covarying positions (the 
 % covar_vec) and the scores are reported.
    covar_vec_zscore = ...
     [ cov_zscore' covar_vec ];

end


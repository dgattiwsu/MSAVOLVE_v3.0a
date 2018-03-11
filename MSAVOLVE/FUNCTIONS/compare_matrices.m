MAT1 = mutcov_COV_bk;
[MAT1_V,MAT1_D] = spectral_2(MAT1,3);
% figure
% scatter3(MAT1_V(:,1),MAT1_V(:,2),MAT1_V(:,3))
ncoef = size(MAT1_V,1);
for i = 1:ncoef
    MAT1_V_norm(i,1) = norm(MAT1_V(i,:));
end
ind = MAT1_V_norm>0.05*max(MAT1_V_norm);
indf = find(ind);
MAT1_V_shell = MAT1_V(indf,:);
% figure
% scatter3(MAT1_V_shell(:,1),MAT1_V_shell(:,2),MAT1_V_shell(:,3))
MAT1_shell = MAT1_V_shell*MAT1_D*MAT1_V_shell';
% figure
% imagesc(MAT1_shell);


MAT2 = cov_COV_bk;
[MAT2_V,MAT2_D] = spectral_2(MAT2,3);
% figure
% scatter3(MAT2_V(:,1),MAT2_V(:,2),MAT2_V(:,3))
% ncoef = size(MAT2_V,1);
% for i = 1:ncoef
%     MAT2_V_norm(i,1) = norm(MAT2_V(i,:));
% end
% r_ind = MAT2_V_norm>0.1;
% rindf = find(r_ind);
MAT2_V_shell = MAT2_V(indf,:);
% figure
% scatter3(MAT2_V_shell(:,1),MAT2_V_shell(:,2),MAT2_V_shell(:,3))
MAT2_shell = MAT2_V_shell*MAT2_D*MAT2_V_shell';
% figure
% imagesc(MAT2_shell);

zeroind = MAT1_shell == 0;
MAT1_shell(zeroind) = NaN;
MAT1_shell_lin = nonzeros(triu(MAT1_shell,1));

zeroind = MAT2_shell == 0;
MAT2_shell(zeroind) = NaN;
MAT2_shell_lin = nonzeros(triu(MAT2_shell,1));

corr(MAT1_shell_lin,MAT2_shell_lin,'rows','pairwise')

% MAT1_recMAT1 = MAT1_shell./MAT2_shell;
% inf_ind = isinf(MAT1_recMAT1);
% MAT1_recMAT1(inf_ind) = 0;
% [~,D] = spectral_real(MAT1_recMAT1);
% norm_MAT1_recMAT1 = max(diag(D))
% figure
% imagesc(MAT1_shell);
% figure
% imagesc(MAT1(ind,ind));





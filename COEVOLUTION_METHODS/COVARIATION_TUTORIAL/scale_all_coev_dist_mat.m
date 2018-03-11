points = 1000;
coev_dist_mat = zeros(points,26);
step = 1/points;
step3 = 3/points;

X1 = (1:3*max(X))';
sX1 = (cumsum(max(X1)*ones(points,1)/points));

Y1 = sorted_MI_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,1) = model(sX1);

Y1 = sorted_logR_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,2) = model(sX1);

Y1 = sorted_ZPX2_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,3) = model(sX1);

Y1 = sorted_DCA_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,4) = model(sX1);

Y1 = sorted_nbZPX2_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,5) = model(sX1);

Y1 = sorted_dbZPX2_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,6) = model(sX1);

Y1 = sorted_dgbZPX2_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,7) = model(sX1);

% Y1 = sorted_PSICOV_1(X1,6);
% [model] = shape_pres_fit(X1,Y1);
% coev_dist_mat(:,8) = model(sX1);

% Y1 = sorted_pcZPX2_1(X1,6);
% [model] = shape_pres_fit(X1,Y1);
% coev_dist_mat(:,8) = model(sX1);

Y1 = sorted_SCA_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,9) = model(sX1);

Y1 = sorted_OMES_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,10) = model(sX1);

Y1 = sorted_McBASC_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,11) = model(sX1);

Y1 = sorted_ELSC_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,12) = model(sX1);

Y1 = sorted_fodorSCA_1(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,13) = model(sX1);

Y1 = sorted_MI_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,14) = model(sX1);

Y1 = sorted_logR_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,15) = model(sX1);

Y1 = sorted_ZPX2_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,16) = model(sX1);

Y1 = sorted_DCA_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,17) = model(sX1);

Y1 = sorted_nbZPX2_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,18) = model(sX1);

Y1 = sorted_dbZPX2_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,19) = model(sX1);

Y1 = sorted_dgbZPX2_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,20) = model(sX1);

% Y1 = sorted_PSICOV_2(X1,6);
% [model] = shape_pres_fit(X1,Y1);
% coev_dist_mat(:,21) = model(sX1);

% Y1 = sorted_pcZPX2_2(X1,6);
% [model] = shape_pres_fit(X1,Y1);
% coev_dist_mat(:,21) = model(sX1);

Y1 = sorted_SCA_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,22) = model(sX1);

Y1 = sorted_OMES_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,23) = model(sX1);

Y1 = sorted_McBASC_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,24) = model(sX1);

Y1 = sorted_ELSC_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,25) = model(sX1);

Y1 = sorted_fodorSCA_2(X1,6);
[model] = shape_pres_fit(X1,Y1);
coev_dist_mat(:,26) = model(sX1);

step = step*100;
step3 = step3*100;

coev_dist_mat = [coev_dist_mat sX1 (step:step:100)' (step3:step3:300)' ];


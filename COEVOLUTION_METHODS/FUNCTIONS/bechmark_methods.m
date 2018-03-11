%%
% This script evaluate the CPU execution time of different methods with a
% reference MSA.

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'MI');
t_MI = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'ZPX2');
t_ZPX2 = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'dbZPX2');
t_dbZPX2 = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'dgbZPX2');
t_dgbZPX2 = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'nbZPX2');
t_nbZPX2 = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'DCA');
t_DCA = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'bayesMI');
t_bayesMI = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'OMES');
t_OMES = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'McBASC');
t_McBASC = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'ELSC');
t_ELSC = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'fodorSCA');
t_fodorSCA = toc(tStart);

tStart = tic;
get_nmsa_covar_vec(REF_nmsa,20,'rsemSCA');
t_ramaSCA = toc(tStart);


coev_cputime = [t_MI t_ZPX2 t_nbZPX2 t_dbZPX2 t_dgbZPX2 t_DCA t_bayesMI ...
    t_OMES t_McBASC t_ELSC t_fodorSCA t_ramaSCA]';

n = 12;
%figure; 
h = bar(coev_cputime);
colormap(jet(n));

ch = get(h,'Children');
fvd = get(ch,'Faces');
fvcd = get(ch,'FaceVertexCData');

[zs, izs] = sortrows(coev_cputime,1);

for i = 1:n
    row = izs(i);
    fvcd(fvd(row,:)) = i;
end
set(ch,'FaceVertexCData',fvcd)


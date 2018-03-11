function [ covar_vec_rand ] = ...
    pick_covarions_2( nmsa,ncov,ncov_mult,conserved_pos )
%
% Here we pick the covarying columns. For example, we want 10% of 
% all the 300 positions of the sequence to be covarying with another 10%, 
% which means a total of 60 columns. We will pick two vectors of 30 columns
% (ncov = 30). 
% In practice, since there can be some overlap in the random choice of the
% positions, we may end up with less than 60 columns. For this reason we
% multiply ncov by ncov_mult. A value of 3 for ncov_mult appears more than
% safe.
[nseq,npos] = size(nmsa);
nmsa25 = nmsa == 25;
nmsa25 = sum(nmsa25,1);
nogaps = nmsa25 == 0;
nogaps_ind = find(nogaps);
allowed = setdiff(nogaps_ind,conserved_pos)';
npos = size(allowed,1);

gncov = ncov;
ncov = round(ncov*ncov_mult);
msa_covar_vec = zeros(ncov,2);
ind1 = randi(npos,ncov,1);
ind2 = randi(npos,ncov,1);
msa_covar_vec(:,1) = allowed(ind1);
msa_covar_vec(:,2) = allowed(ind2);

% The following section is necessary to avoid that a column is called to be
% covarying with more than one other column.

unique1 = unique(msa_covar_vec(:,1));

% Here we pick the unique columns numbers that are present in the first 
% selection, but not in the second.
not_shared = setdiff(unique1, msa_covar_vec(:,2));
not_shared_ind = zeros(length(not_shared),1);

for i = 1:length(not_shared)
    not_shared_ind(i,1) = find(msa_covar_vec(:,1) == not_shared(i),1);   
end

msa_covar_vec_unique1 = msa_covar_vec(not_shared_ind,:);

unique2 = unique(msa_covar_vec_unique1(:,2));

% Here we pick the unique columns numbers that are present in the second 
% selection, but not in the first.
not_shared = setdiff(unique2, msa_covar_vec_unique1(:,1));
not_shared_ind = zeros(length(not_shared),1);

for i = 1:length(not_shared)
    not_shared_ind(i,1) = find(msa_covar_vec_unique1(:,2) == not_shared(i),1);   
end

msa_covar_vec_unique2 = msa_covar_vec_unique1(not_shared_ind,:);

% Here we remove any pairs that have conserved positions.
noncons1 = setdiff(msa_covar_vec_unique2(:,1),conserved_pos');
for i = 1:length(noncons1)
nonconsind1(i,1) = find(msa_covar_vec_unique2(:,1) == noncons1(i),1);
end
msa_covar_vec_unique3 = msa_covar_vec_unique2(nonconsind1,:);

noncons2 = setdiff(msa_covar_vec_unique3(:,2),conserved_pos');
for i = 1:length(noncons2)
nonconsind2(i,1) = find(msa_covar_vec_unique3(:,2) == noncons2(i),1);
end
msa_covar_vec_unique4 = msa_covar_vec_unique3(nonconsind2,:);

% Finally we select randomly only a number of covarions equal to the
% desired fraction of total positions.

cols = size(msa_covar_vec_unique4,1);

% Don't change the following or some entries may appear more than once.
% gncov_ind_all = randi(cols,gncov*5,1);
gncov_ind_all = randperm(cols);
% gncov_ind_unique = unique(gncov_ind_all);
gncov_ind = gncov_ind_all(1:gncov);
covar_vec_rand = msa_covar_vec_unique4(gncov_ind,:);

end


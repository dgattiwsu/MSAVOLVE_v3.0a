function [c_distances,sorted_c_distances,sorted_mat1,...
    sorted_mat2,sorted_mat3,sorted_mat4,pdbstruct] = ...
    coev_distance_matrix_3(pdbfile,chain,first_res_no,last_res_no,...
    mat1,mat2,mat3,mat4,radius,near,npreds,plot_dist,fract)
% This is the same as 'coev_distance_matrix_3e'.
% The pdb file must have a header with the correct sequence for this script
% to work. This is a little confusing because different depositors follow 
% different rules to fill the SEQRES field of the header. We suggest that 
% the 'SEQRES' value should be equal to the residue number of the last aa 
% for which there are coordinates. It doesn't matter if there  
% are initial or central gaps in the actual coordinates, because these will 
% correspond to all zero values in the distance matrix. Therefore if you
% have leading or trailing residues (e.g. a HIS-tag) that were not included
% in the msa and therefore in the coevolution matrix, it is crucial to
% specify correctly what are the 1st and last residue numbers in the
% coordinates, which correspond to the msa for which there is a coevolution 
% matrix. If the coordinates start at a residue that comes after the
% beginning of the msa it is better to trim the msa in such a way that
% there is full correspondance between 1st and last residues in the msa and
% in the pdb.
% The bioinformatics toolbox is a requisite.
% chain: 1,2,3,...
% first_res_no: number of the first residue with coordinates in the pdb file
% last-res_no: number of the last residue with coordinates in the pdb file
% plot_dist can have different values:
% 0 : no plot
% 1 : sparsity plot of the distance matrix
% 2 : heat map of the distance matrix with everything beyond the radius zeroed.
% else : like 2, but with the covarions identified by different methods
% overlaid on top.
% radius: threshold distance to select entries in the protein distance
% matrix
% near: minimum separation between residues in sequence to be included in
% the analysis: 1 = consecutive; 2 = separated by 1 intervening residue; 
% 3 = separated by 2 intervening residues; and so on ...
% OUTPUT:
% c_distances: the protein distance matrix based on the centroid of each
% residue
% sorted_mat1...2...3...4: sorted coevolution matrices with 6 columns.
% column 1: coevolution score
% column 2-3: coevolution matrix indices
% column 4: centroid distance of the pair components in angstroms
% column 5: 1 if centroid distance is inside 'radius', 0 if outside
% column 6: cumulative percentage sum of column 5


if nargin == 12
    fract = 1;
end

pdbstruct = pdbread(pdbfile);
modelstruct = pdbstruct.Model(1);

seq = chain;
    % Here we get the IDs for current chain
    chainID =  pdbstruct.Sequence(seq).ChainID;
    
    % Here we get the coords of Atoms and Heteroatoms
    chainCoords = [modelstruct.Atom.chainID]' == chainID;
    
    if isfield(modelstruct, 'HeterogenAtom')
        hetCoords = [modelstruct.HeterogenAtom.chainID]' == chainID;
    % Here we create an Nx3 matrix of 3D coords
        coords = [modelstruct.Atom(chainCoords).X modelstruct.HeterogenAtom(hetCoords).X;
            modelstruct.Atom(chainCoords).Y modelstruct.HeterogenAtom(hetCoords).Y;
            modelstruct.Atom(chainCoords).Z modelstruct.HeterogenAtom(hetCoords).Z; ]';

        order = [modelstruct.Atom(chainCoords).resSeq,...
            modelstruct.HeterogenAtom(hetCoords).resSeq;];

    % Here we reorder the matrix so that Het atoms are in the correct order.
        [order, perm] = sort(order);
        coords = coords(perm,:);
    else
    % Here we create an Nx3 matrix of 3D coords
        coords = [modelstruct.Atom(chainCoords).X;
            modelstruct.Atom(chainCoords).Y;
            modelstruct.Atom(chainCoords).Z; ]';
        
        order = [modelstruct.Atom(chainCoords).resSeq];
    end
            
    % We make sure that we start counting at 1
    if order(1) < 1
        order = order + (1-order(1));
    end
    
    % Here we calculate the pairwise distances. 
    distances = pdist(coords);
        
    % Here we figure out which atoms are within the radius
    closeAtoms = distances<radius;
    [r,c] = find(squareform(closeAtoms));
    
    % This is very tricky: we use 'sparse' to count how many atoms of each 
    % residue are within the radius from the atoms of another residue.
    interactions = sparse(order(r),order(c),1);
    
    % some hetatoms may be bound elements and not part of the actual
    % protein. We will remove them from the plot.    
    residues = pdbstruct.Sequence(seq).NumOfResidues;
    if residues>min(size(interactions))
        interactions(residues,residues) = 0;
    end

interactions = interactions(1:residues,1:residues);
f_interactions = full(interactions);


% Here we determine the centroid of each residue.
res_order_ind = order <= residues;
res_order = order(res_order_ind);
resnum = unique(res_order);
nullnum = setdiff([1:residues],resnum);
numres = numel(resnum);
centroid = zeros(numres,3);
    for i = resnum
        res_ind = order == i;
        centroid(i,:) = nanmean(coords(res_ind,:));
    end 
% And a new distance matrix based on the centroids.
c_distances = squareform(pdist(centroid));
c_distances(nullnum,:) = 0;
c_distances(:,nullnum) = 0;
range = first_res_no:last_res_no;
c_distances = c_distances(range,range);
sorted_c_distances = sort_matrix_descend_2(c_distances,near);
c_distances_ind = c_distances>0 & c_distances<=radius;
sort_c_dist_ind = sorted_c_distances(:,1)>0 & sorted_c_distances(:,1)<=radius;
sort_c_dist_ind2 = sorted_c_distances(sort_c_dist_ind,2:3);

% Now we can apply to each pair in the coevolution matrix the value of that
% pair in the distance matrix.

sorted_mat1 = sort_matrix_descend_2(mat1,near);
for i = 1:size(sorted_mat1,1)
sorted_mat1(i,4) = c_distances(sorted_mat1(i,2),sorted_mat1(i,3));
end
[~,z_ind,~] = intersect(sorted_mat1(:,2:3),sort_c_dist_ind2,'rows');
sorted_mat1(z_ind,5) = 1;
total = sum(sorted_mat1(:,5));
sorted_mat1(:,6) = cumsum(sorted_mat1(:,5))/total;


sorted_mat2 = sort_matrix_descend_2(mat2,near);
for i = 1:size(sorted_mat2,1)
sorted_mat2(i,4) = c_distances(sorted_mat2(i,2),sorted_mat2(i,3));
end
[~,z_ind,~] = intersect(sorted_mat2(:,2:3),sort_c_dist_ind2,'rows');
sorted_mat2(z_ind,5) = 1;
total = sum(sorted_mat2(:,5));
sorted_mat2(:,6) = cumsum(sorted_mat2(:,5))/total;


sorted_mat3 = sort_matrix_descend_2(mat3,near);
for i = 1:size(sorted_mat3,1)
sorted_mat3(i,4) = c_distances(sorted_mat3(i,2),sorted_mat3(i,3));
end
[~,z_ind,~] = intersect(sorted_mat3(:,2:3),sort_c_dist_ind2,'rows');
sorted_mat3(z_ind,5) = 1;
total = sum(sorted_mat3(:,5));
sorted_mat3(:,6) = cumsum(sorted_mat3(:,5))/total;


sorted_mat4 = sort_matrix_descend_2(mat4,near);
for i = 1:size(sorted_mat4,1)
sorted_mat4(i,4) = c_distances(sorted_mat4(i,2),sorted_mat4(i,3));
end
[~,z_ind,~] = intersect(sorted_mat4(:,2:3),sort_c_dist_ind2,'rows');
sorted_mat4(z_ind,5) = 1;
total = sum(sorted_mat4(:,5));
sorted_mat4(:,6) = cumsum(sorted_mat4(:,5))/total;


% Here we plot the distance matrix with superimposed the top scores in the 
% coevolution matrix.

if plot_dist == 0

elseif plot_dist == 1
% Here we visualize the sparsity pattern.    
figure;
spy(interactions(1:residues,1:residues),'.y',8);
set(gca,'Ydir','normal');
title(sprintf('Chain %s: Residues less than %.2f Angstroms apart',...
    chainID,radius));

elseif plot_dist == 2
figure;
imagesc(f_interactions);
set(gca,'Ydir','normal');
% spy(c_distances_ind,'.y',20);
% set(gca,'Ydir','normal');
% figure
% imagesc(c_distances_ind);
% set(gca,'Ydir','normal');

else
% Here we overlay the pairs identified by the coevolution programs.

npreds = round(npreds/fract);
thr_size = size(mat1,1);
thr1 = false(thr_size);
thr2 = false(thr_size);
thr3 = false(thr_size);
thr4 = false(thr_size);
thr1_val = zeros(thr_size);
thr2_val = zeros(thr_size);
thr3_val = zeros(thr_size);
thr4_val = zeros(thr_size);

for n = 1:npreds
    i = sorted_mat1(n,2);
    j = sorted_mat1(n,3);
    thr1(i,j) = true;
    thr1_val(i,j) = sorted_mat1(n,1);
    %thr1(j,i) = true;
    i = sorted_mat2(n,2);
    j = sorted_mat2(n,3);
    thr2(i,j) = true;
    thr2_val(i,j) = sorted_mat2(n,1);    
    %thr2(j,i) = true;
    i = sorted_mat3(n,2);
    j = sorted_mat3(n,3);
    thr3(i,j) = true;
    thr3_val(i,j) = sorted_mat3(n,1);
    %thr3(j,i) = true;
    i = sorted_mat4(n,2);
    j = sorted_mat4(n,3);
    thr4(i,j) = true;
    thr4_val(i,j) = sorted_mat4(n,1);
    %thr4(j,i) = true;
end
    
% cutoff = sorted_mat1(npreds,1);
% thr1 = mat1>=cutoff;
% cutoff = sorted_mat2(npreds,1);
% thr2 = mat2>=cutoff;
% cutoff = sorted_mat3(npreds,1);
% thr3 = mat3>=cutoff;
% cutoff = sorted_mat4(npreds,1);
% thr4 = mat4>=cutoff;


vname = @(x) inputname(1);
name1 = vname(mat1);
name2 = vname(mat2);
name3 = vname(mat3);
name4 = vname(mat4);

X = sort_c_dist_ind2(:,1);
Y = sort_c_dist_ind2(:,2);

CONTACT_MAP_1 = figure; 
    	set(CONTACT_MAP_1,'Units','normalized','Position',[0.3 0.0 0.7 1.0 ],...
    	'Name','CONTACT MAP 1'); clf;
        
ps = 40;
pc = 0.9;
mult = 80;

subplot1 = subplot(2,2,1,'Parent',figure(gcf),...
    'PlotBoxAspectRatio',[residues+1 residues+1 1]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[0 residues+1]);
ylim(subplot1,[0 residues+1]);

% Colormap-----------------------------------------------------------------
S = sorted_mat1(1:npreds,1);
min_S = min(S);
S = S - min_S;
max_S = max(S);
scale = 1/max_S;
S0 = S*scale;
S1 = S0;
S2 = ones(npreds,1)*1;
CM = [S2 S0 S2];
CM = [1 0 1];
%--------------------------------------------------------------------------
plot_dist_map(sorted_mat1,X,Y,npreds,CM,S0,ps,pc,mult);
% freezeColors(gca);

% CONTACT_MAP_2 = figure; 
%     	set(CONTACT_MAP_2,'Units','normalized','Position',[0.6 0.6 0.4 0.4 ],...
%     	'Name','CONTACT MAP 2'); clf;
subplot2 = subplot(2,2,2,'Parent',figure(gcf),...
    'PlotBoxAspectRatio',[residues+1 residues+1 1]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot1,[0 residues+1]);
ylim(subplot1,[0 residues+1]);

% Colormap-----------------------------------------------------------------
S = sorted_mat2(1:npreds,1);
min_S = min(S);
S = S - min_S;
max_S = max(S);
scale = 1/max_S;
S0 = S*scale;
S1 = S0;
S2 = ones(npreds,1)*1;
CM = [S0 S2 S2];
CM = [0 1 1];

%--------------------------------------------------------------------------
plot_dist_map(sorted_mat2,X,Y,npreds,CM,S0,ps,pc,mult);
% freezeColors(gca);

% CONTACT_MAP_3 = figure; 
%     	set(CONTACT_MAP_3,'Units','normalized','Position',[0.2 0.2 0.4 0.4 ],...
%     	'Name','CONTACT MAP 3'); clf;
subplot3 = subplot(2,2,3,'Parent',figure(gcf),...
    'PlotBoxAspectRatio',[residues+1 residues+1 1]);
box(subplot3,'on');
hold(subplot3,'all');
xlim(subplot1,[0 residues+1]);
ylim(subplot1,[0 residues+1]);

% Colormap-----------------------------------------------------------------
S = sorted_mat3(1:npreds,1);
min_S = min(S);
S = S - min_S;
max_S = max(S);
scale = 1/max_S;
S0 = S*scale;
S1 = S0;
S2 = ones(npreds,1)*1;
S3 = S1;
CM = [S0 S2 S0];
CM = [0.3,0.6,0.2];

%--------------------------------------------------------------------------
plot_dist_map(sorted_mat3,X,Y,npreds,CM,S0,ps,pc,mult);
% freezeColors(gca);

% CONTACT_MAP_4 = figure; 
%     	set(CONTACT_MAP_4,'Units','normalized','Position',[0.6 0.2 0.4 0.4 ],...
%     	'Name','CONTACT MAP 4'); clf;
subplot4 = subplot(2,2,4,'Parent',figure(gcf),...
    'PlotBoxAspectRatio',[residues+1 residues+1 1]);
box(subplot4,'on');
hold(subplot4,'all');
xlim(subplot1,[0 residues+1]);
ylim(subplot1,[0 residues+1]);

% Colormap-----------------------------------------------------------------
S = sorted_mat4(1:npreds,1);
min_S = min(S);
S = S - min_S;
max_S = max(S);
scale = 1/max_S;
S0 = S*scale;
S1 = S0;
S2 = ones(npreds,1)*1;
CM = [S1 S2 S1];
CM = [0 1 0];;
%--------------------------------------------------------------------------
plot_dist_map(sorted_mat4,X,Y,npreds,CM,S0,ps,pc,mult);
% freezeColors(gca);

% legend('Dist. matrix',name4,'Location','Best');
% hold off

% spy(thr,'*b',6);

% spy(thr2,'*r',6);
% spy(thr,'*r',6);

% spy(thr3,'sg',6);
% spy(thr,'*g',6);

% spy(thr4,'dm',6);
% spy(thr,'*g',6);


end



end

function plot_dist_map(coords,X,Y,npreds,CM,S0,ps,pc,mult)
%--------------------------------------------------------------------------

% Create plot
scatter(X,Y,ps,'Marker','o','LineWidth',0.1,...
    'MarkerEdgeColor',[pc pc pc],...
    'MarkerFaceColor',[pc pc pc]);hold on
scatter(Y,X,ps,'Marker','o','LineWidth',0.1,...
    'MarkerEdgeColor',[pc pc pc],...    
    'MarkerFaceColor',[pc pc pc]);

X1 = coords(1:npreds,2);
Y1 = coords(1:npreds,3);
C1 = coords(1:npreds,1);
colormap(CM);
% scatter(X1,Y1,S0*100,C1,'filled','Marker','o')
% scatter(Y1,X1,S0*100,C1,'filled','Marker','o')
scatter(X1,Y1,S0*mult,CM,'filled','o')
scatter(Y1,X1,S0*mult,CM,'filled','o')

xlabel('Residue Number');
%--------------------------------------------------------------------------
end

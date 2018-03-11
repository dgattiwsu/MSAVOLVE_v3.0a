function [score,entropy,relentropy,interactions,sfu_interactions_1,...
    sfu_interactions_2,sfu_interactions_3,sfu_interactions_4,pdbstruct] = ...
    coev_distance_matrix(nmsa,smsa,pdbfile,chain,mat1,mat2,...
    mat3,mat4,radius,near,npreds,plot_dist)
% The pdb file must have a header with the correct sequence for this script
% to work. This means that the 'SEQRES' entry of the pdb file must have the 
% same dimensions of the coevolution matrix, and the 'SEQRES' value must be
% equal to the number of aa in the 'SEQRES' field. It doesn't matter if there  
% are flanking or central gaps in the actual coordinates, because these will 
% correspond to all zero values in the interaction matrix. Therefore if you
% have leading or trailing residues (e.g. a HIS-tag) that were not included
% in the msa and therefore in the coevolution matrix, these must be removed
% from both the 'SEQRES' field and the coordinates.
% The bioinformatics toolbox is a requisite.
% chain: 1,2,3,...
% plot_dist can have different values:
% 0 : no plot
% 1 : sparsity plot of the distance matrix
% 2 : heat map of the distance matrix with everything beyond the radius zeroed.
% else : like 2, but with the covarions identified by different methods
% overlaid on top.
% The first column of entropy and relentropy refers to the positions of the
% msa that are within the top 'npred' covarions identified by a method and
% also within the distance matrix radius ('radius'), and the
% closeness cutoff ('near') in sequence.
% The 2nd column of entropy and relentropy refers just to the positions of 
% the msa that are within the top 'npred' covarions identified by a method.

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
    % protein. We will remove them from the plot as they tend to skew it,
    % however these may be of interest.
    residues = pdbstruct.Sequence(seq).NumOfResidues;
    if residues>min(size(interactions))
        interactions(residues,residues) = 0;
    end

interactions = interactions(1:residues,1:residues);
f_interactions = full(interactions);


% Here we calculate entropy and relative entropy for the nmsa.
% Here we calculate the background probabilities at each position.
REF_hmm_model = nmsa_to_hmm_model(nmsa);
bg_frequencies = REF_hmm_model.NullEmission';
% Next we calculate the profile of the MSA without including gaps.
[REF_profile,~] = seqprofile(smsa,'Gaps','none');
% Here we determine the level of conservation using the traditional 
% definitions of entropy and relative entropy.
REF_entropy = Entropy(nmsa);
[~,REF_rel_entropy] = rel_entropy_nmsa(nmsa,...
    bg_frequencies,REF_profile);


sorted_mat1 = sort_matrix_descend(mat1);
cutoff = sorted_mat1(npreds,1);
thr1 = mat1>=cutoff;
sf_interactions = f_interactions;
sf_interactions(~thr1) = 0;
sfu_interactions_1 = triu(sf_interactions,near);
score(1,1) = sum(sfu_interactions_1(:));
[row,col]=find(sfu_interactions_1);
pick_vec = unique([row;col]);
pick_vec_2 = unique([sorted_mat1(1:npreds,2);sorted_mat1(1:npreds,3);]);
entropy(1,1) = mean(REF_entropy(pick_vec));
entropy(1,2) = mean(REF_entropy(pick_vec_2));
relentropy(1,1) = mean(REF_rel_entropy(pick_vec));
relentropy(1,2) = mean(REF_rel_entropy(pick_vec_2));

sorted_mat2 = sort_matrix_descend(mat2);
cutoff = sorted_mat2(npreds,1);
thr2 = mat2>cutoff;
sf_interactions = f_interactions;
sf_interactions(~thr2) = 0;
sfu_interactions_2 = triu(sf_interactions,near);
score(2,1) = sum(sfu_interactions_2(:));
[row,col]=find(sfu_interactions_2);
pick_vec = unique([row;col]);
pick_vec_2 = unique([sorted_mat2(1:npreds,2);sorted_mat2(1:npreds,3);]);
entropy(2,1) = mean(REF_entropy(pick_vec));
entropy(2,2) = mean(REF_entropy(pick_vec_2));
relentropy(2,1) = mean(REF_rel_entropy(pick_vec));
relentropy(2,2) = mean(REF_rel_entropy(pick_vec_2));

sorted_mat3 = sort_matrix_descend(mat3);
cutoff = sorted_mat3(npreds,1);
thr3 = mat3>cutoff;
sf_interactions = f_interactions;
sf_interactions(~thr3) = 0;
sfu_interactions_3 = triu(sf_interactions,near);
score(3,1) = sum(sfu_interactions_3(:));
[row,col]=find(sfu_interactions_3);
pick_vec = unique([row;col]);
pick_vec_2 = unique([sorted_mat3(1:npreds,2);sorted_mat3(1:npreds,3);]);
entropy(3,1) = mean(REF_entropy(pick_vec));
entropy(3,2) = mean(REF_entropy(pick_vec_2));
relentropy(3,1) = mean(REF_rel_entropy(pick_vec));
relentropy(3,2) = mean(REF_rel_entropy(pick_vec_2));

sorted_mat4 = sort_matrix_descend(mat4);
cutoff = sorted_mat4(npreds,1);
thr4 = mat4>cutoff;
sf_interactions = f_interactions;
sf_interactions(~thr4) = 0;
sfu_interactions_4 = triu(sf_interactions,near);
score(4,1) = sum(sfu_interactions_4(:));
[row,col]=find(sfu_interactions_4);
pick_vec = unique([row;col]);
pick_vec_2 = unique([sorted_mat4(1:npreds,2);sorted_mat4(1:npreds,3);]);
entropy(4,1) = mean(REF_entropy(pick_vec));
entropy(4,2) = mean(REF_entropy(pick_vec_2));
relentropy(4,1) = mean(REF_rel_entropy(pick_vec));
relentropy(4,2) = mean(REF_rel_entropy(pick_vec_2));


if plot_dist == 0

elseif plot_dist == 1
% Here we visualize the sparsity pattern.    
figure;
spy(interactions(1:residues,1:residues),'.y',8);
title(sprintf('Chain %s: Residues less than %.2f Angstroms apart',...
    chainID,radius));

elseif plot_dist == 2
figure;
imagesc(f_interactions);set(gca,'Ydir','normal');

else
% Here we overlay the pairs identified by the coevolution programs.

figure;
imagesc(f_interactions);set(gca,'Ydir','normal');

hold on
spy(thr1,'oy',6);
%spy(thr,'*b',6);

spy(thr2,'*r',6);
%spy(thr,'*r',6);

spy(thr3,'sg',6);
%spy(thr,'*g',6);

spy(thr4,'dm',6);
%spy(thr,'*g',6);

vname = @(x) inputname(1);
name1 = vname(mat1);
name2 = vname(mat2);
name3 = vname(mat3);
name4 = vname(mat4);

legend(name1,name2,name3,name4,'Location','BestOutside');

end



end


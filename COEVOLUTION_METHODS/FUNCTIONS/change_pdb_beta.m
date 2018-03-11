function [mod_pdb] = change_pdb_beta(pdb,beta_array,scale)
% This function replaces the temperature factors of a pdb with those of an
% array supplied by the user, which contains the residues numbers in the
% 1st columns and the beta factors in the 2nd column. A scale factor can be
% used to scale the beta's in the array.

mod_pdb = pdb;
natoms = length(mod_pdb.Model.Atom);
nres = size(beta_array,1);
if length(scale) == 1
    scale = ones(nres,1)*scale;
end
for i = 1:natoms
    resnum = mod_pdb.Model.Atom(1,i).resSeq;
    % mod_pdb.Model.Atom(1,i).tempFactor = scale*10^(1-beta_array(resnum,2));
%     mod_pdb.Model.Atom(1,i).tempFactor = ...
%         scale*(max(beta_array(:,2))-beta_array(resnum,2));
    resind = find(beta_array(:,1) == resnum);
    mod_pdb.Model.Atom(1,i).tempFactor = ...
        scale(resind).*beta_array(resind,2);
end
    

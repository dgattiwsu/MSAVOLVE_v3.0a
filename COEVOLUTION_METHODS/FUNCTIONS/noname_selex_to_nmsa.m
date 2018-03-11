function [nmsa,nrows] = noname_selex_to_nmsa(aln_name)
aln_varname = ['s_' aln_name];
import_noname_selex([aln_name '.aln'],aln_varname);
nrows = size(eval(aln_varname),1);
for i = 1:nrows
    nmsa(i,:) = aa2int(eval([aln_varname '{i,:}']));
end
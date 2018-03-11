function [ aln ] = nmsa_to_aln( nmsa,out_file_name,header )

smsa = nmsa_to_smsa(nmsa);

multialignwrite(out_file_name,smsa,...
     'Header',header,'Format','ALN','WriteCount','true');
end
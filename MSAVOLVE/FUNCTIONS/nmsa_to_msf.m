function [ aln ] = nmsa_to_msf( nmsa,out_file_name,header )

smsa = nmsa_to_smsa(nmsa);

multialignwrite(out_file_name,smsa,...
     'Header',header,'Format','MSF','WriteCount','true');
end
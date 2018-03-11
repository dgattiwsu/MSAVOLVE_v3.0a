function [ aln ] = nmsa_to_faln( nmsa,out_file_name )

smsa = nmsa_to_smsa(nmsa);

fastawrite(out_file_name,smsa);

end
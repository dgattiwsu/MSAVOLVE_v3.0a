function [ output_args ] = printselex( nmsa )
% This function converts the sequence from matlab numeric to selex and
% writes it out as an ascii file.
rows = size(nmsa,1);
cmsa_select = int2aa(nmsa);

    jcmsa0 = ['Sequence_' int2str(0) '     '];
    for i = 1:rows
    jcmsa = ['Sequence_' int2str(i)];
    jcmsa0 = char(jcmsa0,jcmsa);
    end
    jcmsa0 = jcmsa0(2:end,:);  
    jcmsa_select = [jcmsa0 cmsa_select];
    
    fid = fopen('msa.selex','w');
    for i = 1:rows
        fprintf(fid, '%s\n', jcmsa_select(i,:));
    end
    fclose(fid);

    clear jcmsa jcmsa0    
end




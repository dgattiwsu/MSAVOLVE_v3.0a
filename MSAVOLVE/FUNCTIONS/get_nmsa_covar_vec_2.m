function [coev_mat,covar_vec] = get_nmsa_covar_vec_2(nmsa,fcov,method)
% This function calculates various coevolution matrices and returns the
% coevolution vector for the fraction of sequence positions requested. It
% differs from 'get_nmsa_covar_vec' because it uses explicitly the function
% NMSA_to_ngcMI, which does not apply a correction for gaps when it 
% calculates simple MI.
npos = size(nmsa,2);
ncov = round(npos*fcov/100);

    switch method

        case 'MI'
            MI = NMSA_to_ngcMI(nmsa);
            sorted_MI = sort_matrix_descend(MI);
            covar_vec = sorted_MI(1:ncov,2:3);            
            coev_mat = MI;

        case 'MIP'
            MI = NMSA_to_ngcMI(nmsa);
            MIP = MI_to_MIP(MI);
            sorted_MIP = sort_matrix_descend(MIP);
            covar_vec = sorted_MIP(1:ncov,2:3);
            coev_mat = MIP;
 
        case 'ZPX'
            MI = NMSA_to_ngcMI(nmsa);
            MIP = MI_to_MIP(MI);
            [ZPX,~] = MI_to_ZPX(MIP);
            sorted_ZPX = sort_matrix_descend(ZPX);
            covar_vec = sorted_ZPX(1:ncov,2:3);
            coev_mat = ZPX;
 
        case 'ZPX2'
            MI = NMSA_to_ngcMI(nmsa);
            MIP = MI_to_MIP(MI);
            [~,ZPX2] = MI_to_ZPX(MIP);
            sorted_ZPX2 = sort_matrix_descend(ZPX2);
            covar_vec = sorted_ZPX2(1:ncov,2:3); 
            coev_mat = ZPX2;
            
        case 'ZRES'
            MI = NMSA_to_ngcMI(nmsa);
            ZRES = MI_to_ZRES(MI);
            sorted_ZRES = sort_matrix_descend(ZRES);
            covar_vec = sorted_ZRES(1:ncov,2:3); 
            coev_mat = ZRES;

        case 'ZRES2' % This is the ZRES matrix calculated applying the 
                     % ZPX algorithm on the MI matrix rather than the MIP
                     % matrix.
            MI = NMSA_to_ngcMI(nmsa);
            [~,ZRES2] = MI_to_ZPX(MI);
            sorted_ZRES2 = sort_matrix_descend(ZRES2);
            covar_vec = sorted_ZRES2(1:ncov,2:3); 
            coev_mat = ZRES2;
                        
        case 'ZMI'
            MI = NMSA_to_ngcMI(nmsa);
            ZMI = MI_to_ZNMI(MI);
            sorted_ZMI = sort_matrix_descend(ZMI);
            covar_vec = sorted_ZMI(1:ncov,2:3); 
            coev_mat = ZMI;

        case 'NMI'
            [~,NMI] = NMSA_to_NMI(MSA_select);
            sorted_NMI = sort_matrix_descend(NMI);
            covar_vec = sorted_NMI(1:ncov,2:3); 
            coev_mat = NMI;                       

        case 'ZNMI'
            [~,NMI] = NMSA_to_NMI(MSA_select);
            ZNMI = MI_to_ZNMI(NMI);
            sorted_ZNMI = sort_matrix_descend(ZNMI);
            covar_vec = sorted_ZNMI(1:ncov,2:3); 
            coev_mat = ZNMI;           
  
        case 'OMES'
        printselex(nmsa);
        !java covariance.algorithms.OmesCovariance msa.selex OMES.matrix
        import_omes_matrix_caller('OMES.matrix');
        data(:,1:2) = data(:,1:2) + 1;
        OMES = zeros(npos,npos);
        for k = 1:size(data,1)
           i = data(k,1);
           j = data(k,2);
               OMES(i,j) = data(k,3);
        end
        OMES = OMES + OMES';
        for i = 1:npos
            OMES(i,i) = NaN;
        end
        sorted_OMES = sort_matrix_descend(OMES);
        covar_vec = sorted_OMES(1:ncov,2:3);
        coev_mat = OMES;

        case 'ELSC'
        printselex(nmsa);
        !java covariance.algorithms.ELSCCovariance msa.selex ELSC.matrix
        import_omes_matrix_caller('ELSC.matrix');
        data(:,1:2) = data(:,1:2) + 1;
        ELSC = zeros(npos,npos);
        for k = 1:size(data,1)
            i = data(k,1);
            j = data(k,2);
                ELSC(i,j) = data(k,3);
        end
        ELSC = ELSC + ELSC';
        for i = 1:npos
            ELSC(i,i) = NaN;
        end      
        sorted_ELSC = sort_matrix_descend(ELSC);
        covar_vec = sorted_ELSC(1:ncov,2:3);
        coev_mat = ELSC;
 
        case 'McBASC'
        printselex(nmsa);
        !java covariance.algorithms.McBASCCovariance msa.selex McBASC.matrix
        import_omes_matrix_caller('McBASC.matrix');
        data(:,1:2) = data(:,1:2) + 1;
        McBASC = zeros(npos,npos);
        for k = 1:size(data,1)
            i = data(k,1);
            j = data(k,2);
                McBASC(i,j) = data(k,3);
        end
        McBASC = McBASC + McBASC';
        for i = 1:npos
             McBASC(i,i) = NaN;
        end
        sorted_McBASC = sort_matrix_descend(McBASC);
        covar_vec = sorted_McBASC(1:ncov,2:3);
        coev_mat = McBASC;
 
        case 'fodorSCA'
        printselex(nmsa);
        !java covariance.algorithms.JavaSCA msa.selex fodorSCA.matrix
        import_omes_matrix_caller('fodorSCA.matrix');
        data(:,1:2) = data(:,1:2) + 1;
        fodorSCA = zeros(npos,npos);
        for k = 1:size(data,1)
           i = data(k,1);
           j = data(k,2);
               fodorSCA(i,j) = data(k,3);
        end
        fodorSCA = fodorSCA + fodorSCA';
        for i = 1:npos
             McBASC(i,i) = NaN;
        end
        sorted_fodorSCA = sort_matrix_descend(fodorSCA);
        covar_vec = sorted_fodorSCA(1:ncov,2:3);
        coev_mat = fodorSCA;
        
        case 'ssemSCA'
        cnmsa = int2aa(nmsa);
        [SSEM_SCA_DDG_mat,~] = SMSA_to_SSEM_SCA(cnmsa);
        [~,~,SSEM_SCA] = get_sca(SSEM_SCA_DDG_mat);
        sorted_SSEM_SCA = sort_matrix_descend(SSEM_SCA);
        covar_vec = sorted_SSEM_SCA(1:ncov,2:3);
        coev_mat = SSEM_SCA;
        
        case 'rsemSCA'        
        cnmsa = int2aa(nmsa);
        [RSEM_SCA_DDG_mat,~] = SMSA_to_RSEM_SCA(cnmsa,1000,0.5);
        [~,~,RSEM_SCA] = get_sca(RSEM_SCA_DDG_mat);
        sorted_RSEM_SCA = sort_matrix_descend(RSEM_SCA);
        covar_vec = sorted_RSEM_SCA(1:ncov,2:3);
        coev_mat = RSEM_SCA;

        case 'DCA'
        % Here we write out the msa in fasta format.   
        nmsa_to_faln(nmsa,'msa.faln');
        % Here we perform DCA using the original code from Andrea Pagnani.
        dca('msa.faln','dca_mat.matrix');
        % Here we read in the dca matrix
        import_dca_matrix('dca_mat.matrix');
        % Here we convert the dca_mat.matrix into our internal matrix
        [ ~,DCA,~,~ ] = dca_to_wMI_DCA( dca_mat,fcov );
        % Here we remove the external files generated at each cycle.
        clear dca_mat
        !rm msa.faln dca_mat.matrix
        sorted_DCA = sort_matrix_descend(DCA);
        covar_vec = sorted_DCA(1:ncov,2:3);            
        coev_mat = DCA;

        case 'bayesMI'
        % First we remove all external files.
        !rm msa.aln msa.faln msa.fasta msa.post postmi_mat.matrix
        !rm msa.mapping msa_ungapped.aln msa_ungapped.post
        % Here we write out the msa in fasta format.
        nmsa_to_faln(nmsa,'msa.faln');
        % Here we perform bayesMI using the original code from Lukas Burger:
        % lukas.burger@fmi.ch
        !t_coffee -other_pg seq_reformat -in msa.faln -output fasta_aln > msa.fasta
        !perl runContactPredictions.pl msa.fasta
        !mv msa.post postmi_mat.matrix 
        % Here we read in the bayesMI matrix
        import_postmi_matrix_caller('postmi_mat.matrix');
        % Here we convert the postmi_mat.matrix into our internal matrix.
        [ bayesMI,~] = postmi_to_bayesMI( postmi_mat,npos,fcov );
        ones_ind = bayesMI == 1;
        bayesMI(ones_ind) = NaN;
        sorted_bayesMI = sort_matrix_descend(bayesMI);
        covar_vec = sorted_bayesMI(1:ncov,2:3);            
        coev_mat = bayesMI;
        % Here we remove the internal and external files.
        clear postmi_mat
        !rm msa.aln msa.faln msa.fasta postmi_mat.matrix
        !rm msa.mapping msa_ungapped.aln msa_ungapped.post
        
        case 'gbZPX2'
        [gbZPX2] = NMSA_to_gbZPX2(nmsa);
        sorted_gbZPX2 = sort_matrix_descend(gbZPX2);
        covar_vec = sorted_gbZPX2(1:ncov,2:3);            
        coev_mat = gbZPX2;
        
        case 'nbZPX2'
        [nbZPX2] = NMSA_to_ngc_nbZPX2(nmsa);
        sorted_nbZPX2 = sort_matrix_descend(nbZPX2);
        covar_vec = sorted_nbZPX2(1:ncov,2:3);            
        coev_mat = nbZPX2;
               
        case 'dbZPX2'
        [dbZPX2] = NMSA_to_dbZPX2(nmsa);
        sorted_dbZPX2 = sort_matrix_descend(dbZPX2);
        covar_vec = sorted_dbZPX2(1:ncov,2:3);            
        coev_mat = dbZPX2;
        
        case 'fgbZPX2'
        [fgbZPX2] = NMSA_to_fgbZPX2(nmsa);
        sorted_fgbZPX2 = sort_matrix_descend(fgbZPX2);
        covar_vec = sorted_fgbZPX2(1:ncov,2:3);            
        coev_mat = fgbZPX2;
                    
        case 'fdgbZPX2'
        [fdgbZPX2] = NMSA_to_fdgbZPX2(nmsa);
        sorted_fdgbZPX2 = sort_matrix_descend(fdgbZPX2);
        covar_vec = sorted_fdgbZPX2(1:ncov,2:3);            
        coev_mat = fdgbZPX2;
        
        case 'dgbZPX2'
        [dgbZPX2] = NMSA_to_dgbZPX2(nmsa);
        sorted_dgbZPX2 = sort_matrix_descend(dgbZPX2);
        covar_vec = sorted_dgbZPX2(1:ncov,2:3);            
        coev_mat = dgbZPX2;
                    
        case 'none'
        covar_vec = [];
        coev_mat = [];
        
 % Clean up external files.
    delete('*.matrix');
    delete('*.selex');

    end

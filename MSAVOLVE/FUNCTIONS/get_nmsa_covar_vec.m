function [coev_mat,covar_vec] = get_nmsa_covar_vec(nmsa,fcov,method,param1)
% This function calculates various coevolution matrices and returns the
% coevolution vector for the fraction of sequence positions requested.
npos = size(nmsa,2);
ncov = round(npos*fcov/100);

    switch method

        case 'MI'
            MI = NMSA_to_fastMI(nmsa);
            sorted_MI = sort_matrix_descend(MI);
            covar_vec = sorted_MI(1:ncov,2:3);            
            coev_mat = MI;

        case 'MIP'
            MI = NMSA_to_fastMI(nmsa);
            MIP = MI_to_MIP(MI);
            sorted_MIP = sort_matrix_descend(MIP);
            covar_vec = sorted_MIP(1:ncov,2:3);
            coev_mat = MIP;
 
        case 'ZPX'
            MI = NMSA_to_fastMI(nmsa);
            MIP = MI_to_MIP(MI);
            [ZPX,~] = MI_to_ZPX(MIP);
            sorted_ZPX = sort_matrix_descend(ZPX);
            covar_vec = sorted_ZPX(1:ncov,2:3);
            coev_mat = ZPX;
 
        case 'ZPX2'
            MI = NMSA_to_fastMI(nmsa);
            MIP = MI_to_MIP(MI);
            [~,ZPX2] = MI_to_ZPX(MIP);
            sorted_ZPX2 = sort_matrix_descend(ZPX2);
            covar_vec = sorted_ZPX2(1:ncov,2:3); 
            coev_mat = ZPX2;
            
        case '3D_ZPX2'
            [~,md3_ZPX2] = ...
            NMSA_to_mdMI(nmsa,'GAPS','3D','FULL',0.9,1,22,0,0,3,12);
            sorted_md3_ZPX2 = sort_matrix_descend(md3_ZPX2);
            covar_vec = sorted_md3_ZPX2(1:ncov,2:3);            
            coev_mat = md3_ZPX2;
            
        case '4D_ZPX2'
            [~,~,md4_ZPX2] = ...
            NMSA_to_mdMI(nmsa,'GAPS','4D','FULL',0.9,1,22,0,0,3,12);
            sorted_md4_ZPX2 = sort_matrix_descend(md4_ZPX2);
            covar_vec = sorted_md4_ZPX2(1:ncov,2:3);            
            coev_mat = md4_ZPX2;
            
        case 'ZRES'
            MI = NMSA_to_fastMI(nmsa);
            ZRES = MI_to_ZRES(MI);
            sorted_ZRES = sort_matrix_descend(ZRES);
            covar_vec = sorted_ZRES(1:ncov,2:3); 
            coev_mat = ZRES;

        case 'ZRES2' % This is the ZRES matrix calculated applying the 
                     % ZPX algorithm on the MI matrix rather than the MIP
                     % matrix.
            MI = NMSA_to_fastMI(nmsa);
            [~,ZRES2] = MI_to_ZPX(MI);
            sorted_ZRES2 = sort_matrix_descend(ZRES2);
            covar_vec = sorted_ZRES2(1:ncov,2:3); 
            coev_mat = ZRES2;
                        
        case 'ZMI'
            MI = NMSA_to_fastMI(nmsa);
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

        case 'ramaSCA'
        % addpath(genpath('SCA_4.5'));
        [RAMA_SCA,~,~] = NMSA_to_ramaSCA(nmsa);
        sorted_RAMA_SCA = sort_matrix_descend(RAMA_SCA);
        covar_vec = sorted_RAMA_SCA(1:ncov,2:3);
        coev_mat = RAMA_SCA;
        % rmpath(genpath('SCA_4.5'));     
              
        case 'DCA'
        % Here we write out the msa in fasta format.
        !rm msa.faln dca_mat.matrix
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
              
        case 'plmDCA_sym'
        % Here we write out the msa in fasta format.
        !rm msa.faln dca_mat.matrix
        nmsa_to_faln(nmsa,'msa.faln');
        plmDCA_symmetric('msa.faln','dca_mat.matrix',0.01,0.01,0.1,12)        
        % Here we read in the dca matrix
        import_dca_matrix('dca_mat.matrix');
        % Here we convert the dca_mat.matrix into our internal matrix
        [ plmDCA,~] = plm_dca_to_wMI_DCA( dca_mat,fcov );
        % Here we remove the external files generated at each cycle.
        clear dca_mat
        !rm msa.faln dca_mat.matrix
        sorted_plmDCA = sort_matrix_descend(plmDCA);
        covar_vec = sorted_plmDCA(1:ncov,2:3);            
        coev_mat = plmDCA;

        case 'plmDCA_asym'
        % Here we write out the msa in fasta format.
        !rm msa.faln dca_mat.matrix
        nmsa_to_faln(nmsa,'msa.faln');
        plmDCA_asymmetric('msa.faln','dca_mat.matrix',0.1,12)        
        % Here we read in the dca matrix
        import_dca_matrix('dca_mat.matrix');
        % Here we convert the dca_mat.matrix into our internal matrix
        [ plmDCA,~] = plm_dca_to_wMI_DCA( dca_mat,fcov );
        % Here we remove the external files generated at each cycle.
        clear dca_mat
        !rm msa.faln dca_mat.matrix
        sorted_plmDCA = sort_matrix_descend(plmDCA);
        covar_vec = sorted_plmDCA(1:ncov,2:3);            
        coev_mat = plmDCA;
               
        case 'gplmDCA_asym'
        % Here we write out the msa in fasta format.
        !rm msa.faln dca_mat.matrix
        nmsa_to_faln(nmsa,'msa.faln');
        gplmDCA_asymmetric('msa.faln','dca_mat.matrix',0.01,0.01,0.001,0.1,12,-1)        
        % Here we read in the dca matrix
        import_dca_matrix('dca_mat.matrix');
        % Here we convert the dca_mat.matrix into our internal matrix
        [ gplmDCA,~] = plm_dca_to_wMI_DCA( dca_mat,fcov );
        % Here we remove the external files generated at each cycle.
        clear dca_mat
        !rm msa.faln dca_mat.matrix
        sorted_gplmDCA = sort_matrix_descend(gplmDCA);
        covar_vec = sorted_gplmDCA(1:ncov,2:3);            
        coev_mat = gplmDCA;
               
        case 'hpPCA'
        % Here we write out the msa in fasta format.
        !rm msa.faln dca_mat.matrix
        nmsa_to_faln(nmsa,'msa.faln');
        % Here we perform hopfield-pottsPCA using the original code from Martin Weigt.
        [Lambda,Vtilde,N,q] = inverse_hopfield_potts1('msa.faln',0.2,0.5);
        lambda_diag = real(diag(Lambda));
        nlambda = length(lambda_diag);
        p = nlambda-sum(lambda_diag <= 1+0.2 & lambda_diag >= 1-0.2);
        %
        % where input.fasta is the input multiple-sequence alignment.
        % 0.2 is a typical value for the sequence distance used for reweighting
        % (for more than 20,000 sequences in the MSA, no reweighting is used).
        % 0.5 is a typical value for the relative weight of the pseudocount.
        %
        if exist('param1','var')
            p = param1;
        end
        fprintf('Using %i patterns. \n', p);
        F_apc = inverse_hopfield_potts2(Vtilde,Lambda,N,q,p);
        %
        % where p has to be replaced by the wanted number of patterns.
        % The output matrix F_apc is the score matrix for all position pairs
        %
        for i = 1:npos
            F_apc(i,i) = NaN;
        end
        !rm msa.faln 
        sorted_F_apc = sort_matrix_descend(F_apc);
        covar_vec = sorted_F_apc(1:ncov,2:3);            
        coev_mat = F_apc;
        
        case 'logR'
        % First we remove all external files.
        !rm msa.aln msa.faln msa.fasta msa.post postmi_mat.matrix
        !rm msa.mapping msa_ungapped.aln msa_ungapped.post
        % Here we write out the msa in fasta format.
        nmsa_to_faln(nmsa,'msa.faln');
        % Here we perform logR using the original code from Lukas Burger:
        % lukas.burger@fmi.ch
        !t_coffee -other_pg seq_reformat -in msa.faln -output fasta_aln > msa.fasta
        !perl runContactPredictions.pl msa.fasta
        !mv msa.post postmi_mat.matrix 
        % Here we read in the logR matrix
        import_postmi_matrix_caller('postmi_mat.matrix');
        % Here we convert the postmi_mat.matrix into our internal matrix.
        [ logR,~] = postmi_to_logR( postmi_mat,npos,fcov );
        ones_ind = logR == 1;
        logR(ones_ind) = NaN;
        sorted_logR = sort_matrix_descend(logR);
        covar_vec = sorted_logR(1:ncov,2:3);            
        coev_mat = logR;
        % Here we remove the internal and external files.
        clear postmi_mat
        !rm msa.aln msa.faln msa.fasta postmi_mat.matrix
        !rm msa.mapping msa_ungapped.aln msa_ungapped.post
        
        case 'GREMLIN'
        [rows,cols] = size(nmsa);
        c_MSA_select = int2aa(nmsa);          
        fid = fopen('MSA_select.selex','w');
        for i = 1:rows
            fprintf(fid, '%s\n', c_MSA_select(i,:));
        end
        fclose(fid);        
        gremlin('MSA_select.selex','gremlin.matrix')            
        GREMLIN = importdata('gremlin.matrix');
        for i = 1:npos
            GREMLIN(i,i) = NaN;
        end
        sorted_GREMLIN = sort_matrix_descend(GREMLIN);
        covar_vec = sorted_GREMLIN(1:ncov,2:3);            
        coev_mat = GREMLIN;
        
        case 'PSICOV'
        rows = size(nmsa,1);
        c_MSA_select = int2aa(nmsa);          
        fid = fopen('MSA_select.selex','w');
        for i = 1:rows
            fprintf(fid, '%s\n', c_MSA_select(i,:));
        end
        fclose(fid);        
        !./psicov -p -d 0.03 -j 0 MSA_select.selex > psicov.matrix
        % !./psicov -d 0.03 -r 0.005 -i 62 MSA_select.selex > psicov.matrix
        % We run here in some relaxed condition to favor convergence        
        % !./psicov -a -r 0.005 -i 62 MSA_select.selex > psicov.matrix
        import_psicov_matrix('psicov.matrix');
        PSICOV = zeros(npos,npos);
        for k = 1:size(psicov,1)
           i = psicov(k,1);
           j = psicov(k,2);
           PSICOV(i,j) = psicov(k,5);
        end
        PSICOV = PSICOV + PSICOV';
        clear psicov
        for i = 1:npos
            PSICOV(i,i) = NaN;
        end
        sorted_PSICOV = sort_matrix_descend(PSICOV);
        covar_vec = sorted_PSICOV(1:ncov,2:3);            
        coev_mat = PSICOV;
        !rm MSA_select.selex psicov.matrix

        case 'localPSICOV'
        [localPSICOV] = NMSA_to_PSICOV(nmsa,...
            'NOGAPS',0.8,1.0,'ZERO','FRO',...
            0.006,0.0,0.0001,'QUIC','RHO',20);
        sorted_localPSICOV = sort_matrix_descend(localPSICOV);
        covar_vec = sorted_localPSICOV(1:ncov,2:3);            
        coev_mat = localPSICOV;
       
        case 'localPSICOV_2'
        [localPSICOV_2] = NMSA_to_PSICOV_2(nmsa,...
            'NOGAPS',0.8,1.0,21,20,'ZERO','KEEP','FRO',...
            0.006,0.0,0.0001,'QUIC','RHO');
        sorted_localPSICOV_2 = sort_matrix_descend(localPSICOV_2);
        covar_vec = sorted_localPSICOV_2(1:ncov,2:3);            
        coev_mat = localPSICOV_2;

        case 'slPSICOV'        
        [slPSICOV] = NMSA_to_slPSICOV(nmsa,'NOGAPS',0.9,1.0,...
            20,'SUM','NONE','fro',0.0001,0.35,0.015,...
            100,0.0,0.0001,'RHO',0,1);        
        sorted_slPSICOV = sort_matrix_descend(slPSICOV);
        covar_vec = sorted_slPSICOV(1:ncov,2:3);            
        coev_mat = slPSICOV;
       
        case 'nbZPX2'
        [nbZPX2] = NMSA_to_nbZPX2(nmsa);
        sorted_nbZPX2 = sort_matrix_descend(nbZPX2);
        covar_vec = sorted_nbZPX2(1:ncov,2:3);            
        coev_mat = nbZPX2;

        case 'nb2_ZPX2'
        [nb2_ZPX2] = NMSA_to_nb2_ZPX2(nmsa,...
            'GAPS','DISTANCE','INCLUDE');
        sorted_nb2_ZPX2 = sort_matrix_descend(nb2_ZPX2);
        covar_vec = sorted_nb2_ZPX2(1:ncov,2:3);            
        coev_mat = nb2_ZPX2;

        case 'nb3_ZPX2'
        [nb3_ZPX2] = NMSA_to_nb3_ZPX2(nmsa,...
            'GAPS','DISTANCE',1.0,0.5,21,'INCLUDE');
        sorted_nb3_ZPX2 = sort_matrix_descend(nb3_ZPX2);
        covar_vec = sorted_nb3_ZPX2(1:ncov,2:3);            
        coev_mat = nb3_ZPX2;

        case 'nb4_ZPX2' % nb4 is the same as nbZPX2_original
        [nb4_ZPX2] = NMSA_to_nb4_ZPX2(nmsa);
        sorted_nb4_ZPX2 = sort_matrix_descend(nb4_ZPX2);
        covar_vec = sorted_nb4_ZPX2(1:ncov,2:3);            
        coev_mat = nb4_ZPX2;

        case 'fnb_ZPX2' % like nb, but with the 'fastMI' function
        [fnb_ZPX2] = NMSA_to_fnb_ZPX2(nmsa,'GAPS',0.9,1,0,0,3);
        sorted_fnb_ZPX2 = sort_matrix_descend(fnb_ZPX2);
        covar_vec = sorted_fnb_ZPX2(1:ncov,2:3);            
        coev_mat = fnb_ZPX2;
               
        case 'dbZPX2'
        [dbZPX2] = NMSA_to_dbZPX2(nmsa);
        sorted_dbZPX2 = sort_matrix_descend(dbZPX2);
        covar_vec = sorted_dbZPX2(1:ncov,2:3);            
        coev_mat = dbZPX2;
        
        case 'db2_ZPX2'
        [db2_ZPX2] = NMSA_to_db2_ZPX2(nmsa,'NOGAPS','DISTANCE','EXCLUDE');
        sorted_db2_ZPX2 = sort_matrix_descend(db2_ZPX2);
        covar_vec = sorted_db2_ZPX2(1:ncov,2:3);            
        coev_mat = db2_ZPX2;
                        
        case 'dgbZPX2'
        [dgbZPX2] = NMSA_to_dgbZPX2(nmsa);
        sorted_dgbZPX2 = sort_matrix_descend(dgbZPX2);
        covar_vec = sorted_dgbZPX2(1:ncov,2:3);            
        coev_mat = dgbZPX2;
        
        case 'dgb2_ZPX2'
        [dgb2_ZPX2] = ...
            NMSA_to_dgb2_ZPX2(nmsa,'GAPS','DISTANCE','INCLUDE','ZERO','FRO');
        sorted_dgb2_ZPX2 = sort_matrix_descend(dgb2_ZPX2);
        covar_vec = sorted_dgb2_ZPX2(1:ncov,2:3);            
        coev_mat = dgb2_ZPX2;
               
        case 'swPSIMI'
        [swPSIMI] = NMSA_to_swPSIMI(nmsa,'NOGAPS',...
            1.0,0.5,21,0.025,0.0,0.0001,'RHO');
        sorted_swPSIMI = sort_matrix_descend(swPSIMI);
        covar_vec = sorted_swPSIMI(1:ncov,2:3);            
        coev_mat = swPSIMI;
        
        case 'fastPSIMI'
        [fastPSIMI] = NMSA_to_fastPSIMI(nmsa,0.025,0.0,0.0001,'RHO');
        sorted_fastPSIMI = sort_matrix_descend(fastPSIMI);
        covar_vec = sorted_fastPSIMI(1:ncov,2:3);            
        coev_mat = fastPSIMI;
                   
        case 'fpcZPX2'
        [fpcZPX2] = NMSA_to_fpcZPX2(nmsa,'NOGAPS',0.90,...
            'SDP',1,'AUTO','QUIC',0.005,...
            0.0,0.001,'INVERSE',21,'NONE',0.05);
        sorted_fpcZPX2 = sort_matrix_descend(fpcZPX2);
        covar_vec = sorted_fpcZPX2(1:ncov,2:3);            
        coev_mat = fpcZPX2;
         
        case 'none'
        covar_vec = [];
        coev_mat = [];
        
 % Clean up external files.
    delete('*.matrix');
    delete('*.selex');

    end

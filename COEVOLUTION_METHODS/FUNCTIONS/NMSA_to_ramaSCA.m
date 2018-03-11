function [ RAMA_SCA,RAMA_simpleSCA,RAMA_eigSCA ] = NMSA_to_ramaSCA( nmsa )
% This function calculates SCA coevolution matrices using the SCA v4.5
% Matlab Toolbox from Rama Ranganathan.


[nrows,ncols] = size(nmsa);

algn = int2aa(nmsa);

%%%%%%%% original code %%%%%%%%%

[C_sca_sn]=sca_mat(algn,'spectral');

%%%%%% end original code %%%%%%%
      
% As described in Halabi et al. (Cell (2009), 138, 774-786), the true
% evolutionary correlations between sequence positions (the "signals") are
% mixed with correaltions arising from two sources of noise: (1) finite
% sampling of sequences (so-called "statistical noise"), and (2) biased
% sampling of sequences due to pure phylogenetic releationships between
% them ("historical noise").  An important goal is to separate signal from
% noise in the SCA matrix.  
% To correct for noise due to finite sampling, we compute the SCA
% correlation matrix for many trials of randomized alignments in which each
% column in the sequence alignment is independently scrambled.  This
% manipulation preserves the composition of amino acids are each position
% but removes all non-random correlation between positions. Thus, the
% average correlation scores for each pair of positions arising from these
% ramdomized alignments provides an estimate of the correlations expected
% just due to finite sampling alone and are subtracted from the SCA matrix.
     
[N_seq,N_pos] = size(algn);

%%%%%%%% original code %%%%%%%%%

N_samples=100;
C_rnd=zeros(N_pos,N_pos);
C_tmp=zeros(N_pos,N_pos,N_samples);

% fprintf('%s\n',' ');
% disp(['Computing random trials...']);

for s=1:N_samples
%    tic;
    for pos=1:N_pos
        perm_seq=randperm(N_seq); 
        algn_rnd(:,pos)=algn(perm_seq(:),pos);
    end
    [C_tmp(:,:,s)]=sca_mat(algn_rnd,'spectral','nv'); 
    C_rnd=C_rnd+squeeze(C_tmp(:,:,s));
%    disp([num2str(s) '    ' num2str(toc) ' seconds']);
%    exec_time(s)=toc;
end
C_rnd=C_rnd./N_samples;

% fprintf('%s\n',' ');
% disp(['Total time    ' num2str(sum(exec_time)) ' seconds']);

% We redefine the final SCA matrix for the actual alignment as the
% noise-subtracted version of C_sca_sn, and also redefine each of the
% random SCA matrices by the same subtraction for further comparisons
% below.

C_sca_final=C_sca_sn-C_rnd;

for s=1:N_samples
    C_tmp(:,:,s)=C_tmp(:,:,s)-C_rnd;
end

%%%%%% end original code %%%%%%%

% Before leaving SCA 4.5 we perform a spectral analysis of the coevolution 
% matrix.
    
%%%%%%%% original code %%%%%%%%%

% As described in Halabi et al. (Cell (2009), 138, 774-786), a more
% powerful approach for sector identification is spectral (or eigenvalue)
% decomposition of the SCA matrix.  In this approach, we first determine
% the number of significant eigenmodes through comparison with randomized
% alignments and then examine the pattern of positional weights along the
% significant eigenvectors to define sectors. Sectors are evident as
% groupings of sequence positions in these top eigenmodes - clusters of
% positions that significantly co-evolve.

[ev,lbd]=eigenvect(C_sca_final);

% Comparison of eigenvalue spectra of the real alignment (bars) with that
% expected from the randomized trials (red line):

N_ev=5;
lbd_rnd=zeros(N_samples,N_pos);
ev_rnd=zeros(N_samples,N_pos,N_ev);
lbd_rnd_tmp=zeros(N_pos,1);
ev_rnd_tmp=zeros(N_pos,N_pos);
for s=1:N_samples
    [ev_rnd_tmp,lbd_rnd_tmp]=eigenvect(C_tmp(:,:,s));  
    lbd_rnd(s,:)=lbd_rnd_tmp;   
   ev_rnd(s,:,:)=ev_rnd_tmp(:,1:N_ev);
end
    
% h_eigenspectrum=figure; 
% set(h_eigenspectrum,'Name','Eigenmodes of SCA matrix'); clf; 
[yhist,xhist]=hist(lbd,N_pos); 
% bar(xhist,yhist,'b'); 
% axis([min(lbd) 1.1*max(lbd) 0 1.4*max(yhist)]);
% hold on
% xlabel('eigenvalue','FontSize',10,'FontWeight','bold');
% ylabel('number','FontSize',10,'FontWeight','bold');
[nn]=hist(lbd_rnd(:),xhist);
% plot(xhist,nn/N_samples,'r','LineWidth',1.5); 
% axis([min(lbd) 1.1*max(lbd) 0 1.4*max(yhist)]);
% grid on;
% hold off

%%%%%% end original code %%%%%%%

% Based on the spectral analysis, the first "neigmod' eigenmodes of the SCA   
% matrix are not-explainable by finite size effects and are used to obtain  
% an eigenfiltered SCA matrix. We take a 'conservative' approach to 
% calculating neigmod, with a minimum value of 5.

nn_ind = find(nn == 0 & xhist > 0);
neigmod = sum(yhist(nn_ind)) - 3;
neigmod_ALL(n,1) = neigmod;
if neigmod<5
    neigmod = 5;
end    
size_lbd = size(lbd,1);
lbd_mat = zeros(size_lbd);

for i = 1:size_lbd
lbd_mat(i,i) = lbd(i,1);
end

eig_C_sca = ev(:,1:neigmod)*lbd_mat(1:neigmod,1:neigmod)*(ev(:,1:neigmod))';

% Here we clear up some memory.

clear algn algn_rnd neigmod
clear ev lbd ev_rnd lbd_rnd ev_rnd_tmp lbd_rnd_tmp lbd_mat 
clear N_ev N_samples N_pos N_seq C_rnd C_tmp perm_seq

% Here we store the various matrices obtained with the SCA 4.5 package in
% the usual way for consistency with the other algorithms.

    RAMA_simpleSCA = C_sca_sn;
    RAMA_SCA = C_sca_final; 
    RAMA_eigSCA = eig_C_sca; 
    for i = 1:npos;
        RAMA_simpleSCA(i,i) = NaN;
        RAMA_SCA(i,i) = NaN;
        RAMA_eigSCA(i,i) = NaN;
    end

end


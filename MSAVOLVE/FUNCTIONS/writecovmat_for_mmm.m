% This script converts coevolution matrices into a format that can be used
% to compare the matrices using Tillier's MMM algorithm.

COV = COV_ALL(:,:,n);
cov_COV = cov_COV_ALL(:,:,n);
mut_COV = mut_COV_ALL(:,:,n);
recomb_COV = recomb_COV_ALL(:,:,n);
MI = MI_ALL(:,:,n);
NMI = NMI_ALL(:,:,n);
ZMI = ZMI_ALL(:,:,n);
MIP = MIP_ALL(:,:,n);
ZPX = ZPX_ALL(:,:,n);
ZPX2 = ZPX2_ALL(:,:,n);
ZRES2 = ZRES2_ALL(:,:,n);
OMES = OMES_ALL(:,:,n);
ELSC = ELSC_ALL(:,:,n);
McBASC = McBASC_ALL(:,:,n);
omesSCA = omesSCA_ALL(:,:,n);
RAMA_SCA = RAMA_SCA_ALL(:,:,n);
RSEM_SCA = RSEM_SCA_ALL(:,:,n);
SSEM_SCA = SSEM_SCA_ALL(:,:,n);


neg_ind = COV<0;
COV_pos = COV;
COV_pos(neg_ind) = 0;
COV_pos = nantozero(COV_pos);
COV_pos = COV_pos/max(COV_pos(:));
% COV_pos = 1-nantozero(COV_pos);
[nrows,ncols] = size(COV);
fid = fopen('COV.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,COV_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            COV_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind COV_pos nrows ncols

neg_ind = recomb_COV<0;
recomb_COV_pos = recomb_COV;
recomb_COV_pos(neg_ind) = 0;
recomb_COV_pos = nantozero(recomb_COV_pos);
recomb_COV_pos = recomb_COV_pos/max(recomb_COV_pos(:));
% recomb_COV_pos = 1-nantozero(recomb_COV_pos);
[nrows,ncols] = size(recomb_COV);
fid = fopen('recomb_COV.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,recomb_COV_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            recomb_COV_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind recomb_COV_pos nrows ncols

neg_ind = cov_COV<0;
cov_COV_pos = cov_COV;
cov_COV_pos(neg_ind) = 0;
cov_COV_pos = nantozero(cov_COV_pos);
cov_COV_pos = cov_COV_pos/max(cov_COV_pos(:));
% cov_COV_pos = 1-nantozero(cov_COV_pos);
[nrows,ncols] = size(cov_COV);
fid = fopen('cov_COV.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,cov_COV_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            cov_COV_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind cov_COV_pos nrows ncols

neg_ind = mut_COV<0;
mut_COV_pos = mut_COV;
mut_COV_pos(neg_ind) = 0;
mut_COV_pos = nantozero(mut_COV_pos);
mut_COV_pos = mut_COV_pos/max(mut_COV_pos(:));
% mut_COV_pos = 1-nantozero(mut_COV_pos);
[nrows,ncols] = size(mut_COV);
fid = fopen('mut_COV.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,mut_COV_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            mut_COV_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind mut_COV_pos nrows ncols

neg_ind = MI<0;
MI_pos = MI;
MI_pos(neg_ind) = 0;
MI_pos = nantozero(MI_pos);
MI_pos = MI_pos/max(MI_pos(:));
% MI_pos = 1-nantozero(MI_pos);
[nrows,ncols] = size(MI);
fid = fopen('MI.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,MI_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            MI_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind MI_pos nrows ncols

neg_ind = NMI<0;
NMI_pos = NMI;
NMI_pos(neg_ind) = 0;
NMI_pos = nantozero(NMI_pos);
NMI_pos = NMI_pos/max(NMI_pos(:));
% NMI_pos = 1-nantozero(NMI_pos);
[rows,cols] = size(NMI);
fid = fopen('NMI.matrix','w');
fprintf(fid,'%5d\n',rows);
    for i = 1:rows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,NMI_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            NMI_pos(i,6:cols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind NMI_pos nrows ncols

neg_ind = ZMI<0;
ZMI_pos = ZMI;
ZMI_pos(neg_ind) = 0;
ZMI_pos = nantozero(ZMI_pos);
ZMI_pos = ZMI_pos/max(ZMI_pos(:));
% ZMI_pos = 1-nantozero(ZMI_pos);
[nrows,ncols] = size(ZMI);
fid = fopen('ZMI.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,ZMI_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            ZMI_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind ZMI_pos nrows ncols

neg_ind = MIP<0;
MIP_pos = MIP;
MIP_pos(neg_ind) = 0;
MIP_pos = nantozero(MIP_pos);
MIP_pos = MIP_pos/max(MIP_pos(:));
% MIP_pos = 1-nantozero(MIP_pos);
[nrows,ncols] = size(MIP);
fid = fopen('MIP.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,MIP_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            MIP_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind MIP_pos nrows ncols

neg_ind = ZPX<0;
ZPX_pos = ZPX;
ZPX_pos(neg_ind) = 0;
ZPX_pos = nantozero(ZPX_pos);
ZPX_pos = ZPX_pos/max(ZPX_pos(:));
% ZPX_pos = 1-nantozero(ZPX_pos);
[nrows,ncols] = size(ZPX);
fid = fopen('ZPX.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,ZPX_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            ZPX_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind ZPX_pos nrows ncols

neg_ind = ZPX2<0;
ZPX2_pos = ZPX2;
ZPX2_pos(neg_ind) = 0;
ZPX2_pos = nantozero(ZPX2_pos);
ZPX2_pos = ZPX2_pos/max(ZPX2_pos(:));
% ZPX2_pos = 1-nantozero(ZPX2_pos);
[nrows,ncols] = size(ZPX2);
fid = fopen('ZPX2.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,ZPX2_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            ZPX2_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind ZPX2_pos nrows ncols

neg_ind = ZRES2<0;
ZRES2_pos = ZRES2;
ZRES2_pos(neg_ind) = 0;
ZRES2_pos = nantozero(ZRES2_pos);
ZRES2_pos = ZRES2_pos/max(ZRES2_pos(:));
% ZRES2_pos = 1-nantozero(ZRES2_pos);
[nrows,ncols] = size(ZRES2);
fid = fopen('ZRES2.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,ZRES2_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            ZRES2_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind ZRES2_pos nrows ncols

neg_ind = OMES<0;
OMES_pos = OMES;
OMES_pos(neg_ind) = 0;
OMES_pos = nantozero(OMES_pos);
OMES_pos = OMES_pos/max(OMES_pos(:));
% OMES_pos = 1-nantozero(OMES_pos);
[nrows,ncols] = size(OMES);
fid = fopen('OMES.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,OMES_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            OMES_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind OMES_pos nrows ncols

neg_ind = ELSC<0;
ELSC_pos = ELSC;
ELSC_pos(neg_ind) = 0;
ELSC_pos = nantozero(ELSC_pos);
ELSC_pos = ELSC_pos/max(ELSC_pos(:));
% ELSC_pos = 1-nantozero(ELSC_pos);
[nrows,ncols] = size(ELSC);
fid = fopen('ELSC.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,ELSC_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            ELSC_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind ELSC_pos nrows ncols

neg_ind = McBASC<0;
McBASC_pos = McBASC;
McBASC_pos(neg_ind) = 0;
McBASC_pos = nantozero(McBASC_pos);
McBASC_pos = McBASC_pos/max(McBASC_pos(:));
% McBASC_pos = 1-nantozero(McBASC_pos);
[nrows,ncols] = size(McBASC);
fid = fopen('McBASC.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,McBASC_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            McBASC_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind McBASC_pos nrows ncols

neg_ind = omesSCA<0;
omesSCA_pos = omesSCA;
omesSCA_pos(neg_ind) = 0;
omesSCA_pos = nantozero(omesSCA_pos);
omesSCA_pos = omesSCA_pos/max(omesSCA_pos(:));
% omesSCA_pos = 1-nantozero(omesSCA_pos);
[nrows,ncols] = size(omesSCA);
fid = fopen('omesSCA.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,omesSCA_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            omesSCA_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind omesSCA_pos nrows ncols

neg_ind = RAMA_SCA<0;
RAMA_SCA_pos = RAMA_SCA;
RAMA_SCA_pos(neg_ind) = 0;
RAMA_SCA_pos = nantozero(RAMA_SCA_pos);
RAMA_SCA_pos = RAMA_SCA_pos/max(RAMA_SCA_pos(:));
% RAMA_SCA_pos = 1-nantozero(RAMA_SCA_pos);
[nrows,ncols] = size(RAMA_SCA);
fid = fopen('RAMA_SCA.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,RAMA_SCA_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            RAMA_SCA_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind RAMA_SCA_pos nrows ncols

neg_ind = RSEM_SCA<0;
RSEM_SCA_pos = RSEM_SCA;
RSEM_SCA_pos(neg_ind) = 0;
RSEM_SCA_pos = nantozero(RSEM_SCA_pos);
RSEM_SCA_pos = RSEM_SCA_pos/max(RSEM_SCA_pos(:));
% RSEM_SCA_pos = 1-nantozero(RSEM_SCA_pos);
[nrows,ncols] = size(RSEM_SCA);
fid = fopen('RSEM_SCA.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,RSEM_SCA_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            RSEM_SCA_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind RSEM_SCA_pos nrows ncols

neg_ind = SSEM_SCA<0;
SSEM_SCA_pos = SSEM_SCA;
SSEM_SCA_pos(neg_ind) = 0;
SSEM_SCA_pos = nantozero(SSEM_SCA_pos);
SSEM_SCA_pos = SSEM_SCA_pos/max(SSEM_SCA_pos(:));
% SSEM_SCA_pos = 1-nantozero(SSEM_SCA_pos);
[nrows,ncols] = size(SSEM_SCA);
fid = fopen('SSEM_SCA.matrix','w');
fprintf(fid,'%5d\n',nrows);
    for i = 1:nrows    
    fprintf(fid,'%05d|%03d|%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
        i,i,SSEM_SCA_pos(i,1:5));
        fprintf(fid, '%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n',...
            SSEM_SCA_pos(i,6:ncols));
        fprintf(fid, '\n');
    end
fclose(fid);
clear neg_ind SSEM_SCA_pos nrows ncols

% neg_ind = MI<0;
% MI_pos = MI;
% MI_pos(neg_ind) = 0;
% dlmwrite('MI_standard.matrix',nantozero(MI_pos), 'delimiter', '', ...
%          'precision', '%10.6f');
% neg_ind = NMI<0;
% NMI_pos = NMI;
% NMI_pos(neg_ind) = 0;     
% dlmwrite('NMI_standard.matrix',nantozero(NMI_pos), 'delimiter', '', ...
%          'precision', '%10.6f');
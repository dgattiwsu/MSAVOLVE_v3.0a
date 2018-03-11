% This script corrects the effect of gaps on the coevolution 
% matrices calculated by various methods. Three progressively stronger
% corrections (gapW,gapW2,gapW3) can be applied. MI, bayesMI, and SCA do
% not need a correction.
%%
    gapW = correct_coevmat_forgaps(nmsa);
    gapW2 = gapW.^2;
    gapW3 = gapW.^3;
    gW = gapW3;

 MI_orig = MI;
 MIP_orig = MIP;
 ZRES_orig = ZRES;
 ZPX2_orig = ZPX2;
 nbZPX2_orig = nbZPX2;
 gbZPX2_orig = gbZPX2;
 dbZPX2_orig = dbZPX2;
 fgbZPX2_orig = fgbZPX2;
 dgbZPX2_orig = dgbZPX2;
 DCA_orig = DCA;
 
 MI = MI_orig.*gW; 
 MIP = MIP_orig.*gW;
 ZRES = ZRES_orig.*gW;
 ZPX2 = ZPX2_orig.*gW;
 nbZPX2 = nbZPX2_orig.*gW;
 gbZPX2 = gbZPX2_orig.*gW;
 dbZPX2 = dbZPX2_orig.*gW;
 fgbZPX2 = fgbZPX2_orig.*gW;
 dgbZPX2 = dgbZPX2_orig.*gW;
 DCA = DCA_orig.*gW;

% contents.m
% 
%--------------------------------------------------------------------------
% Directory MSAVOLVE:
% 
%     There should never be a need to change anything in this directory. It 
%     contains 3 subdirectories:
%     
%     Subdirectory MSAVOLVE:  it contains the simulation program 
%                             MSAvolve_v3_0.m and the script RUN_MSAvolve.m 
%                             that drives it. Earlier versions of the
%                             program are also stored here.
%                             
%     Subdirectory FUNCTIONS: it contains all the functions used by the 
%                             simulation program MSAvolve_v3_0.m.
%                             
%     Subdirectory PLOTS:     it contains the script PLOT_VARIOUS_EVOLVED.m,
%                             which gives some examples of how to plot the 
%                             results obtained with MSAvolve_v3_0.m.
%                             
%--------------------------------------------------------------------------
%  Directory COEVOLUTION_METHODS:
%  
%     It contains all the coevolution detection methods. There are 3 
%     subdirectories:
%     
%     Subdirectory FUNCTIONS: it contains all the functions used by different
%                             coevolution detection methods. Some of these 
%                             functions require additional functions that are 
%                             part of separate toolboxes. These toolboxes 
%                             appear as separate subdirectories with their 
%                             own licenses:
%                             
%                             DCA: the DCA method.
%                             plmDCA: the plmDCA method.
%                             gplmDCA: variation of plmDCA method with
%                                      special treatment of gaps
%                             Hopfield_Potts_PCA: the hpPCA method.
%                             GREMLIN: the GREMLIN method
%                             FastICA: Independent Component Analysis (ICA).
%                             Graphical_Lasso: the GLASSO sparse inverse matrix
%                                              algorithm.
%                             logR: the logR method. This directory is empty as
%                                   we do not have permission to redistribute
%                                   the logR programs.
%                             OMES: Fodor's package of coevolution detection 
%                                   methods.
%                             PSICOV: Jones' psicov method.
%                             QUIC_MEX: the QUIC sparse inverse algorithm.
%                             RPCA: the Robust Principal Component Analysis 
%                                   toolbox.
%                             SCA_Sharma: this directory is empty because we 
%                                         do not have permission to redistribute
%                                         this toolbox.
%                             SCA_4.5: this directory is empty because we 
%                                      do not have permission to redistribute
%                                      this toolbox.
%                                      
%     Subdirectory METHODS:   it contains scripts to analyze on the fly with 
%                             various coevolution detection methods an MSA
%                             simulation carried out with MSAvolve.
%                                    
%--------------------------------------------------------------------------
% Directory TUTORIAL: it contains a simple tutorial on how to use MSAvolve to
%                     simulate MSAs and to analyze the simulations with 
%                     different coevolution detection methods.
%--------------------------------------------------------------------------

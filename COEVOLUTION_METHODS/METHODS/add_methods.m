% This scripts adds coevolution analysis by multiple methods to a MSAvolve
% run. Comment out what you don't want. Remember that MI and MIP are
% required before running ZPX, ZPX2, ZRES, ZRES2. A license from Rama
% Ranganthan is required to run the ramaSCA method. Original code to run 
% the logR method can be obtained from Lukas Burger (lukas.burger@fmi.ch), 
% and must be placed in the directory from which you are running.

run add_MI_method
run add_MIP_method
% run add_ZPX_method
run add_ZPX2_method
% run add_ZRES_method
% run add_ZRES2_method
run add_logR_method
run add_DCA_method
run add_dbZPX2_method
run add_dgbZPX2_method
run add_nbZPX2_method
run add_OMES_method
% run add_PSICOV_method
run add_localPSICOV_method
run add_fpcZPX2_method
% run add_sharmaSCA_method
% run add_ramaSCA_method

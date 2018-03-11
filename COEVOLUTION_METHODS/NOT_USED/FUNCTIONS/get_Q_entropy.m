% In this script we process a set of simulated MSAs producing for each msa
% various types of differential binary msa. Then various forms of entropy
% are calculated for each original and binary msa. The most interesting forms
% are contained in columns 1-4 of entropy_vec. The same four are
% repeated as zscores in columns 14-17:
% entropy_vec(:,1) = entropy of the regular MSA
% entropy_vec(:,2) = entropy of the differential binary MSA based on 0's
% entropy_vec(:,3) = entropy of the differential binary MSA based on 1's
% entropy_vec(:,4) = Shannon binary entropy of the differential binary MSA

end_cycle = size(MSA_select_ALL,3);
entropy_vec = zeros(end_cycle,14);
for i = 1:end_cycle
    nmsa = MSA_select_ALL(:,:,i);
    entropy_vec(i,1) = Entropy(nmsa(:));
    [dif_binmsa,cum_dif_binmsa,dif_binmsa_2,cum_dif_binmsa_2,...
        dif_binmsa_3,cum_dif_binmsa_3] = ...
        nmsa_to_dif_binmsa(nmsa);
    [ entropy_vec(i,2),entropy_vec(i,3),entropy_vec(i,4) ] = ...
        binmsa_to_Q_entropy( dif_binmsa );
    [ entropy_vec(i,5),entropy_vec(i,6),~ ] = ...
        binmsa_to_Q_entropy( cum_dif_binmsa );
    [ entropy_vec(i,7),entropy_vec(i,8),~ ] = ...
        binmsa_to_Q_entropy( cum_dif_binmsa_2 );
    [ entropy_vec(i,9),entropy_vec(i,10),~ ] = ...
        binmsa_to_Q_entropy( cum_dif_binmsa_3 );
    entropy_vec(i,11) = Entropy(cum_dif_binmsa(:));
    entropy_vec(i,12) = Entropy(cum_dif_binmsa_2(:));
    entropy_vec(i,13) = Entropy(cum_dif_binmsa_3(:));
    
end
entropy_vec(:,14:26) = zscore(entropy_vec(:,1:13));

% plot(entropy_vec(:,14),'k'); hold on
% plot(entropy_vec(:,18),'r');
% plot(entropy_vec(:,19),'g');
% plot(entropy_vec(:,20),'b'); 
% plot(entropy_vec(:,21),'y');
% plot(entropy_vec(:,22),'c');
% plot(entropy_vec(:,23),'m');
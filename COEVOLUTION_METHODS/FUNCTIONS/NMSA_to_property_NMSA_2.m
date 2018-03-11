function [pnmsa,A_hat] = NMSA_to_property_NMSA_2(nmsa,dist_method,threshold,...
    psc_method,psc_lambda,psc_blosum,nsymbols,recover_method,...
    recover_level,property,replace_gaps)
% This function converts every numeric symbol in the msa to a property of
% that symbol in each column. Possible properties are:
%
% 'FREQ' = frequency of each residue type in each column.
% 'RELH' = relative entropy of each residue type in each column without 
%           any weights or pseudocounts (recommended).
% 'MV' = mean residue volume.
% 'VDWV' = van der Waals volume.
% 'MZ' = Moret-Zebende scale of polarity based on relative SASA.
% 'ASA' = side chain accessible surface area.
% 'rSASA' = relative side chain surface area.
% 'VR' = average volume of buried residues (A^3).
% 'CHG' = charge
% 'AROM' = aromatic side chain interaction.
% 'BRANCH' = branched side chain interaction.
% 'OMIS' = Stephen White's Octanol-Interface scale of free energy of
%      transfer from water to interface or to octanol and their difference.
%
% 'Weights and pseudocounts are applied to the
% frequencies to correct for the similarities between sequences.
% dist_method = 'GAPS'|'NOGAPS' (include or not gaps in the distance matrix).
% threshold = (similarity threshold: e.g. 0.9).
% psc_method = 'DCA'|'SDP' (pseudocount method DCA or SDP style as described 
% in "SDPpred ..." Kalinina, OV et al. Nucleic Acid Research 2004, Vol32, 
% W424-W428, and in "H2r ..." Merkl, R. and Zwick, M., BMC Bioinformatics 
% 2008, 9,151).
% psc_lambda = pseudocount scale (e.g. 0.5 for DCA; 1.0 for SDP).
% psc_blosum = 'AUTO'|30|35|40|45|50|55|60|62|65|70|75|80|85|90|100. Blosum  
% matrix used to determine the pseudocounts in the SDP method. If set to 
% 'AUTO' the best blosum matrix for the mean value of the distance matrix
% is chosen automatically.
% nsymbols = 20|21 (if set to 20 gaps frequencies are zeroed; 20 is 
% recommended if using psc_method = 'DCA'.
% If 'replace_gaps' is set to 'REPLACE' (as opposed to 'KEEP') gaps in each
% column are replaced with the mean value of the chosen property in that
% column. However, this statement does not affect 'FREQ' and 'RELH'.
% recover_method = 'ALM'|'iALM'|'APG'|'pAPG'|'NONE'; it is the robust PCA
% method used to recover the low rank matrix corresponding to the frequency
% converted msa with noise removed. It uses one of the methods developed by
% the Perception and Decision Lab at the University of Illinois
% (http://perception.csl.uiuc.edu/). The options are: 'exact ALM' (ALM)[can
% be very slow], 'inexact ALM' (iALM)[strongly recommended], 'accelerated
% proximal gradient' (APG) , 'partial accelerated proximal gradient' (pAPG)
% [recommended], or no recovering (NONE). 
% recover_level = specifies the degree of noise removed: smaller values remove 
% more noise, larger values leave the msa more similar to the original one.
% Useful values are in the range 0.2-0.025. If using recovering (e.g.
% recover_method = iALM), the recommended value is recover_level = 0.05.

% Possible usage:
% [pnmsa] = NMSA_to_property_NMSA(nmsa,'NOGAPS',0.90,'SDP',0.5,'AUTO',21,...
%                                  'NONE',0.05,'CHG');
% [pnmsa,noise_red_pnmsa] = NMSA_to_property_NMSA(nmsa,'GAPS',0.90,...
%                                  'SDP',1.0,'AUTO',21,'iALM',0.05,'FREQ');

%-------------------------------------------------------------------------- 
% Declare the property vectors.
% AA  = [ 'A' 'R' 'N' 'D' 'C' 'Q' 'E' 'G' 'H' 'I' 'L' 'K' 'M' 'F' 'P' ...
% 'S' 'T' 'W' 'Y' 'V']';
% Amino acid polarity.
% AP = [9 15 16 19 7 17 18 11 10 1 3 20 5 2 13 14 12 6 8 4]';
% BF (ASA rescaled to MZ)
% BF = [0.156 0.112 0.107 0.087 0.230 0.103 0.103 0.164 0.173 0.217 0.204 ...
%     0.059 0.204 0.217 0.112 0.121 0.138 0.203 0.164 0.208]';
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

[nrows,ncols] = size(nmsa);
ncols_o = ncols;

% First we replace unusual symbols possibly present in the nmsa

ind_25 = nmsa == 25;
nmsa(ind_25) = 21;
ind_23 = nmsa == 23;
nmsa(ind_23) = 21;
ind_22 = nmsa == 22;
nmsa(ind_22) = 21;
ind_0 = nmsa == 0;
nmsa(ind_0) = 21;

switch dist_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(nmsa);
        dist = bin_ordered*bin_ordered';

        sdist = zeros(nrows,nrows);      
        for i = 1:nrows
            sdist(i,:) = dist(i,:)/dist(i,i);
        end

        udist = triu(sdist,1);
        ldist = tril(sdist,-1)';
        lind = udist < ldist;
        udist(lind) = 0;
        ldist(~lind) = 0;
        mdist = udist + ldist;
        dist = mdist + mdist' + eye(nrows);
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;
end

% Here we calculate the weights for the similarity between sequences.

if threshold == 1
    W = ones(nrows,1);
    Meff = nrows;
else
    dist_threshold = dist >= threshold;
    W = 1./sum(dist_threshold)';
    Meff=round(sum(W));
end

fprintf('Meff = %d \n', Meff);

% Here we create a 3d matrix in which every layer has the dimensions of the
% nmsa and represents one of 20 symbols (if gaps are excluded). Then every
% row of each layer is scaled by the weight of that sequence.

nmsa_3 = zeros(nrows,ncols,nsymbols);
for i = 1:nsymbols    
    nmsa_3(:,:,i) = nmsa == i;        
end

% Here we add a total extra fraction of 1 aa equal to the pseudocount
% value to the columns of each slice that count 0; the fraction is added
% only to the rows which are not gaps. Then every row of each slice is
% scaled by the weight of that sequence, based on the similarity threshold.
% Finally we determine the mean of each column counting only the rows
% in which there are no gaps, and we subtract this mean from each column.

profile = zeros(nsymbols,ncols);

for i = 1:ncols
    slice = squeeze(nmsa_3(:,i,:));
    slice_zeros = slice == 0;
    row_gap_ind = nmsa(:,i) == 21;
    row_nongap_ind = nmsa(:,i) ~= 21;
    slice_zeros(row_gap_ind,:) = false;
    col_sum = sum(slice,1);
    col_nongap_ind = col_sum ~= 0;
    slice_zeros(:,col_nongap_ind) = false;
    loq = psc_lambda/sum(slice_zeros(:));
    slice(slice_zeros) = loq;
    slice = slice.*W(:,ones(1,nsymbols));
    % The following line replaces dividing by the no. of observations.
    slice = slice/sum(slice(:));
%     slice_mean = mean(slice(row_nongap_ind,:));
%     slice_mean_mat = zeros(nrows,nsymbols);
%     row_non_gap_ind = find(row_nongap_ind);
%     for k = 1:length(row_non_gap_ind)
%         k_ind = row_non_gap_ind(k);
%         slice_mean_mat(k_ind,:) = slice_mean;
%     end
%     slice = slice - slice_mean_mat;
    nmsa_3(:,i,:) = slice;
    profile(:,i) = sum(slice)';
end

% nmsa = zeros(nrows,0);
% for i = 1:ncols
%     nmsa = [nmsa squeeze(nmsa_3(:,i,:))];
% end

[~,ncols] = size(nmsa);

%--------------------------------------------------------------------------

% Create an msa in which we replace the aa integer representation
% with a property representation.  

switch property
%-------------------------------------------------------------------------
    case 'FREQ'
% Here we create a new msa in which every symbol is replaced by the
% weigthed frequency of that symbol in that column of the msa.

pnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        pnmsa(i,j) = profile(row(j),j);
        end
end

%-------------------------------------------------------------------------
    case 'ddG'
% Here we create a new msa in which every symbol is replaced by the
% anticipated DeltaDeltaG contribution of that symbol in that column of the
% msa to the stability of the protein. The substitution is based on the 
% approximation that ddG ~ -ln(frequency)

% Here we calculate the profile including gaps.
%     profile = zeros(nsymbols,ncols);
%     
%     for i = 1:nsymbols
%         layer = nmsa == i;
%         profile(i,:) = mean(layer); 
%     end
% Here we approximate the ddG contributions.
pnmsa = zeros(nrows,ncols);

if nsymbols == 20
    profile = profile - profile/nrows;
    profile = [profile ; ones(1,ncols)/nrows];
end
    
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
            freq_ij = profile(row(j),j);
            pnmsa(i,j) = -freq_ij*log(freq_ij);
        end
end

% Here we scale each row for the similarity between sequences.
% Wmat = W(:,ones(1,ncols));
% pnmsa = pnmsa.*Wmat;
% pnmsa = nantozero(pnmsa);

%--------------------------------------------------------------------------
    case 'RELH'
        
%     profile = zeros(21,ncols);
%     
%     for i = 1:21
%         layer = nmsa == i;
%         profile(i,:) = mean(layer); 
%     end
    
    bg_frequencies = mean(profile,2);
            
    [relh_aa,relh] = rel_entropy_nmsa_gaps(profile,bg_frequencies,ncols);

%
RELH_profile = figure; clf; 
set(RELH_profile,'Units','normalized','Position',[0 0.3 0.5 0.7],'Name',...
    'REF Sequence Conservation');

subplot(2,1,1);
imagesc(relh_aa,'DisplayName','REF_rel_entropy_aa');figure(gcf)
xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
ylabel('AAs Relative Entropy','FontSize',12,'FontWeight','n');
set(gca,'Xlim',[1 ncols],'YTickLabel',{'A','R','N','D','C','Q',...
    'E','G','H','I','L','K',...
    'M','F','P','S','T','W','Y','V','-'},...
    'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21],...
    'YDir','reverse',...
    'TickDir','in',...
    'Layer','top');

subplot(2,1,2);
stairs(relh,'-r','DisplayName','REF_rel_entropy_sum');figure(gcf)
xlabel('Seq. Position','FontSize',14,'FontWeight','n'); 
ylabel('Total Relative Entropy','FontSize',14,'FontWeight','n'); 
set(gca,'Xlim',[1 ncols]);

% Here we create a new msa in which every symbol is replaced by the
% relative entropy value of that symbol in that column of the msa.

pnmsa = zeros(nrows,ncols);
for i = 1:nrows
        row = nmsa(i,:);
        for j = 1:ncols
        pnmsa(i,j) = relh_aa(row(j),j);
        end
end

% % Here we scale each row for the similarity between sequences:
% Wmat = W(:,ones(1,ncols));
% pnmsa = pnmsa.*Wmat;
% pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'MZ'
% Moret-Zebende polarity scale based on relative SASA.
MZ = [0.157 0.078 0.113 0.087 0.246 0.105 0.094 0.156 0.152 0.222 0.197 ...
    0.069 0.221 0.218 0.121 0.100 0.135 0.174 0.222 0.238 NaN]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = MZ(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'klMZ'
% Moret-Zebende polarity scale based on relative SASA.
MZ = [0.157 0.078 0.113 0.087 0.246 0.105 0.094 0.156 0.152 0.222 0.197 ...
    0.069 0.221 0.218 0.121 0.100 0.135 0.174 0.222 0.238 NaN]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = MZ(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
bg_value = nanmean(pnmsa);

for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = bg_value(i);
end

% Kullback-Leibler divergence of individual residues property from the
% mean value in each column: here we treat each property value as a
% probability;
for i = 1:ncols
    pnmsa(:,i) = pnmsa(:,i).*log(pnmsa(:,i)/bg_value(i));
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'klMZ2'
% Moret-Zebende polarity scale based on relative SASA.
MZ = [0.157 0.078 0.113 0.087 0.246 0.105 0.094 0.156 0.152 0.222 0.197 ...
    0.069 0.221 0.218 0.121 0.100 0.135 0.174 0.222 0.238 NaN]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = MZ(i);
end

% Kullback-Leibler divergence of individual residues property from the mean
% value in each column: here we treat each property value as a probability;
% then we replace any gaps in a column with the mean property value in the
% ungapped rows of that column;

bg_value = nanmean(pnmsa(:));
for i = 1:ncols
    nan_ind = isnan(pnmsa(:,i));
    pnmsa(~nan_ind,i) = pnmsa(~nan_ind,i).*log(pnmsa(~nan_ind,i)/bg_value);
    pnmsa(nan_ind,i) = mean(pnmsa(~nan_ind,i));
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'VDWV'
% Van der Waal volume (A^3)
VDWV = [67 148 96 91 86 114 109 48 118 124 124 135 124 135 90 73 93 163 ...
    141 105 0]';
        
pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = VDWV(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'MV'
% Mean volume (A^3) Gerstein 1999. ***
MV = [89.3 190.3 122.4 114.4 102.5 146.9 138.8 63.8 157.5 163.0 163.1 ...
    165.1 165.8 190.8 121.6 94.2 119.6 226.1 194.6 138.2 0]';
        
pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = MV(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'klMV'
% Mean volume (A^3) Gerstein 1999.
MV = [89.3 190.3 122.4 114.4 102.5 146.9 138.8 63.8 157.5 163.0 163.1 ...
    165.1 165.8 190.8 121.6 94.2 119.6 226.1 194.6 138.2 NaN]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = MV(i);
end
bg_value = nanmean(pnmsa);

% Here we replace any gaps in a column with the mean property value of that
% column:
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = bg_value(i);
end

% Kullback-Leibler divergence of individual residues property from the
% mean value in each column: here we treat each property value as a
% probability;
for i = 1:ncols
    pnmsa(:,i) = pnmsa(:,i).*log(pnmsa(:,i)/bg_value(i));
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%-------------------------------------------------------------------------
    case 'VR'
% Average volume of buried residues (A^3).
VR = [92 225 135 125 106 161 155 66 167 169 168 171 171 203 129 99 122 ...
    240 203 142 0]';
        
pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = VR(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%--------------------------------------------------------------------------
    case 'ASA'
% Side chain Accessible Surface Area (ASA A^2).
ASA = [67 196 113 106 104 144 138 0 151 140 137 167 160 175 105 80 102 ...
    217 187 117 0]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = ASA(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%--------------------------------------------------------------------------
    case 'rSASA'
% Side chain Accessible Surface Area (ASA A^2) (based on Durham et al: Solvent 
% accessible surface area approximation for rapid and accurate protein 
% structure prediction .
rSASA = [209.02 335.73 259.85 257.99 240.50 286.76 285.03 185.15 290.04 ...
    273.46 278.44 303.43 291.52 311.30 235.41 223.04 243.55 350.68 ...
    328.82 250.09 0]';
        
pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = rSASA(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%--------------------------------------------------------------------------
    case 'CHG'
% Charge
CHG = [0 1 0 -1 0 0 -1 0 0.1 0 0 1 0 0 0 0 0 0 0 0 0]';
        
pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = CHG(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%-------------------------------------------------------------------------
    case 'klCHG' % NOTICE: it cannot be used with noise reduction because it 
                 % produces complex numbers.
% Charge
CHG = [0 1 0 -1 0 0 -1 0 0.1 0 0 1 0 0 0 0 0 0 0 0 NaN]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = CHG(i);
end
% bg_value = nanmean(pnmsa);

% Here we replace any gaps in a column with the mean property value of that
% column:
bg_value = zeros(1,ncols);
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    if i == 1
    bg_value(i) = mean(nanmean(pnmsa(:,i:i+1)));
    elseif i == ncols 
    bg_value(i) = mean(nanmean(pnmsa(:,i-1:i)));
    else    
    bg_value(i) = mean(nanmean(pnmsa(:,i-1:i+1)));
    end
    pnmsa(nanind,i) = bg_value(i);
end

% Kullback-Leibler divergence of individual residues property from the
% mean value in each column: here we treat each property value as a
% probability;
for i = 1:ncols
    pnmsa(:,i) = pnmsa(:,i).*log(pnmsa(:,i)/bg_value(i));
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

%--------------------------------------------------------------------------
    case 'OMIS'
% Octanol - Interface Scale Stephen White scale of free energy of transfer 
% from water to interface or to octanol and their difference.
OMIS = [0.33 1.00 0.43 2.41 0.22 0.19 1.61 1.14 -0.06 -0.81 -0.69 1.81 ...
    -0.44 -0.58 -0.31 0.33 0.11 -0.24 0.23 -0.53]';

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = OMIS(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%-------------------------------------------------------------------------
    case 'AROM'
% Aromatic interactions.
AROM = [0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 0 0];

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = AROM(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%-------------------------------------------------------------------------
    case 'BRANCH'
% Branched chain interactions.
BRANCH = [0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 0];

pnmsa = nan(nrows,ncols);
for i = 1:nsymbols
    nmsa_ind = nmsa == i;
    pnmsa(nmsa_ind) = BRANCH(i);
end

% Here we replace any gaps in a column with the mean property value of that
% column:
if strcmp(replace_gaps,'REPLACE')
for i = 1:ncols
    nanind = isnan(pnmsa(:,i));
    pnmsa(nanind,i) = nanmean(pnmsa(:,i));
end
end

% Here we scale each row for the similarity between sequences:
Wmat = W(:,ones(1,ncols));
pnmsa = pnmsa.*Wmat;
pnmsa = nantozero(pnmsa);

% Finally we scale each column for the pseudocount.
% for i = 1:ncols
%      for j = 1:nsymbols
%          ind = nmsa(:,i) == j;
%          pnmsa(ind,i) = pnmsa(ind,i) + psc_scale(j,i)*pnmsa(ind,i);
%      end
% end

%-------------------------------------------------------------------------

end

% Noise reduction.
% First we recover the msa matrix with noise removed using a technique of
% robust PCA.
switch recover_method
    case 'ALM'
    [A_hat E_hat iter] = exact_alm_rpca(pnmsa,recover_level);
    case 'iALM'        
    [A_hat E_hat iter] = inexact_alm_rpca(pnmsa,recover_level);
    case 'APG'
    [A_hat E_hat iter] = proximal_gradient_rpca(pnmsa,recover_level);
    case 'pAPG'
    [A_hat E_hat iter] = partial_proximal_gradient_rpca(pnmsa,recover_level);
    case 'NONE'
    A_hat = pnmsa;
end

end

function [dist] = get_distance_matrix(nmsa,dist_method)

[nrows,ncols] = size(nmsa);

switch dist_method
    case 'NOGAPS'
        bin_ordered = nmsa_to_binmsa_20q(nmsa);        
        dist = bin_ordered*bin_ordered';

        sdist = zeros(nrows,nrows);      
        for i = 1:nrows
            sdist(i,:) = dist(i,:)/dist(i,i);
        end

        udist = triu(sdist,1);
        ldist = tril(sdist,-1)';
        lind = udist < ldist;
        udist(lind) = 0;
        ldist(~lind) = 0;
        mdist = udist + ldist;
        mdist = mdist + mdist' + eye(nrows);

        dist = mdist;
        
    case 'GAPS'
        bin_ordered = nmsa_to_binmsa_21q(nmsa);
        dist = (bin_ordered*bin_ordered')/ncols;

end

end

function [rel_entropy_aa,rel_entropy] = ...
    rel_entropy_nmsa_gaps(profile,bg_prob,ncols)
% This function returns a traditional expression for relative entropy 
% including the gaps count.

rel_entropy_aa = zeros(21,ncols);

for i = 1:ncols
    rel_entropy_aa(:,i) = ...
        nantozero(profile(:,i) .* log(profile(:,i)./bg_prob));
end

rel_entropy = sum(rel_entropy_aa);
end
function [ZPX2,vZPX2,tZPX2,zZPX2,MI,vMI,vMIP,vMI_3d_orig,vMI_3d,tMI,tMIP,...
    tMI_3d_orig,tMI_3d,zMI,zMIP,zMI_3d_orig,zMI_3d,Hi,gW,lgW,...
    v_profile,t_profile,z_profile] = NMSA_to_mdZPX2(nmsa,gapW_exp)

% The direct coupling vMI_3d is the negative of the 'interaction
% information' minus the mutual information. tMI_3d is the 'total
% interaction' between the three variables. zMI_3d is a weighted sum of 
% vMI_3d and tMI_3d. 'gapW_exp' is the exponent of the weights on the 
% 3rd dimension: 0 gives a weight of 1 to everything; 1 leaves the weight
% as calculated by the 'correct_coevmat_forgaps' function; higher powers
% make the weight progressively smaller.

    [~,ncols] = size(nmsa);
    
    % MI, vMI, wvMI 
    [Hi,MI,vMI_3d,tMI_3d,zMI_3d] = NMSA_to_vMI_3d(nmsa,ncols);
    
    % Correct for gaps along the 3rd dimension-----------------------------
    vMI_3d_orig = vMI_3d;
    tMI_3d_orig = tMI_3d;
    zMI_3d_orig = zMI_3d;
    gapW = correct_coevmat_forgaps(nmsa);
    gW = gapW.^gapW_exp;
    lgW = mean(gW);
    for k = 1:ncols
    vMI_3d(:,:,k) = vMI_3d(:,:,k)*lgW(k);
    tMI_3d(:,:,k) = tMI_3d(:,:,k)*lgW(k);
    zMI_3d(:,:,k) = zMI_3d(:,:,k)*lgW(k);
    end
    % ---------------------------------------------------------------------
   
    vMI = squeeze(nanmean(vMI_3d,3));
    tMI = squeeze(nanmean(tMI_3d,3));
    zMI = squeeze(nanmean(zMI_3d,3));
   
    % MIP
    MIP = MI_to_MIP(MI,ncols);
    vMIP = MI_to_MIP(vMI,ncols);
    tMIP = MI_to_MIP(tMI,ncols);
    zMIP = MI_to_MIP(zMI,ncols);
    
    % ZPX2
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    vZPX2 = MIP_to_ZPX2(vMIP,ncols);
    tZPX2 = MIP_to_ZPX2(tMIP,ncols);
    zZPX2 = MIP_to_ZPX2(zMIP,ncols);
    
    % Mean interaction profile based on MIs
    v_profile = zeros(ncols,1);
    t_profile = zeros(ncols,1);
    z_profile = zeros(ncols,1);
    for k = 1:ncols
        v_level = (MI-vMI_3d_orig(:,:,k))*lgW(k);
        v_profile(k) = nanmean(v_level(:));
        t_level = (MI-tMI_3d_orig(:,:,k))*lgW(k);
        t_profile(k) = nanmean(t_level(:));
        z_level = (MI-zMI_3d_orig(:,:,k))*lgW(k);
        z_profile(k) = nanmean(z_level(:));
    end
                    
end


function [Hi,MI,vMI_3d,tMI_3d,zMI_3d] = NMSA_to_vMI_3d(nmsa,ncols)
%--------------------------------------------------------------------------
% MI, vMI_3d, tMI_3d calculation.

MI = NaN(ncols,ncols);
vMI_3d = NaN(ncols,ncols,ncols);
tMI_3d = NaN(ncols,ncols,ncols);
zMI_3d = NaN(ncols,ncols,ncols);

Hi = NaN(ncols,1);
for i = 1:ncols
        Hi(i) = Entropy(nmsa(:,i));
end

Hij = NaN(ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols
        
        Hij(i,j) = JointEntropy(nmsa(:,[i j]));
        Hij(j,i) = Hij(i,j);
        
        MI(i,j) = Hi(i) + Hi(j) - Hij(i,j);
        MI(j,i) = MI(i,j);
        
    end
end

Hijk = NaN(ncols,ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols
        for k = j+1:ncols
        
            Hijk(i,j,k) = JointEntropy(nmsa(:,[i j k]));

            % We store also all the permutations. 
        
            Hijk(i,k,j) = Hijk(i,j,k);
            Hijk(j,i,k) = Hijk(i,j,k);
            Hijk(j,k,i) = Hijk(i,j,k);
            Hijk(k,i,j) = Hijk(i,j,k);
            Hijk(k,j,i) = Hijk(i,j,k);
            
        end        
    end
end

all_ind = (1:ncols);
for i = 1:ncols
    for j = i+1:ncols
        k_ind = setdiff(all_ind,[i j]);
        for k = k_ind
        
            % The following is the interaction information minus the MI.
            % vMI_3d(i,j,k) = Hi(k) - Hij(j,k) - Hij(i,k) + Hijk(i,j,k);
            % We take directly the negative of it, which is equal to 
            % Sunil Srinivasa expression.
            % vMI_3d(i,j,k) = Hij(i,k) - Hijk(i,k,j) - Hi(k) + Hij(k,j);
            % Since we are going to average over the 3rd dimension, here
            % only i and j can be permuted.
            
            vMI_3d(i,j,k) = Hij(j,k) + Hij(i,k) - Hi(k) - Hijk(i,j,k);
            vMI_3d(j,i,k) = vMI_3d(i,j,k);
            tMI_3d(i,j,k) = Hi(i) + Hi(j) + Hi(k) - Hijk(i,j,k);
            tMI_3d(j,i,k) = tMI_3d(i,j,k);
            zMI_3d(i,j,k) = vMI_3d(i,j,k) + 0.5*tMI_3d(i,j,k);
            zMI_3d(j,i,k) = zMI_3d(i,j,k);
                    
        end        
    end
end

end


function [pMIP] = MI_to_MIP(pMI,ncols)
%--------------------------------------------------------------------------
% MIP calculation
mean_mat = nanmean(pMI(:));
mean_row = zeros(ncols,1);
MCA_mat = zeros(ncols,ncols);

% Here  we calculate the MCA matrix.

for m = 1:ncols
    mean_row(m) = nanmean(pMI(m,:));
end

for m = 1:ncols
    for n = m:ncols
    MCA_mat(m,n)=(mean_row(m)*mean_row(n))/mean_mat;
    MCA_mat(n,m) = MCA_mat(m,n);    
    end
MCA_mat(m,m) = NaN;
end

% Finally we subtract the MCA matrix from the MI matrix.
pMIP = pMI-MCA_mat;

end


function [pZPX2] = MIP_to_ZPX2(pMIP,ncols)
%--------------------------------------------------------------------------
% ZPX2 calculation
mean_row=zeros(ncols,1);
std_row=zeros(ncols,1);
pZPX2=zeros(ncols,ncols);

for m=1:ncols
    mean_row(m)=nanmean(pMIP(m,:));
    std_row(m)=nanstd(pMIP(m,:));   
end
for m=1:ncols
    for n=m:ncols

    ZPX2_i=(pMIP(m,n)-mean_row(m))/std_row(m);
    ZPX2_j=(pMIP(m,n)-mean_row(n))/std_row(n);
        
    pZPX2(m,n)=ZPX2_i*ZPX2_j;

% Here we correct for the product of two negative ZPX2_i and ZPX2_j, which
% would give the wrong MI. Comment: I am not sure the following three lines
% make a difference. Change of sign is not in the original ZPX2 algorithm by
% Gloor, but is included in the ZRES algorithm by Chen.
    
    if (ZPX2_i<0&&ZPX2_j<0)
        pZPX2(m,n)=-pZPX2(m,n);
    end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end




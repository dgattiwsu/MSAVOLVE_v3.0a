function [ZPX2,vZPX2,MI,vMI,Hi,vMI_3d,profile] = ...
    NMSA_to_vZPX2_2(nmsa)
% The direct coupling vMI_3d is the negative of the 'interaction
% information' minus the mutual information. An additional weighing
% (wvMI_3d) is possible based on the entropy of the 3rd position, on the
% assumption that more conserved 3rd positions are more likely to affect
% the interactions between the 1st two positions.

    [~,ncols] = size(nmsa);


    % MI, vMI, wvMI 
    [Hi,MI,vMI_3d] = NMSA_to_vMI_3d(nmsa,ncols);
    
    gapW0 = ones(ncols,ncols);
    gapW1 = correct_coevmat_forgaps(nmsa);
    gapW2 = gapW1.^2;
    gapW3 = gapW1.^3;    
    gW = gapW3;
    MI = MI.*gW;
    for k = 1:ncols
    vMI_3d(:,:,k) = vMI_3d(:,:,k).*gW;
    end
    
    vMI = MI - squeeze(nanmean(vMI_3d,3));
   
    % MIP
    MIP = MI_to_MIP(MI,ncols);
    vMIP = MI_to_MIP(vMI,ncols);
    
    % ZPX2
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    vZPX2 = MIP_to_ZPX2(vMIP,ncols);
    
    % Mean interaction profile based on MIs
    profile = zeros(ncols,1);
    for k = 1:ncols
        level = MI-vMI_3d(:,:,k);
        profile(k) = nanmean(level(:));
    end
                    
end


function [Hi,MI,vMI_3d] = NMSA_to_vMI_3d(nmsa,ncols)
%--------------------------------------------------------------------------
% MI and vMI_3d calculation.

MI = NaN(ncols,ncols);
vMI_3d = NaN(ncols,ncols,ncols);

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
        
            % The following is the interaction information.
            vMI_3d(i,j,k) = ...
                MI(i,j) + Hi(k) - Hij(j,k) - Hij(i,k) + Hijk(i,j,k);
            % We take directly the negative of it, which is equal to 
            % Sunil Srinivasa expression.
            % vMI_3d(i,j,k) = Hij(i,k) - Hijk(i,k,j) - Hi(k) + Hij(k,j);
            
            % vMI_3d(i,j,k) = Hij(j,k) + Hij(i,k) - Hi(k) - Hijk(i,j,k);
            vMI_3d(j,i,k) = vMI_3d(i,j,k); 
                    
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
% Gloor, but is included in the ZRES algorithm by Chen. At any rate the
% change of sign seems to affect only i,j positions with very small counts,
% and therefore the final effect of the change is marginal.
    
%     if (ZPX2_i<0&&ZPX2_j<0)
%         pZPX2(m,n)=-pZPX2(m,n);
%     end

% Symmetrize.

    pZPX2(n,m)=pZPX2(m,n);
    
    end
    
    pZPX2(m,m)=NaN;

end
end




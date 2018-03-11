function [v1_ZPX2,v2_ZPX2,MI,MI_3d,mean_MI_3d,ZPX2,mean_ZPX2_3d] = ...
    NMSA_to_v3_ZPX2(nmsa)
% The MI_3d and ZPX2_3d represent the 'interaction information' between
% three positions. In order to get the direct coupling between two
% positions we subtract the interaction information from the mutual
% information.

% Here we calculate the ZPX2 matrix in which the effect of a 3rd variable 
% on the 1st 2 variables has been removed.

    [~,ncols] = size(nmsa);
    
    % MI and MI_3d matrix. 
    [MI,MI_3d] = NMSA_to_MI_3d(nmsa,ncols);
    mean_MI_3d = squeeze(nanmean(MI_3d,3));
    vMI = MI - mean_MI_3d;
    
    % MIP3d matrix
    vMIP = MI_to_MIP(vMI,ncols);
    
    % and finally vZPX2 matrix
    v1_ZPX2 = MIP_to_ZPX2(vMIP,ncols);
    
    % Next we consider the standard ZPX2
    MIP = MI_to_MIP(MI,ncols);
    ZPX2 = MIP_to_ZPX2(MIP,ncols);
    
    % Then the mean ZPX2_3d 
    MIP_3d = MI_to_MIP(mean_MI_3d,ncols);    
    mean_ZPX2_3d = MIP_to_ZPX2(MIP_3d,ncols);
    v2_ZPX2 = ZPX2 - mean_ZPX2_3d;
            
end


function [MI,MI_3d] = NMSA_to_MI_3d(nmsa,ncols)
%--------------------------------------------------------------------------
% MI and MI_3d calculation.

MI = NaN(ncols,ncols);
MI_3d = NaN(ncols,ncols,ncols);
% MI_3d_unique = NaN(ncols,ncols,ncols);

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

all_ind = (1:ncols);
Hijk = NaN(ncols,ncols,ncols);
for i = 1:ncols
    for j = i+1:ncols
        k_ind = setdiff(all_ind,[i j]);
        for k = k_ind
        
            Hijk(i,j,k) = JointEntropy(nmsa(:,[i j k]));
            Hijk(j,i,k) = Hijk(i,j,k);

            % We store also all the permutations. 
        
%             Hijk(i,k,j) = Hijk(i,j,k);
%             Hijk(j,i,k) = Hijk(i,j,k);
%             Hijk(j,k,i) = Hijk(i,j,k);
%             Hijk(k,i,j) = Hijk(i,j,k);
%             Hijk(k,j,i) = Hijk(i,j,k);
            
            
%             MI_3d(i,j,k) = MI(i,j) + Hi(k) ...
%                 - Hij(j,k) - Hij(i,k) + Hijk(i,j,k);

            MI_3d(i,j,k) = Hi(k) ...
                - Hij(j,k) - Hij(i,k) + Hijk(i,j,k);
            
            MI_3d(j,i,k) = MI_3d(i,j,k); 
            
        % This is the unique part of the matrix.
        % MI_3d_unique(i,j,k) = MI_3d(i,j,k);
        
        
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




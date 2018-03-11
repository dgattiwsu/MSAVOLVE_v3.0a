function [ md_MI ] = multidim_MI( MI_mat )

% This function calculates the multidimensional MI starting from an MI
% matrix.

TRACE = trace(MI_mat);

% Here we correct for situations in which there is only one MI value or the
% trace of the MI matrix is 0.

if size(MI_mat,1) == 1 || TRACE == 0

    md_MI = 1;

else
    
    
EIG = eig(MI_mat);

% Here we correct for possibly negative eigenvalues that would give origin
% to complex values of md_MI

EIG(EIG<0)=0;

EIG_TRACE_RATIO = EIG/TRACE;

    if sum(EIG_TRACE_RATIO) == 0
           md_MI = 1;
    else
        
% Remove values of the ratio that are equal to zero to prevent NaN's.

EIG_TRACE_RATIO = nonzeros(EIG_TRACE_RATIO); 

S = sum(EIG_TRACE_RATIO .* log2(EIG_TRACE_RATIO));

S1 = S/log2(length(EIG));

md_MI = 1 + S1;

    end

end 

end


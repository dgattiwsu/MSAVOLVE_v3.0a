function [ t_mut_prob_mat ] = get_blosum_prob( ident,bg_frequencies )

% Here we decide which conservation matrix to use for representing the
% probability of mutation. We recall that the entries Sij in the blosum 
% matrix are Sij=lambda*log2(Mij/pj), where Mij is the probability of 
% i mutating to j and pj is the frequency of j. Lambda is a scaling factor 
% which is equal to 3 for blosum50 and 2 for blosum62 if we use log2. 
% Please notice that if you use the equivalent expression with ln instead 
% of log2 the equation becomes Sij=(1/lambda)*ln(Mij/pj), and lambda
% becomes ln(2)/3 for blosum50.
% Thus, we multiply the ratios by the vector of background
% probabilities to obtain the mutation frequency.
% In calling the substitution matrix, we make sure you are
% taking only the values of the matrix corresponding to the standard aa's by
% adding the flag "'Extended',false".
% All substitution matrices can be downloaded from 
% ftp://ftp.ncbi.nih.gov/blast/matrices/

[blosum_S0,info] = blosum(ident,'Extended',false);
scale = info.Scale;

% Be very careful with the scaling: it changes from matrix to
% matrix! If we are using the matrices downloaded from ncbi the scaling is 
% in units of bits/(scale number). For example, the blosum50 matrix is given  
% in 1/3 bit units. Please note that in principle we could
% use a different substitution matrix at different stages depending on the
% mutation rate we are imposing. For every matrix you can get the actual
% matrix values and the scale (for log2 calculations) with the command:
% [B,MATRIXINFO] = blosum(ident). Remember that Lambda is 1/scale, therefore:

blosum_S1 = blosum_S0*scale;
blosum_S2 = pow2(blosum_S1);

% If the matrix was in nats rather than bits we would use instead:
% blosum_S1 = blosum_S0*log(2)*scale;
% blosum_S2 = exp(blosum_S1);

blosum_S3 = zeros(20,20);
for i = 1:20
   blosum_S3(i,:) = blosum_S2(i,:) .* bg_frequencies';
end

% blosum_S4 = round(blosum_S3*10000);

t_mut_prob_mat = blosum_S3';

%clear blosum* i

end


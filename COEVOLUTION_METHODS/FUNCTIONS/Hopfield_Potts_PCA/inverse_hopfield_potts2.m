function F_apc = inverse_hopfield_potts2( Vtilde, Lambda, N, q, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hopfield-Potts approach to Direct-Coupling Aanalysis 
%
% Copyright for this implementation: 
%             2013 - Martin Weigt
%                    martin.weigt@upmc.fr
% 
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Any publication resulting from applications of DCA should cite:
%
%     S Cocco, R Monasson and M Weigt (2013)
%     From principal component analysis to direct-coupling analysis
%     of coevolution in proteins: Low-eigenvalue modes are needed
%     for structure prediction, PLoS Computational Biology XX, XXXX
%
% Usage
% 
% >> [ Lambda, Vtilde, N, q ] = inverse_hopfield_potts1('imput.fasta',0.2,0.5);
%
% where input.fasta is the input multiple-sequence alignment.
% 0.2 is a typical value for the sequence distance used for reweighting
% (for more than 20,000 sequences in the MSA, no reweighting is used).
% 0.5 is a typical value for the relative weight of the pseudocount.
%
% >> F_apc = inverse_hopfield_potts2( Vtilde, Lambda, N, q, p);
%
% where p has to be replaced by the wanted number of patterns.
% The output matrix F_apc is the score matrix for all position pairs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lambda1 = eye( N*(q-1));

% log-likelihood contribution of patterns 

ll = diag(Lambda) - ones( N*(q-1),1 ) - log( diag( Lambda ) );

[~,b] = sort( ll, 'descend');

% consider only p patterns of highest ll

for i=1:p
    Lambda1( b(i), b(i) ) = Lambda( b(i), b(i) );
end

Lambda2 = ( Lambda1 - eye(N*(q-1)) )/Lambda1;

invC = -Vtilde * Lambda2 * Vtilde';

% sampling-corrected Frobenius norm of the couplings

F_apc = calc_norm_apc( invC, N, q);

end

function A=mapkey(i,alpha,q)
A = (q-1)*(i-1)+alpha;
end

function W=ReturnW(C,i,j,q)

W = zeros(q,q);
W(1:q-1,1:q-1) = C(mapkey(i,1:q-1,q),mapkey(j,1:q-1,q)) ;

end

function F_apc = calc_norm_apc( invC, N, q)

F = zeros(N);
F_apc = zeros(N);

for i = 1:(N-1)
    for j = (i+1):N
        J_mf = ReturnW(invC,i,j,q);
        J_j = mean(J_mf,1);
        J_i = mean(J_mf,2);
        J = mean(J_i);
     
        for a = 1:q
            for b = 1:q
                J_mf(a,b) = J_mf(a,b) - J_i(a) - J_j(b) + J;
            end
        end
        
        F(i,j) = norm( J_mf, 'fro' );
        F(j,i) = F(i,j);
        
    end
end

F_i = mean(F);
F_av = mean(F_i);

for i = 1:(N-1)
    for j = (i+1):N
        F_apc(i,j) = F(i,j) - F_i(i)*F_i(j)/F_av;
        F_apc(j,i) = F_apc(i,j);
    end
end

end

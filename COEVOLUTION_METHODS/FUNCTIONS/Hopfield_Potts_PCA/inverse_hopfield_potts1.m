function [ Lambda, Vtilde, N, q ] = inverse_hopfield_potts1(nomefile, theta, pseudocount_weight)

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

[Pij_true,Pi_true,N,q] = read_alignment(nomefile, theta);
[Pij,Pi] = with_pc(Pij_true,Pi_true,pseudocount_weight,N,q);
clear Pij_true;
clear Pi_true;

[C,C_self] = Compute_C(Pij,Pi,N,q);  % covariance matrix
clear Pij;
clear Pi;

D = sqrtm(C_self);

Gamma = (D\C)/D;     % Pearson correlation matrix

[V,Lambda] = eig(Gamma);  

Vtilde = D\V;    % Patterns up to prefactor sqrt(abs(1-1/lambda))

end

function [Pij_true,Pi_true,N,q] = read_alignment(nomefile, theta)

[N,M,align] = read_alignment_fasta(nomefile); 

% reweighting (only if M<20000 since pdist becomes large)

W = ones(1,M);

if( theta > 0.0 & M<20000)
   W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
end

Meff=sum(W);

q = max(max(align));

[N,M,Meff,q]

% frequencies of amino acid occurrence

Pij_true = zeros(N,N,q,q);
Pi_true = zeros(N,q);

for j=1:M
    for i=1:N
        Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
    end
end
Pi_true = Pi_true/Meff;

for l=1:M
    for i=1:N-1
        for j=i+1:N
            Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
            Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
        end
    end
end
Pij_true = Pij_true/Meff;

scra = eye(q,q);
for i=1:N
    for alpha=1:q
        for beta=1:q
            Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
        end
    end
end
end

function [N,M,Z] = read_alignment_fasta(nomefile)
X = fastaread(nomefile);
N = size(X(1).Sequence,2);
M = size(X,1);
Z = zeros(M,N);
for i=1:M
    for j=1:N
        Z(i,j)=letter2number(X(i).Sequence(j));
    end
end
end

function x=letter2number(a)
switch(a)

    case '-'
         x=1;
    case 'A'    
        x=2;    
    case 'C'    
        x=3;
    case 'D'
        x=4;
    case 'E'  
        x=5;
    case 'F'
        x=6;
    case 'G'  
        x=7;
    case 'H'
        x=8;
    case 'I'  
        x=9;
    case 'K'
        x=10;
    case 'L'  
        x=11;
    case 'M'
        x=12;
    case 'N'  
        x=13;
    case 'P'
        x=14;
    case 'Q'
        x=15;
    case 'R'
        x=16;
    case 'S'  
        x=17;
    case 'T'
        x=18;
    case 'V'
        x=19;
    case 'W'
        x=20;
    case 'Y'
        x=21;
    otherwise
        x=1;
end
end

function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)

Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

scra = eye(q);

for i=1:N
    for alpha = 1:q
        for beta = 1:q
           Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
        end
    end
end 

end

function [C,C_self] = Compute_C(Pij,Pi,N,q)
C=zeros(N*(q-1),N*(q-1));
C_self=zeros(N*(q-1),N*(q-1));
for i=1:N
    for j=1:N
        for alpha=1:q-1
            for beta=1:q-1
                 C(mapkey(i,alpha,q),mapkey(j,beta,q)) = Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
                 if ( i == j )
                     C_self(mapkey(i,alpha,q),mapkey(j,beta,q)) = C(mapkey(i,alpha,q),mapkey(j,beta,q));
                 end
            end
        end
    end
end
end

function A=mapkey(i,alpha,q)
A = (q-1)*(i-1)+alpha;
end

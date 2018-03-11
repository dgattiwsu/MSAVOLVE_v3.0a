% NormMutualInformation: returns mutual information (in bits) of the 'X' and 'Y'
% by Will Dwinnell modified by DG.
%
%
% MI  = calculated mutual information (in bits)
% NMI  = calculated normalized mutual information (in bits)
% X  = variable(s) to be analyzed (column vector)
% Y  = variable to be analyzed (column vector)
%
% Note: Multiple variables may be handled jointly as columns in matrix 'X'.
% Note: Requires the 'Entropy' and 'JointEntropy' functions.
%
% Last modified: Nov-12-2006

function [MI,NMI] = NormMutualInformation(X,Y)

if (size(X,2) > 1)  % More than one predictor?
    % Axiom of information theory
    JE = JointEntropy([X Y]);
%    MI = JointEntropy(X) + Entropy(Y) - JointEntropy([X Y]);
%    NMI = I/JointEntropy([X Y]);
    MI = JointEntropy(X) + Entropy(Y) - JE;
    NMI = MI/JE;
    
else
    % Axiom of information theory
    JE = JointEntropy([X Y]); 
%    I = Entropy(X) + Entropy(Y) - JointEntropy([X Y]);
%    J = I/JointEntropy([X Y]);
    MI = Entropy(X) + Entropy(Y) - JE;
    NMI = MI/JE;    
    
end


% God bless Claude Shannon.

% EOF



% MutualInformation: returns mutual information (in bits) of the 'X' and 'Y'
% by Will Dwinnell, modified by DG.
%
% Note: Multiple variables may be handled jointly as columns in matrix 'X'.
% Note: Requires the 'Entropy' and 'JointEntropy' functions.
%
% Last modified: Nov-12-2006

function SU = SymmetricUncertainty(X,Y)

if (size(X,2) > 1)  % More than one predictor?
    % Axiom of information theory
    JH_X = JointEntropy(X);
    H_Y = Entropy(Y);
    JH_XY = JointEntropy([X Y]);
    SU = 2*(JH_X + H_Y - JH_XY)/(JH_X + H_Y);
else
    % Axiom of information theory
    H_X = Entropy(X);
    H_Y = Entropy(Y);
    JH_XY = JointEntropy([X Y]);
    SU = 2*(H_X + H_Y - JH_XY)/(H_X + H_Y);
end


% EOF



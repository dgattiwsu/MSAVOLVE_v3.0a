% Joint Probability Distribution
% 

function alphabet = JointProbDistr(X)

% Define the alphabet including symbols from 1 to 20 

alphabet = zeros(400,3);

for n = 0:19
k = n+1;
for i = 1:20
    j = (20*n)+i;
        alphabet(j,1) = k;
        alphabet(j,2) = i;
end
end

% Loop through the alphabet.

for i = 1:400;
    
A = X == alphabet(i,1);
B = X == alphabet(i,2);
C = [A(:,1) B(:,2)];
D = C(:,1) == 1 & C(:,2) == 1;
    
alphabet(i,3) = sum(D)/length(X);
end





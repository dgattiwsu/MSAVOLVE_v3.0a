% Probability Distribution
%

function alphabet = ProbDistr(X)

% Define the alphabet including symbols from 1 to 20 

alphabet = zeros(20,2);


for i = 1:20

        alphabet(i,1) = i;
end


% Loop through the alphabet
for i = 1:20
    
A = X == alphabet(i,1);
B = A(:,1) == 1 ;
    
alphabet(i,2) = sum(B)/length(X);
end





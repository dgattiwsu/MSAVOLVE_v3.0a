function H = Entropy20(X)

% Entropy: Returns entropy (in bits) of each column of 'X'

% Establish size of data
[~, m] = size(X);

% Housekeeping
H = zeros(1,m);
const = log2(20);

for Column = 1:m,
    % Assemble observed alphabet
    Alphabet = unique(X(:,Column));
	
    % Housekeeping
    Frequency = zeros(size(Alphabet));
	
    % Calculate sample frequencies
    for symbol = 1:length(Alphabet)
        Frequency(symbol) = sum(X(:,Column) == Alphabet(symbol));
    end
	
    % Calculate sample class probabilities
    P = Frequency / sum(Frequency);
	
    % Calculate entropy in bits
    % Note: floating point underflow is never an issue since we are
    % dealing only with the observed alphabet
    H(Column) = -sum(P .* (log2(P)/const));
end

end
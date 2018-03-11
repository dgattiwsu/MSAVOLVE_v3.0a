%the ln(factorial) function, two ways: by gamma function if <170, and
%by stirlings approximation if >170.
function [m]=rr_fact(n)
%This checks for size and calculates the ln(factorial)
if n<=170
m=gammaln(n+1);
else
m=n*log(n)-n;
end


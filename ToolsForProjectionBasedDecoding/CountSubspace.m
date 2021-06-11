function num = CountSubspace(n,k,s)
%% Count the number of k-dimensional subspaces of EG(n, 2^s)
% Input:
%   n: n-dimensional EG over GF(2^s)
%   k: k-dimensional subspace
%   s: over GF(2^s)
% Output:
%   num: subspace number
% Eg:
% Count the 3-dimensional subspaces of EG(4, 2^2):
% num = CountSubspace(4,3,2);

%%

num = 1;
for i = 1:k
    num = num * (2^(s*(n-i+1))-1) / (2^(s*(k-i+1))-1);
end



end
function [message, vEsti] = hardDecisionDecoding(G, y)
% see the ref:
% Soft-Decision Decoding of Linear Block Codes Based on Ordered Statistics
% Trans. Information Theory 1995.
%%
[k, n] = size(G);
yAbs = abs(y);
% first permutation
[~, indexPermutation1] = sort(yAbs, 'descend');
yReOrder1 = y(indexPermutation1);
GReOrder1 = G(:,indexPermutation1);
% second permutation
[GReOrder1RowEchelon, indexColPivot, ~] = ...
    getEchelonMatrix(GReOrder1);
GReOrder1Eliminated =  backSubstitution(GReOrder1RowEchelon, ...
    indexColPivot, k);
indexColFree = setdiff((1:n),indexColPivot);
indexPermutation2 = nan(1, n);
indexPermutation2(1:k) = indexColPivot;
indexPermutation2(k+1:n) = indexColFree;
GReOrder2 = GReOrder1Eliminated(:,indexPermutation2);
yReOrder2 = yReOrder1(indexPermutation2);
% hard decision decoding
message = (yReOrder2(1:k) < 0);
vEsti = mod(message * GReOrder2, 2);
% inv permutation
vEsti(indexPermutation2) = vEsti;
vEsti(indexPermutation1) = vEsti;

end
function [uEsti, vEsti] = osdDecoding(G, orderL, shiftIndexCell, ...
    listSizeSegment, listSize, y)
% see the ref:
% Soft-Decision Decoding of Linear Block Codes Based on Ordered Statistics
% Trans. Information Theory 1995.
%% Perform hard decision decoding to find the permutaed generator matrix
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
vHardDecision = mod(message * GReOrder2, 2);
% inv permutation
% vEsti(indexPermutation2) = vEsti;
% vEsti(indexPermutation1) = vEsti;
%% perform list decoding
messageList = nan(listSize, k);
metricList = nan(listSize, 1);
codewordList = nan(listSize, n);
pList = 1;
messageList(pList, :) = message;
codewordList(pList, :) = vHardDecision;
metricList(pList, :) = sum((1-2*codewordList(pList, :)) .* yReOrder2);
for iL = 1:orderL
    shiftIndexArray = shiftIndexCell{1,iL};
    for iList = 1:listSizeSegment(iL)
        shiftIndex = shiftIndexArray(iList,:);
        pList = pList+1;
        messageList(pList,:) = message;
        messageList(pList,shiftIndex) = (~message(shiftIndex));
        codewordList(pList, :) = ...
            mod(vHardDecision + sum(GReOrder2(shiftIndex,:),1),2);
        metricList(pList, :) = sum((1-2*codewordList(pList, :)) .* yReOrder2);
    end
end
[~, indexMLinList] = max(metricList);

uEsti = messageList(indexMLinList,:);
vEsti = codewordList(indexMLinList,:);
vEsti(indexPermutation2) = vEsti;
vEsti(indexPermutation1) = vEsti;
end

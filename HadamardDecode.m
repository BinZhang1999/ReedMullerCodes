function [uHat, vHat] = HadamardDecode(hadamardSetting, y)

hadamardMatrix = hadamardSetting.hadamardMatrix;
m = hadamardSetting.m;
N = hadamardSetting.N;
G = hadamardSetting.G;

[nVec, ~] = size(y);
corr = hadamardMatrix * y';
[~, idx] = max(abs(corr),[],1); 

% select the most likely vector
uHat = zeros(nVec, m+1);
uHat(:, m:-1:1) = myDec2Bin(idx-1, m);
idx = N.*(0:nVec-1)+idx';
corr = corr';
uHat(:, m+1) = (corr(idx)<0);
vHat = mod(uHat*G, 2);  

end
function binMatrix = myDec2Bin(decVector, m)
%% Tranform (num)_dec form to binary extension form (num)_bin
% input:
%       decVector: An unsigned int column vector(NX1) 
%       m: The bit width of the binary extension
% output: 
%       binMatrix: An NXm matrix;
%                the i-th row is the binary extension of i-th element of
%                decVector. The left is the higher order.
% Eg:
% % To get the binary representation of vec=[32, 31, 17]'
% binMatrix = myDec2Bin([32, 31, 17]', 6);

%% 
n = length(decVector);
binMatrix = zeros(n, m);
binBasis = 2.^((m-1):-1:0);
binMatrix(:) = mod(floor(decVector ./ binBasis),2);
end
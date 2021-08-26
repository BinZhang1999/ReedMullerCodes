function decVector = myBin2Dec(binMatrix, bitWidth)
%% Tranform (num)_bin form to binary extension form (num)_dec
% input:
%       binMatrix: An NXm matrix;
%                the i-th row is the binary extension of i-th element of
%                decVector. The left is the higher order.
%       bitWidth: The bit width of the binary extension
% output: 
%       decVector: An unsigned int column vector(NX1) 

% Eg:
% % To get the binary representation of vec=[[0 0 0];[0 1 0];[0 1 1]]
% decVector = myBin2Dec([[0 0 0];[0 1 0];[0 1 1]], 3);

%%
binBasis = 2.^((bitWidth-1):-1:0);
binBasis = binBasis';
decVector = binMatrix * binBasis;
end

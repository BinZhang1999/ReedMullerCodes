function num_bin = MyDec2Bin(num_dec, m)
%% Tranform (num)_dec form to binary extension form (num)_bin
% input:
%       num_dec: An unsigned int column vector(NX1) 
%       m: The bit width of the binary extension
% output: 
%       num_bin: An NXm matrix;
%                the i-th row is the binary extension of i-th element of
%                num_dec. The left is the higher order.
% Eg:
% % To get the binary representation of vec=[32, 31, 17]'
% num_bin = MyDec2Bin([32, 31, 17]', 6);

%% 
N = length(num_dec);
num_bin = zeros(N, m);
bin_list = 2.^((m-1):-1:0);
num_bin(:) = mod(floor(num_dec ./ bin_list),2);
end
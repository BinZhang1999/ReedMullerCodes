function num_dec = MyBin2Dec(num_bin, m)
%% Tranform (num)_bin form to binary extension form (num)_dec
% input:
%       num_bin: An NXm matrix;
%                the i-th row is the binary extension of i-th element of
%                num_dec. The left is the higher order.
%       m: The bit width of the binary extension
% output: 
%       num_dec: An unsigned int column vector(NX1) 

% Eg:
% % To get the binary representation of vec=[[0 0 0];[0 1 0];[0 1 1]]
% num_dec = MyBin2Dec([[0 0 0];[0 1 0];[0 1 1]], 3);

%%
bin_list = 2.^((m-1):-1:0);
bin_list = bin_list';
num_dec = num_bin*bin_list;
end

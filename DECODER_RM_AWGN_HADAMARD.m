classdef DECODER_RM_AWGN_HADAMARD
%% Decode first order RM codes using fast hadamard transform
properties
      m;
      K;
      N;
      G;
      hadamard_matrix;
end
   
methods
        function obj = Init(obj, m)
            obj.m = m;
            obj.K = m + 1;
            obj.N = 2^m;
            obj.hadamard_matrix = hadamard(obj.N);
            
            obj.G = zeros(m+1, obj.N);
            obj.G(m:-1:1,:) = (MyDec2Bin((0:(2^m-1))',m))';
            obj.G(m+1,:) = ones(1,obj.N);
        end
        
        function [u_hat, v_hat] = Decode(obj, rx, sigma)
         %% Perform maximum likely hood decoding for RM(1,m)
         % rx: num_rx X N matrix
         % sigma: sigma^2 is the noise power
         %%
            % perform hadamard transfer
             [num_rx, ~] = size(rx);
             me = rx * obj.hadamard_matrix;
             [~, idx] = max(abs(me),[],2); 
             % select the most likely vector
             u_hat = zeros(num_rx,obj.m+1);
             u_hat(:,obj.m:-1:1) = MyDec2Bin(idx-1, obj.m);
             idx = obj.N.*(0:num_rx-1)+idx';
             me = me';
             u_hat(:,obj.m+1) = (me(idx)<0);
             v_hat = mod(u_hat*obj.G, 2);    
        end
    end
    
end
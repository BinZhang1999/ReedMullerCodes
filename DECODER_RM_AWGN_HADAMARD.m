classdef DECODER_RM_AWGN_HADAMARD
%% Decode first order RM codes using fast hadamard transform
properties
      hadamardSetting;
end
   
methods
        function obj = Init(obj, m)
            obj.hadamardSetting = HadamardSetting(m);
        end
        
        function [uHat, vHat] = decode(obj, y, sigma)
         %% Perform maximum likely hood decoding for RM(1,m)
         % rx: num_rx X N matrix
         % sigma: sigma^2 is the noise power
         %%
            [uHat, vHat] = HadamardDecode(obj.hadamardSetting, y);
        end
    end
    
end
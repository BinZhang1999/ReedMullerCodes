classdef DECODER_RM_AWGN_CXA
%% Generate Collapsed Punctured Aggregation decoder
% The algorithm is described in this paper: 
% Decoding Reed-Muller codes using redundant constrain
% by Mengke Lian ,Henry D. Pfister 
% published in ISIT 2020
%%

properties
   r;m; % Decoder for RM(r,m) codes
   K;N; % (N,K) code
   
   puncture_idx; % puncture_idx
   % it is generated from projection_idx
   % [num_subspace X num_overlap, num_coset] = size(puncture_idx);
   % [num_cn, 
end

methods
    
    
end
    
end
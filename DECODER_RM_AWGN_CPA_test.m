classdef DECODER_RM_AWGN_CPA_test
%% Generaete Collapsed Projection-Aggregation decoder
% The algorithm is described in this paper: 
% Decoding Reed-Muller codes using redundant constrain
% by Mengke Lian ,Henry D. Pfister 
% published in ISIT 2020
%% 
    properties
        r;m; % Decoder for RM(r,m) codes
        K;N; % (N,K) code
        
        projection_idx; % projection idx
        % All the (r-1)-dimensional subspaces of EG(m,2)
        projection_extrinsic_idx; % extrinsic node idx for each variable node
        num_iter; % iteration number 
        
        num_subspace; 
        num_coset;
        num_overlap;
        % [num_overlap, num_coset, num_subspace] = size(projection_idx)
        
        decoder_hadamard;
        sf; % a scaling factor
        decoder_scl;
    end
    methods
        function obj = Init(obj, r, m, projection_idx, ...
                projection_extrinsic_idx, num_iter, decoder_hadamard, sf, decoder_scl)
            obj.r = r;
            obj.m = m;
            obj.K=0;
            obj.N = 2^m;
            for r_i = 0:r
                obj.K = obj.K+nchoosek(m,r_i);
            end
            obj.projection_idx = projection_idx;
            obj.projection_extrinsic_idx = projection_extrinsic_idx;
            
            [obj.num_overlap, obj.num_coset, obj.num_subspace] = ...
                size(projection_idx);
            
            obj.num_iter = num_iter;
            
            obj.decoder_hadamard = decoder_hadamard;
            obj.sf = sf;
            
            obj.decoder_scl = decoder_scl;
        end
        
        
        function [u_hat, v_hat] = Decode(obj, rx, sigma)
            sigma_2 = sigma^2;
            llr = 2*rx / sigma_2;
            u_hat = zeros(1,obj.K);
            v_hat = zeros(1,obj.N);
            
                llr_x = llr(obj.projection_idx);
                llr_p = zeros(obj.num_coset, obj.num_subspace); 
                % Init projected llr
                llr_s = zeros(obj.N, obj.num_subspace);
                % Init the summation of the projected llr
                
                llr_p(:) = 2*atanh(prod(tanh(llr_x/2),1));
                llr_p = llr_p';
                
                [u_sub_hat,v_sub_hat] = obj.decoder_hadamard.Decode(llr_p); 
                
                idx = find(sum(u_sub_hat(:,1:6),2)==0);
                if isempty(idx)
                    return; % decoding error
                end
                
                me = sum( llr_p(idx,:).*(1-2*v_sub_hat(idx,:)), 2);
                [val, idx_max] = max(me);
                idx = idx(idx_max);
                
                llr(obj.projection_idx(1,:,idx))=(1-2*u_sub_hat(idx,end)).*...
                    llr(obj.projection_idx(1,:,idx));
                llr_p_sub = llr(obj.projection_idx(1,:,idx))+...
                    llr(obj.projection_idx(2,:,idx));
                
                [~, v_sub] = obj.decoder_scl.Decode(llr_p_sub, sqrt(2));
                
                v_hat(obj.projection_idx(1,:,idx)) = mod((v_sub+v_sub_hat(idx,:)),2);
                v_hat(obj.projection_idx(2,:,idx)) = v_sub;

        end
    
    end
end
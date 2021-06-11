classdef DECODER_RM_AWGN_CPA
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
    end
    methods
        function obj = Init(obj, r, m, projection_idx, ...
                projection_extrinsic_idx, num_iter, decoder_hadamard, sf)
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
        end
        
        
        function [u_hat, v_hat] = Decode(obj, rx, sigma)
            sigma_2 = sigma^2;
            llr = 2*rx / sigma_2;
            u_hat = zeros(1,obj.K);
            v_hat = zeros(1,obj.N);
            
            
            
            % start iteration decoding
            for iter = 1:obj.num_iter
                llr_x = llr(obj.projection_idx);
                llr_p = zeros(obj.num_coset, obj.num_subspace); 
                % Init projected llr
                llr_s = zeros(obj.N, obj.num_subspace);
                % Init the summation of the projected llr
                
                llr_p(:) = 2*atanh(prod(tanh(llr_x/2),1));
                llr_p = llr_p';
                
                [~,v_sub_hat] = obj.decoder_hadamard.Decode(llr_p); 
                change_vote = zeros(obj.num_subspace, obj.N);
                % The below codes is ugly but works.
                % Many thx for you, if you can provide better solution 
                % in the issue. 
                for i = 1:obj.num_subspace
                for j = 1:obj.num_overlap
                    change_vote(i, obj.projection_idx(j,:,i))=v_sub_hat(i,:);
                end
                end
                
                llr_z = llr(obj.projection_extrinsic_idx);
                llr_s(:) = 2*atanh(prod(tanh(llr_z/2),1));
                llr_s = llr_s';
                
                alpha = obj.sf*1/obj.num_subspace; % scaling factor
                llr_nxt = sum(alpha.*(1-2*change_vote).*llr_s, 1);
                
                llr_delta = abs(llr_nxt-llr);
                if ~any(llr_delta > 0.05)
                    break;
                end
                llr = llr_nxt;
            end
            v_hat(:) = llr<0;
        end
    
    end
end
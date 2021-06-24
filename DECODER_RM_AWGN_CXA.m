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
   % [num_cn, num_coset] = size(puncture_idx);
   num_cn;
   num_coset;
   
   num_iter;
   alpha;
   % parameters for map decoder
   code_dual;
   N_p;
end

methods
    function obj = Init(obj, r,m,puncture_idx, num_iter, alpha)
       obj.r = r;
       obj.m = m;
        obj.K=0;
        obj.N = 2^m;
        for r_i = 0:r
            obj.K = obj.K+nchoosek(m,r_i);
        end
       obj.puncture_idx = puncture_idx;
       [obj.num_cn, obj.num_coset] = size(puncture_idx);
       obj.num_iter = num_iter;
       obj.alpha = alpha;
       % parameters for map decoder
       obj.N_p = 2^(r+2);
       matrix_hadamard = hadamard(obj.N_p);
        obj.code_dual = zeros(obj.N_p, 2*obj.N_p);
        obj.code_dual(:,1:obj.N_p) = matrix_hadamard<0;
        obj.code_dual(:,obj.N_p+1:2*obj.N_p) = matrix_hadamard>0;
    end
    
       function llr_out = MAP_SISO(obj, rx_in, sigma)
        % the sigma_2 can be modified to lead different performance
        sigma_2 = sigma^2;
        rho = (1 - (exp(-(rx_in+1).^2/sigma_2 ) ./ exp(-(rx_in-1).^2/sigma_2 )) ) ./ ...
            (1 + (exp(-(rx_in+1).^2/sigma_2 ) ./ exp(-(rx_in-1).^2/sigma_2 )) );
        
        tau = log(abs(rho));
        
        beta_0 = (rho<0);
        
        s = tau * (1-2*obj.code_dual(:,1:obj.N_p));
        
        t_0 = (s(:,1)+s) / 2;
        t_1 = (s(:,1)-s) / 2;
        
        sigma_0 = mod(beta_0*(obj.code_dual(:,1:obj.N_p)==0), 2);
        % sigma_1 = mod(obj.N - sigma_0, 2);
        sigma_1 = mod(beta_0*(obj.code_dual(:,1:obj.N_p)==1), 2);
                                                                                                                                            
        h_0 = (1-2*sigma_0) .* exp(t_0);
        h_1 = (1-2*sigma_1) .* exp(t_1);
            
        x = h_0 - h_1;
        rj = x * (1-2*obj.code_dual(:,1:obj.N_p)); 
        % use rj represent the vector r 
        % cause we already used r to represent the order of RM code
        
        q = sum(h_0+h_1,2);
        
        omega_0 = (q-rj) / 2;
        omega_1 = (q+rj) / 2;
        
        z = (omega_0.*rho + omega_1./rho) ./ ...
            (omega_0 + omega_1);
        llr_out = z;

        app_0 = (z+1)/2;
        app_1 = (1-z)/2;
        llr_out = log2(app_0./app_1);
    end
    
    function [u_hat, v_hat] = Decode(obj, rx, sigma, v)
        v_hat = zeros(1,obj.N);
        u_hat = zeros(1,obj.K);
        % init lambda_v2c
        lambda_v2c = rx(obj.puncture_idx);
        % start iteration
        w = obj.alpha*1/obj.num_cn;
        for i = 1:obj.num_iter
            % perform map decoding for RM(r,r+2)
            lambda_c2v = obj.MAP_SISO(lambda_v2c, sigma);
            lambda_cvs = zeros(1, obj.N);
            for j = 1:obj.num_cn
                lambda_cvs(obj.puncture_idx(j,:)) =...
                    lambda_cvs(obj.puncture_idx(j,:))+ w*lambda_c2v(j,:);
            end
            %
            lambda_v2c = zeros(obj.num_cn, obj.N_p);
            for j = 1:obj.num_cn
                lambda_v2c(j,:) = ...
                    lambda_cvs(obj.puncture_idx(j,:))-w*lambda_c2v(j,:)+...
                    rx(obj.puncture_idx(j,:)); 
            end    
            
        end
        for j = 1:obj.num_cn
             v_hat(obj.puncture_idx(j,:)) =...
                 v_hat(obj.puncture_idx(j,:))+ w*lambda_c2v(j,:);
        end
        v_hat(:) = v_hat(:)<0;
    end
end
    
end
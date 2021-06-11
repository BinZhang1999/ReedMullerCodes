classdef DECODER_RM_AWGN_REED
%% Generaete Reed Decoder for RM(r,m) codes over AWGN channel
% decoder = DECODER_RM_AWGN_REED;
% decoder_rm37 = decoder.Init(3,7);
properties
    r;
    m; % RM(r,m)
    
    K; % code dimension
    N; % code length
    G; % Generaete matrix used for this decoder
    check_sum_idx; % index for the construction of check-sums 

end

methods
    function obj = Init(obj, r, m)
        % Init the basic parameters
        obj.r=r;
        obj.m=m;
        obj.K=0;
        obj.N = 2^m;
        for r_i = 0:r
            obj.K = obj.K+nchoosek(m,r_i);
        end
        % Init check sum
        obj.check_sum_idx = cell(1, obj.r);
        obj.G = zeros(obj.K,2^m); cnt_k = 0;
        
        variable = 1:obj.m;
        
        % generate the check-sums for each order monomials
        for r_i = obj.r:-1:1
            k_r = nchoosek(obj.m,r_i);
            subspace_list = zeros(k_r, 2^r_i);
            u_list = MyDec2Bin([0:2^r_i-1]', r_i);
            % treat the construction of check-sum as projection onto
            % r_i-demensional subspaces
            x = nchoosek(variable, r_i);
            subspace_init = zeros(r_i, obj.m);
            for k_i = 1:k_r
                % generate the r_i-dimensional subspace
               subspace = subspace_init;
               subspace(:,x(k_i,:)) = eye(r_i);
               subspace = subspace(:,end:-1:1);
               subspace = mod(u_list * subspace, 2);
               subspace_dec = sort(MyBin2Dec(subspace, obj.m),'ascend');
               subspace_list(k_i,:) = subspace_dec;
            end
            
            projection_idx = GenerateProjectionIdx(subspace_list, obj.m, r_i);
            obj.check_sum_idx{1,r_i} = projection_idx+1;
            
            for k_i = 1:k_r
                cnt_k = cnt_k+1;
                obj.G(cnt_k,projection_idx(end,:,k_i)+1) = 1;
            end
        end
        obj.G(end,:) = ones(1,2^m);
    end
    
    function [u_hat, v_hat] = Decode(obj, rx, sigma)
        sigma_2 = sigma^2;
        llr = 2*rx/sigma_2;
        u_hat = zeros(1,obj.K);
        v_hat = zeros(1,obj.N);
        cnt_k=0;
        for r_i = obj.r:-1:1
           k_i = nchoosek(obj.m,r_i);
           projection_idx = obj.check_sum_idx{1,r_i};
           llr_x = llr(projection_idx);
           llr_p = 2*atanh(prod(tanh(llr_x/2),1));
           llr_info = sum(llr_p,2);
           u_hat(cnt_k+1:cnt_k+k_i) = llr_info<0;

           llr = (1-2*...
               mod(u_hat(cnt_k+1:cnt_k+k_i)*obj.G(cnt_k+1:cnt_k+k_i,:),2) ...
                  ) .* llr;
           cnt_k = cnt_k+k_i;
        end
        u_hat(end) = (sum(llr)<0);
        v_hat(:) = mod(u_hat * obj.G,2);
    end
end
    
    
end
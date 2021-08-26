classdef DECODER_EHAMMING_AWGN_MAP
    % Decode RM(r,r+2) codes as extended hamming code using MAP decoding
    % the algorithm is described in the paper
    % Simple MAP Decoding of First-Order Reed¨CMuller and Hamming Codes
    % by Alexei Ashikhmin, Simon Litsyn
    % published in IEEE Trans. Info. Theory 
    properties

        m;
        K;
        N;
        code_dual; 
        % N X 2^(N-k) matrix 
        % each column is a codeword of the dual code
    end
    methods
        
    function obj = Init(obj, m)
        obj.m=m;
        obj.N = 2^m;
        matrix_hadamard = hadamard(obj.N);
        obj.code_dual = zeros(obj.N, 2*obj.N);
        obj.code_dual(:,1:obj.N) = matrix_hadamard<0;
        obj.code_dual(:,obj.N+1:2*obj.N) = matrix_hadamard>0;
        % the codewords of RM(1,m) that not containing 1 enumerate column
        % by column. 
        obj.K = obj.N-1-m;
    end
    
    function [u_hat, v_hat] = decode(obj, rx, sigma)
        % the sigma_2 can be modified to lead different performance
        sigma2 = sigma^2;
        rho = (1 - (exp(-(rx+1).^2/sigma2 ) ./ exp(-(rx-1).^2/sigma2 )) ) ./ ...
            (1 + (exp(-(rx+1).^2/sigma2 ) ./ exp(-(rx-1).^2/sigma2 )) );
        
        tau = log(abs(rho));
        
        beta0 = (rho<0);
        
        s = tau * (1-2*obj.code_dual(:,1:obj.N));
        
        t0 = (s(:,1)+s) / 2;
        t1 = (s(:,1)-s) / 2;
        
        sigma0 = mod(beta0*(obj.code_dual(:,1:obj.N)==0), 2);
        % sigma_1 = mod(obj.N - sigma_0, 2);
        sigma1 = mod(beta0*(obj.code_dual(:,1:obj.N)==1), 2);
                                                                                                                                            
        h0 = (1-2*sigma0) .* exp(t0);
        h1 = (1-2*sigma1) .* exp(t1);
            
        x = h0 - h1;
        rj = x * (1-2*obj.code_dual(:,1:obj.N)); 
        % use rj represent the vector r 
        % cause we already used r to represent the order of RM code
        
        q = sum(h0+h1,2);
        
        omega0 = (q-rj) / 2;
        omega1 = (q+rj) / 2;
        
        z = (omega0.*rho + omega1./rho) ./ ...
            (omega0 + omega1);
        v_hat = z<0;
        u_hat = zeros(1,obj.K); % do not offer the message about the info bits
    end
    
    end
    
    
end
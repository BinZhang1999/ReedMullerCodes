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
    
    function [u_hat, v_hat] = Decode(obj, rx, sigma)
        % the sigma_2 can be modified to lead different performance
        sigma_2 = sigma^2;
        rho = (1 - (exp(-(rx+1).^2/sigma_2 ) ./ exp(-(rx-1).^2/sigma_2 )) ) ./ ...
            (1 + (exp(-(rx+1).^2/sigma_2 ) ./ exp(-(rx-1).^2/sigma_2 )) );
        % sigma^2=2
        
        tau = log(abs(rho));
        
        beta_0 = (rho<0);
        
        s = tau * (1-2*obj.code_dual(:,1:obj.N));
        
        t_0 = (s(1)+s) / 2;
        t_1 = (s(1)-s) / 2;
        
        sigma_0 = mod(beta_0*(obj.code_dual(:,1:obj.N)==0), 2);
        % sigma_1 = mod(obj.N - sigma_0, 2);
        sigma_1 = mod(beta_0*(obj.code_dual(:,1:obj.N)==1), 2);
                                                                                                                                            
        h_0 = (1-2*sigma_0) .* exp(t_0);
        h_1 = (1-2*sigma_1) .* exp(t_1);
            
        x = h_0 - h_1;
        rj = x * (1-2*obj.code_dual(:,1:obj.N)); 
        % use rj represent the vector r 
        % cause we already used r to represent the order of RM code
        
        q = sum(h_0+h_1);
        
        omega_0 = (q-rj) / 2;
        omega_1 = (q+rj) / 2;
        
        z = (omega_0.*rho + omega_1./rho) ./ ...
            (omega_0 + omega_1);
        v_hat = z<0;
        u_hat = zeros(1,obj.K); % do not offer the message about the info bits
    end
    
    end
    
    
end
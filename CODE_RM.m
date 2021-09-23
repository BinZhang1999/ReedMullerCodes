classdef CODE_RM
%% Generate Binary Reed-Muller codes.
% eg: Generate RM(3,7):
% code_rm = CODE_RM;
% rm37 = code_rm.Init(3,7);
% u = (rand(1,rm37.K)>0);
% v = rm37.Encode(u);
properties
   r; % r order
   m; % m variable
   
   N; % Code length
   K; % Code dimension
   R; % Code rate
   
   
   Gm; % m kronecker-power matrix
   frozen_bits; % Frozen bits of RM codes when treat them as polar codes
   G; % Generator matrix
end

methods
    function obj = Init(obj,r,m)
        %% Init RM(r,m) codes
        obj.r = r; 
        obj.m = m;
        
        obj.N = 2^m;
        obj.K = 0;
        for i = 0:r
            obj.K = obj.K + nchoosek(m,i);
        end
        obj.R = obj.K / obj.N;
        
        F = [[1 0]
             [1 1]];
        obj.Gm = F;
        for i = 2:m
            obj.Gm = kron(F, obj.Gm);
        end
        obj.frozen_bits = (sum(obj.Gm,2) < 2^(m-r));
        obj.G = obj.Gm(~obj.frozen_bits,:);
    end
    
    function v = encode(obj, u)
        v = mod(u * obj.G, 2);
    end
    
end


    
    
end
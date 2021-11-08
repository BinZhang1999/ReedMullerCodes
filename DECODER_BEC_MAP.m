classdef DECODER_BEC_MAP
   properties
       H;
       N;
       K
   end
   methods
       function obj = Init(obj, H)
           obj.H = H;
           [R,obj.N]=size(H);
           obj.K = obj.N-R;
       end
       function [uHat, vHat] = decode(obj, r)
           uHat = nan(1,obj.K);
           vHat = nan(1,obj.N);
           isErasure = isnan(r);
           nErasure = sum(isErasure);
           
           if nErasure > obj.N-obj.K % cannot reconstruct
               return;
           end
           He = obj.H(:, isErasure);
           Hne = obj.H(:,~isErasure);
           rne = r(~isErasure);
           b = mod(Hne * rne',2);
           
           % we only need num_erasure rows
           % H_e * r_e' = H_ne * r_ne'
           [x, ~, nFree] = mySolveEquation(He, b);
           if nFree > 0
               return;
           end
           bHat = mod(He*x(1:nErasure),2);
           if any(bHat(1:nErasure)~=b(1:nErasure))
               1;
           end
           vHat(isErasure) = x(1:nErasure);
           vHat(~isErasure) = r(~isErasure);
       end
   end
end
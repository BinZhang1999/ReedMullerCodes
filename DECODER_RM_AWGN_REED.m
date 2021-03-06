classdef DECODER_RM_AWGN_REED
%% Generaete Reed Decoder for RM(r,m) codes over AWGN channel
% decoder = DECODER_RM_AWGN_REED;
% decoderRM37 = decoder.Init(3,7);
properties
    r;
    m; % RM(r,m)
    
    k; % code dimension
    n; % code length
    G; % Generaete matrix used for this decoder
    checkSumIndexCell; % index for the construction of check-sums 
end

methods
    function obj = Init(obj, r, m)
        % Init the basic parameters
        obj.r=r;
        obj.m=m;
        obj.k=0;
        obj.n = 2^m;
        for iOrder = 0:r
            obj.k = obj.k+nchoosek(m,iOrder);
        end
        % Init check sum
        obj.checkSumIndexCell = cell(1, obj.r);
        obj.G = zeros(obj.k,2^m); pk = 0;
        
        variable = 1:obj.m;
        % generate the check-sums for each order monomials
        for iOrder = obj.r : -1 : 1
            iOrderDim = nchoosek(obj.m, iOrder);
            iOrderSubspaceList = nan(iOrderDim, 2^iOrder);
            u_list = myDec2Bin([0:2^iOrder-1]', iOrder);
            % treat the construction of check-sum as projection onto
            % r_i-demensional subspaces
            iOrderMonomialList = nchoosek(variable, iOrder);
            for iDim = 1:iOrderDim
                % generate the r_i-dimensional subspace
               subspaceInBin = zeros(iOrder, obj.m);
               subspaceInBin(:,iOrderMonomialList(iDim,:)) = eye(iOrder);
               subspaceInBin = subspaceInBin(:,end:-1:1);
               subspaceInBin = mod(u_list * subspaceInBin, 2);
               subspaceInDec = sort(myBin2Dec(subspaceInBin, obj.m),'ascend');
               iOrderSubspaceList(iDim,:) = subspaceInDec;
            end
            
            iOrderCheckSumIndex = GenerateProjectionIdx(iOrderSubspaceList, obj.m, iOrder);
            obj.checkSumIndexCell{1,iOrder} = iOrderCheckSumIndex;
            
            for iDim = 1:iOrderDim
                pk = pk+1;
                obj.G(pk,iOrderCheckSumIndex(end,:,iDim)) = 1;
            end
        end
        obj.G(end,:) = ones(1,2^m);
    end
    
    function [uEsti, vEsti] = decode(obj, y, sigma)
        sigma2 = sigma^2;
        llr = 2*y/sigma2;
        uEsti = nan(1,obj.k);
        vEsti = nan(1,obj.n);
        pointK=0;
        llrForNextOrder = llr;
        for iOrder = obj.r:-1:1
           iK = nchoosek(obj.m,iOrder);
           iOrderCheckSumIndex = obj.checkSumIndexCell{1,iOrder};
           llrTemp = llrForNextOrder(iOrderCheckSumIndex);
           
           llrCheckSum = 2*atanh(prod(tanh(llrTemp/2),1));
           llrInfoBit = sum(llrCheckSum,2);
           uEsti(pointK+1:pointK+iK) = llrInfoBit<0;

           llrForNextOrder = (1-2*...
               mod(uEsti(pointK+1:pointK+iK)*obj.G(pointK+1:pointK+iK,:),2) ...
                  ) .* llrForNextOrder;
           pointK = pointK+iK;
        end
        uEsti(end) = (sum(llrForNextOrder)<0);
        vEsti(:) = mod(uEsti * obj.G,2);
    end
end
    
    
end
classdef DECODER_AWGN_OSD
    properties
        G;
        orderL;
        shiftIndexCell;
        listSizeSegment;
        listSize;
    end
    methods
        function obj = Init(obj,G,orderL)
            obj.G = G;
            [k,n]=size(G);
            obj.orderL = orderL;
            obj.listSize = 1;
            obj.listSizeSegment = nan(orderL, 1);
            obj.shiftIndexCell = cell(1,orderL);
            for iL = 1:orderL
                obj.listSizeSegment(iL) = nchoosek(k,iL);
                obj.listSize = obj.listSizeSegment(iL) + obj.listSize;
                obj.shiftIndexCell{1,iL} = nchoosek(1:k,iL);
            end
        end
        function [uEsti, vEsti] = decode(obj, y, sigma)
            %¡¡[uEsti,vEsti] = hardDecisionDecoding(obj.G,y);
            [uEsti, vEsti] = osdDecoding(obj.G, ...
                obj.orderL,...
                obj.shiftIndexCell, ...
                obj.listSizeSegment, ...
                obj.listSize, ...
                y);
        end
        
    end
    
end
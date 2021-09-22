classdef DECODER_AWGN_OSD
    properties
        G;
        dH;
        orderL;
        osdSetting;
        
    end
    methods
        function obj = Init(obj,G,dH,orderL)
            obj.G = G;
            obj.dH = dH;
            obj.orderL = orderL;
            obj.osdSetting = setOsdDecoder(G, dH, orderL);
        end
        function [uEsti, vEsti] = decode(obj, y, sigma)
            %¡¡[uEsti,vEsti] = hardDecisionDecoding(obj.G,y);
            [uEsti, vEsti] = osdDecoding(obj.osdSetting, y);
        end
        
    end
    
end
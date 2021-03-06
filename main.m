%% Code parameters
r=1;m=7;

%% Generate Code
code = CODE_RM;
code = code.Init(r,m);

%% Generate Decoder
decoder_select = 'Hadamard';
switch decoder_select
    case 'Reed'
        decoder = DECODER_RM_AWGN_REED;decoder = decoder.Init(r,m);
        
    case 'Hadamard' % only for RM(1,m) codes
        decoder = DECODER_RM_AWGN_HADAMARD;decoder = decoder.Init(m);
        
    case 'eHammingMAP' % only for RM(m-2,m) codes
        decoder = DECODER_EHAMMING_AWGN_MAP;decoder = decoder.Init(m);
        
    case 'CPA'
        % Load the projection_idx and projection_extrinsic_idx mannually
        % first.
        % Them can be generaeted using ToolsForProjectionDecoding.
        % I generate some xxx_idx and you can load from 
        % the folder ToolsForProjectionDecoding directly.
        decoder = DECODER_RM_AWGN_CPA; num_iter = 3;sf=2;
        decoder_hadamard = DECODER_RM_AWGN_HADAMARD;
        decoder_hadamard = decoder_hadamard.Init(m-r+1);
        decoder = decoder.Init(r,m,projection_idx,...
            projection_extrinsic_idx,num_iter, decoder_hadamard,sf);
        
    case 'CXA'
        decoder = DECODER_RM_AWGN_CXA; num_iter = 3;alpha = 3;
        decoder = decoder.Init(r,m,puncture_idx, num_iter,alpha);
    
    case 'OSD'
        decoder = DECODER_AWGN_OSD;
        decoder = decoder.Init(code.G, 2, 8); % osd-2
        
    otherwise
        disp(['Error: No matched decoder!']);
        return;
end

chase = 0; lise_size_log = 4;
if chase % if perform chase list decoding
    decoder_prototype = decoder;
    decoder = DECODER_CHASE; 
    decoder = decoder.Init(decoder_prototype, lise_size_log);
end

%% Simulation Settings
% over awgn channel
G = code.G;
simulationSetting.EbNoArray = 1.5:0.5:4.0;
simulationSetting.MIN_NUM_ERROR_FRAME = 100;
simulationSetting.displayName = '% OSD-2';
simulationSetting.description = '% OSD-2 Algorithm';
simulationResult = simulationAWGN(simulationSetting, G, decoder);


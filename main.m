%% Code parameters
r=1;m=7;

%% Generate Code
code = CODE_RM;
code = code.Init(r,m);

%% Generate Decoder
decoder_select = 'MAPDecoder';
switch decoder_select
    case 'ReedDecoder'
        decoder = DECODER_RM_AWGN_REED;decoder = decoder.Init(r,m);
    case 'HadamardDecoder' % only for 1-st order RM codes
        decoder = DECODER_RM_AWGN_HADAMARD;decoder = decoder.Init(m);
    otherwise
        disp(['No matched decoder']);
        exit(0);
end

chase = 1; lise_size_log = 4;
if chase % if perform chase list decoding
    decoder_prototype = decoder;
    decoder = DECODER_CHASE; 
    decoder = decoder.Init(decoder_prototype, lise_size_log);
end

%% Simulation Settings
num_least_error_frame = 100;

% over awgn channel
EbNo = 2:0.5:5.5;
sim = SIMULATION_AWGN;
sim = sim.Simulation(code, decoder, EbNo, num_least_error_frame);

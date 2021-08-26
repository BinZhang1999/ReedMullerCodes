function simulationResult = simulation(simulationSetting, G, decoder)

EbNoArray = simulationSetting.EbNoArray;
lengthEbNoArray = length(EbNoArray);
MIN_NUM_ERROR_FRAME = simulationSetting.MIN_NUM_ERROR_FRAME;
[k, n] = size(G);
R = k / n;
% Initialize simulation result
simulationResult.description = simulationSetting.description;
simulationResult.G = G;
simulationResult.decoder = decoder;
simulationResult.EbNoArray = simulationSetting.EbNoArray;
simulationResult.wer = nan(1, lengthEbNoArray);
simulationResult.ber = nan(1, lengthEbNoArray);
simulationResult.errorFrameMatrix = nan(MIN_NUM_ERROR_FRAME, n, lengthEbNoArray);
simulationResult.errorRecievedSymbolMatrix = nan(MIN_NUM_ERROR_FRAME, n, lengthEbNoArray);
simulationResult.errorMessageMatrix = nan(MIN_NUM_ERROR_FRAME, k, lengthEbNoArray);

for iEbNo = 1:lengthEbNoArray
    % Initialize simulation settings
    sigma = 1 / sqrt(2*R)*10^(-EbNoArray(iEbNo)/20);
    nErrorFrame = 0;
    nErrorBit = 0;
    nFrame = 0;
    tStart = tic;
    rng(1);
    % simulation at EbNoArray(iEbNo)
    while (nErrorFrame < MIN_NUM_ERROR_FRAME)
        % transimit codeword
        nFrame = nFrame + 1;
        u = (randn(1, k) > 0);
        v = mod(u * G, 2);
        bpskSymbol = 1 - 2*v;
        noise = randn(1, n);
        recievedSymbol = bpskSymbol + sigma.*noise;
        
        % recievedSymbol = awgn(bpskSymbol, EbNoArray(iEbNo));
        
        [uEsti, vEsti] = decoder.decode(recievedSymbol, sigma);
        % if decoding error
%        isError = ~any(ismember(uEsti(:,1:35),u(1:35),'row'));
%         isError = any(uEsti(:,1:35)~=u(1:35));
           
        isError = any(vEsti ~= v);
        if isError
            nErrorFrame = nErrorFrame + 1;
            nErrorBit = nErrorBit + sum(uEsti~=u);
            simulationResult.wer(iEbNo) = nErrorFrame / nFrame;
            simulationResult.ber(iEbNo) = nErrorBit / nFrame / k;
            
            simulationResult.errorFrameMatrix(nErrorFrame,:,iEbNo) = v;
            simulationResult.errorRecievedSymbolMatrix(nErrorFrame,:,iEbNo) = ...
                recievedSymbol;
            simulationResult.errorMessageMatrix(nErrorFrame,:,iEbNo) = u;
        end
        % print some simulation message
        isPrinting = isError;
        if isPrinting
            tSpan = toc(tStart);
            disp('############################################');
            disp(simulationSetting.description);
            disp(['Running time duration at this EbNo: ' num2str(tSpan) 's']);
            disp(['Error frame number at this EbNo: ' num2str(nErrorFrame)]);
            disp(['Eb/No = ' num2str(EbNoArray(iEbNo)) ' dB']);
            disp(['N = ' num2str(n) ' K = ' num2str(k)]);
            disp('EbNo    wer    ber');
            disp(num2str([simulationResult.EbNoArray' simulationResult.wer' ...
                simulationResult.ber']));
            disp('############################################');
        end
    end % end of the simulation at this frame
end % end of the simulation at this EbNo
end % end of the simulation function
classdef SIMULATION_AWGN
%% Simulation platform for AWGN channel
properties
    CODE;
    DECODER;
    
    num_least_error_frame;
    EbNo;
    wer;
    ber;
end

methods
    function obj = Simulation(obj, CODE, DECODER, EbNo, num_least_error_frame)
        G = CODE.G;
        obj.CODE = CODE;
        [K, N] = size(G);
        R = (K) / N;
        obj.DECODER = DECODER;
        obj.EbNo = EbNo;
        obj.num_least_error_frame = num_least_error_frame;
        obj.wer = zeros(1, length(EbNo));
        obj.ber = zeros(1, length(EbNo));
       
        for i = 1:length(EbNo)
            sigma = 1 / sqrt(2*R)*10^(-EbNo(i)/20);
            num_error_frame = 0;
            num_error_bits = 0;
            num_frame = 0;

            t_start = tic;
            while (num_error_frame < num_least_error_frame)
                u = randn(1, K) > 0.5;
                v = mod(u * G, 2);
                
                tx = 2 * (0.5-v);
                noise = randn(1, N);
                rx = sigma .* noise + tx;
                
                [u_hat, v_hat] = DECODER.Decode(rx, sigma);
                
                if any(v_hat~=v,'all')            
                    num_error_frame = num_error_frame + 1;
                    num_error_bits = num_error_bits+sum(u_hat~=u);
                end

               num_frame = num_frame + 1;

            end
            
            obj.wer(i) = num_error_frame / num_frame;
            obj.ber(i) = num_error_bits / (num_frame * K);

            t_end = toc(t_start);
            disp(['Running Time : ' num2str(t_end) 's']);
            fprintf("Eb/No = %.3f(dB) || wer = %.3g || ber = %.3g\n ",...
                EbNo(i), obj.wer(i), obj.ber(i));

        end
    end
end


    
    
end
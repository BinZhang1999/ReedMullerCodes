m = 9;
r = 4;

rm = CODE_RM;
rm = rm.Init(r,m);
d = 2^(m-r);
G = rm.G;
[k,n] = size(G);

rmDual = rm.Init(m-r-1,m);
H = rmDual.G;

maxRun = 4e5;
epsiArray = 0.05:0.05:1;
nP = zeros(1, length(epsiArray));

decoder = DECODER_BEC_MAP;
decoder = decoder.Init(H); 

tic
parfor iEpsi = 1:length(epsiArray)
    epsi = epsiArray(iEpsi);
    for i = 1:maxRun

        v = zeros(1,n);
        isErasure = logical(bsc(v, epsi));
        if ~isErasure(1) % 1st place not erasure
            continue;
        end
        nE = sum(isErasure);
        if nE < d % bounded distance decoder can even retrieve it
            continue;
        end
        
        % using H*v' = 0
        A = H(:,isErasure);
        b = zeros(n-k, 1);
        [xParticular, xNullspace, nFree] = mySolveEquation(A, b);
        if nFree == 0
            continue;
        end
        if any(xNullspace(1,:))
            nP(iEpsi) = nP(iEpsi)+1;
        end
        
        % using G'*u = v' 
%         A = G(:,~isErasure)';
%         b = zeros(n-nE, 1);
%         [xParticular, xNullspace, nFree] = mySolveEquation(A, b);
%         if nFree == 0 % there is a solution
%             continue;
%         end
%         
%         u = xNullspace';
%         v = mod(u * G,2);
%         if any(v(:,1))
%             nP(iEpsi) = nP(iEpsi)+1;
%         end
        
        
        
    end
end
toc
prm49 = nP./maxRun ./ epsiArray;
figure;
displayName = 'RM(4,9)';
plot(epsiArray, prm49, 'linewidth',1.5,'DisplayName',displayName);
ax = gca;
ax.FontWeight = 'bold';
% ax.GridLineStyle = '--';
% ax.MinorGridLineStyle = '--';
ax.LineWidth = 0.75;
ax.Box = 'on';
ax.GridAlpha = 0.3;
grid(ax,'on');
xticks(0:0.25:1); yticks(0:0.25:1);
legend;
axis([0 1 0 1]);
xlabel('Erasure Probability');
ylabel('Average EXIT function h');

n_k = 4;
[h,g,n,k] = hammgen(n_k);
h_enc = h(:,n_k+1:end);
H = sparse(h);

L2P = @(x) 1./(1+exp(x));

% ldpcEnc = comm.LDPCEncoder('ParityCheckMatrix',H);
ldpcDec = comm.LDPCDecoder('ParityCheckMatrix',H, 'DecisionMethod' , 'Soft decision');

qpskMod = comm.QPSKModulator('BitInput',true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true,...
    'DecisionMethod','Approximate log-likelihood ratio', ...
    'VarianceSource','Input port');
errorCnt = comm.ErrorRate;

noisebD = 3;
noiseVar = 1/10^(noisebD/10);
data = randi([0 1],k,1);
%% Hard Pass

hardData = logical(data);

% hard encode
encData = mod(h_enc * data, 2);
outData = [data; encData; 0];

% % channel
modSig = qpskMod(outData);
rxSig = awgn(modSig,noisebD);
demodSig = qpskDemod(rxSig,noiseVar);

% soft decode
rxL = ldpcDec(demodSig(1:end-1));
rxP = L2P(rxL)
% disp(errorCnt(data, rxData));

%% Soft Pass
softData = data + 0.1*randn(size(data));

% soft encode
for i=1:n_k
    encData(i) = softXorMod2(h_enc(i,:) .* softData');
end
outData = [softData; encData; 0];

% % channel
modSig = symbMap(outData);
rxSig = ChannelPass(modSig, [1], noiseVar);
demodSig = -2*(2*symbDemap(rxSig)-1);

% soft decode
rxL = ldpcDec(demodSig(1:end-1));
rxP = L2P(rxL)
% disp(errorCnt(data, rxData));

%%

for i=1:n_k
    encData(i) = softXorMod2(h_enc(i,:) .* rxP');
end
outData = [hardData; encData; 0];


%%
n_k = 4;
[h,g,n,k] = hammgen(n_k);

softIn = randn(k,1);
hardIn = round(softIn);

% mod(g' * softIn,2)
% softRes
for i=1:n
    softRes(i) = softXorMod2(softIn .* g(:,i));
end

for i=1:n_k
    softBack(i) = softXorMod2(softRes .* h(i,:));
end
softBack
% mod(g' * hardIn,2)' + softRes





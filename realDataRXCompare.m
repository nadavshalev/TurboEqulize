% close all;
% clear;
% clc;

%% settings

% load('./LiatModem/11-Mar-2019 09_12_39.mat');
% load('./LiatModem/11-Mar-2019 09_11_33.mat');
load('./LiatModem/11-Mar-2019 10_28_43.mat');


ecc.type = 'ldpc';

inputs.num_msg_enc = 64800;

ecc.n = inputs.num_msg_enc;
ecc.r = 0.75;
ecc.k = ecc.r * ecc.n;
inputs.num_msg_bits = ecc.k;

%% init

% set ECC
ecc.H = dvbs2ldpc(ecc.r);
ecc.encoder = comm.LDPCEncoder('ParityCheckMatrix',ecc.H);
ecc.decoder = comm.LDPCDecoder('ParityCheckMatrix',ecc.H, ...
                                       'DecisionMethod' , 'Soft decision', ...
                                       'OutputValue', 'Whole codeword');
                                   
padd_bits = record.TX.codeBits_zp(ecc.n:end);
                                   
% set Intrlv
intrlv.row = record.intrlv.row;
intrlv.col = record.intrlv.col;
intrlv.step = 3; % hstep is the slope of the diagonal

pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4);
L2P = @(x) 1./(1+exp(x));

%% comppadd_bits

% diff in size from tx to rx
cut_factor = (length(record.RX.res.dfe_out_hard) - length(record.TX.dataSymbls));
% cut the beggining of rx
dfe_out = record.RX.res.dfe_out_hard(cut_factor+1:end);

dmd = pskDemod(dfe_out);

bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
bits_dintrlv = bits_dintrlv(1:ecc.n);

decoded_llr = ecc.decoder(bits_dintrlv);
decoded_bits = L2P(decoded_llr);
Decoded = double(decoded_bits(1:ecc.k) > 0.5);

%% res

demapErr = sum(abs(dmd/2 - record.RX.res.after_demap(:)));
dintrlvErr = sum(abs(bits_dintrlv/2 - record.RX.res.after_deintlv(:)));
bitErr = sum(abs(Decoded - record.RX.res.res_bits(:)));

fprintf('map error: %e, dinterlive error: %e, bit error: %e\n', demapErr, dintrlvErr, bitErr)
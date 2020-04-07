close all;
clear;
clc;

%% settings

% load('./LiatModem/11-Mar-2019 09_12_39.mat');
% load('./LiatModem/11-Mar-2019 09_11_33.mat');
% load('./LiatModem/11-Mar-2019 10_28_43_c1.mat');
load('./LiatModem/11-Mar-2019 10_28_43_c2.mat');
% load('./LiatModem/11-Mar-2019 10_28_43_c3.mat');


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
                                   
% set Intrlv
intrlv.row = record.intrlv.row;
intrlv.col = record.intrlv.col;
intrlv.step = 3; % hstep is the slope of the diagonal

pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4,...
                                 'Variance',1);
L2P = @(x) 1./(1+exp(x));

%% comppadd_bits

% diff in size from tx to rx
cut_factor = (length(record.RX.res.dfe_out_hard) - length(record.TX.dataSymbls));
% cut the beggining of rx
dfe_out = record.RX.res.dfe_out_soft(cut_factor+1:end);

dmd = pskDemod(dfe_out);
dmdHard = (2*double(dmd/2>0)-1);

bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
bits_dintrlv = bits_dintrlv(1:ecc.n);
bits_dintrlvHard = (2*double(bits_dintrlv/2>0)-1);


decoded_llr = ecc.decoder(bits_dintrlv);
decoded_bits = L2P(decoded_llr);
Decoded = double(decoded_bits(1:ecc.k) > 0.5);

% res

demapErr = sum(abs(dmdHard - record.RX.res.after_demap(:)));
dintrlvErr = sum(abs(bits_dintrlvHard - record.RX.res.after_deintlv(:)));
bitErr = sum(abs(Decoded - record.RX.res.res_bits(:)));

fprintf('map error: %f, dinterlive error: %f, bit error: %d\n', demapErr, dintrlvErr, bitErr)

%% complete loop
padd_bits = record.TX.codeBits_zp(ecc.n+1:end)';
dfe_out_soft = record.RX.res.dfe_out_soft(cut_factor+1:end);






% [dn_, Decoded] = DecoderPath_real(dfe_out_soft, ecc, intrlv, pskDemod, padd_bits);


%     demodulate
    dmd = pskDemod(dfe_out_soft);
    
%     deinterlive
    bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
    bits_dintrlv = bits_dintrlv(1:ecc.n);
    
%     decode
    decoded_llr = ecc.decoder(bits_dintrlv);
    decoded_bits = L2P(decoded_llr);
    Decoded = double(decoded_bits(1:ecc.k) > 0.5);
    
    % re-encrypt
    decoded_bits_pd = [decoded_bits;padd_bits];
    bits_intrv = helscanintrlv(decoded_bits_pd,intrlv.row,intrlv.col,intrlv.step);
    dn_ = symbMap(bits_intrv);








loopbitErr = sum(abs(Decoded - record.TX.info_msg_with_CRC(:)));
loopSymbErr = sum(abs(dn_ - record.TX.dataSymbls(:)));

fprintf('CHAIN - BER: %d, SER: %d\n', loopbitErr, loopSymbErr)

figure;plt3(dn_)
hold on;plt3(dfe_out_soft)


%% Iter

% [dn_1, Decoded1] = DecoderPath_real(dfe_out_soft, ecc, intrlv, pskDemod, padd_bits);
% 
% loopbitErr = sum(abs(Decoded1 - record.RX.res.res_bits(:)));
% loopSymbErr = sum(abs(dn_1 - dfe_out));
% fprintf('CHAIN1 - BER: %d, SER: %d\n', loopbitErr, loopSymbErr)
% 
% [dn_2, Decoded2] = DecoderPath_real(dn_1, ecc, intrlv, pskDemod, padd_bits);
% 
% loopbitErr = sum(abs(Decoded2 - record.RX.res.res_bits(:)));
% loopSymbErr = sum(abs(dn_2 - dfe_out));
% fprintf('CHAIN2 - BER: %d, SER: %d\n', loopbitErr, loopSymbErr)
% 
% [dn_3, Decoded3] = DecoderPath_real(dn_2, ecc, intrlv, pskDemod, padd_bits);
% 
% loopbitErr = sum(abs(Decoded3 - record.RX.res.res_bits(:)));
% loopSymbErr = sum(abs(dn_3 - dfe_out));
% fprintf('CHAIN3 - BER: %d, SER: %d\n', loopbitErr, loopSymbErr)
% 
% 
% figure;hold on;
% plt3(dn_1);
% plt3(dn_2);
% plt3(dn_3);


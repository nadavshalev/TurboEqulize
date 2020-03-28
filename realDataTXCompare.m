close all;
clear;
% clc;

%%
% load('./LiatModem/rec1.mat');
load('./LiatModem/11-Mar-2019 09_12_39.mat');
% load('./LiatModem/11-Mar-2019 09_11_33.mat');
% load('./LiatModem/11-Mar-2019 10_28_43');

inputs.msg_bits = record.TX.info_msg_with_CRC(:);

%% settings

chann.overSamp = 2;

inputs.num_msg_enc = 64800;

ecc.type = 'ldpc';
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
inputs.enc_bits = ecc.encoder(inputs.msg_bits);


% inputs.enc_bits_zp = [inputs.enc_bits;zeros(record.TX.zp_len,1)];
padd_bits = record.TX.codeBits_zp(ecc.n:end);
inputs.enc_bits_zp = [inputs.enc_bits;padd_bits(:)];
inputs.enc_bits_zp = record.TX.codeBits_zp(:);

% set Intrlv
intrlv.row = record.intrlv.row;
intrlv.col = record.intrlv.col;
intrlv.step = 3; % hstep is the slope of the diagonal

% create symbs
% inputs.msg_symb = EncryptorPath2(inputs.enc_bits_zp, intrlv);
bits_intrv = helscanintrlv(inputs.enc_bits_zp,intrlv.row,intrlv.col,intrlv.step);
inputs.msg_symb = symbMap(bits_intrv);

%% compare

codedErr = sum(abs(inputs.enc_bits - record.TX.codeBits(:)));
intrlvErr = sum(abs(bits_intrv - record.TX.tx_bits(:)));
symbErr = sum(abs(inputs.msg_symb - record.TX.dataSymbls(:)));

fprintf('coded error: %e, interlive error: %e, symble error: %e\n', codedErr, intrlvErr, symbErr)

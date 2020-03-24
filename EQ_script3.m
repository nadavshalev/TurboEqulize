close all;
clear;
clc;

%% settings

% set params
h = [2-0.4j 1.5+1.8j 1 1.2-1.3j 0.8+1.6j];
nvar = 0.01; % noise var
overSamp = 2;

Ld = length(h);
Lr = overSamp*length(h);

ldpc.k = 48600;
ldpc.n = 64800;

inputs.num_msg_bits = ldpc.k;
inputs.num_train_bits = 2^13;

disp('========= Turbo Simulation =========');
fprintf('SNR: %f, chan L: %d\n', 10*log10(1/nvar), length(h));
fprintf('train num: %d, msg num %d (%f)\n',inputs.num_train_bits,inputs.num_msg_bits, inputs.num_train_bits/inputs.num_msg_bits);

%% init

% set LDPC
ldpc.r = ldpc.k/ldpc.n;
ldpc.H = dvbs2ldpc(ldpc.r);
ldpc.encoder = comm.LDPCEncoder('ParityCheckMatrix',ldpc.H);
ldpc.decoder = comm.LDPCDecoder('ParityCheckMatrix',ldpc.H, ...
                               'DecisionMethod' , 'Soft decision', ...
                               'OutputValue', 'Whole codeword');
% set Intrlv
intrlv.row = 10;
intrlv.col = ldpc.n/intrlv.row;
intrlv.step = 3; % hstep is the slope of the diagonal

% create bits
inputs.msg_bits = randi([0 1], 1, inputs.num_msg_bits)';
inputs.train_bits = randi([0 1], 1, inputs.num_train_bits)';
inputs.enc_bits = ldpc.encoder(inputs.msg_bits);

% create symbs
inputs.msg_symb = EncryptorPath2(inputs.enc_bits, intrlv);
inputs.train_symb = symbMap(inputs.train_bits);
inputs.num_msg_symb = length(inputs.msg_symb);
inputs.num_train_symb = length(inputs.train_symb);

errorCnt = comm.ErrorRate;
L2P = @(x) 1./(1+exp(-x));

%% channel

sn = inputs.msg_symb + 0.01*randn(size(inputs.msg_symb));

bits = symbDemap2(sn);
bits_dintrlv = helscandeintrlv(bits,intrlv.row,intrlv.col,intrlv.step);
decoded_llr = ldpc.decoder(bits_dintrlv);
decoded_bits = double(L2P(decoded_llr) > 0.5);

e = errorCnt(decoded_bits(1:inputs.num_msg_bits), inputs.msg_bits);
fprintf('ldpc errors: %d (%f)\n', e(2), e(2) / e(3));

hardBits = double(bits_dintrlv(1:inputs.num_msg_bits)>0);
e = errorCnt(hardBits(1:inputs.num_msg_bits), inputs.msg_bits);
fprintf('ref  errors: %d (%f)\n', e(2), e(2) / e(3));
% % oversample
% os_el = ones(overSamp,1);
% symb_chan_pre = kron([symb_in_train;symb_in],os_el);
% 
% % channel
% channel_out = ChannelPass(symb_chan_pre, h, nvar);
% symb_chan_train = channel_out(1:overSamp*num_symb_train+(Ld-1)); % train symbs after channel with memory of L-1
% symb_chan = channel_out(overSamp*num_symb_train+1:end);



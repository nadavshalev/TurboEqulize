% close all;
clear;
% clc;

% channel
% h = [2-0.4j 1.5+1.8j 1 1.2-1.3j 0.8+1.6j];
% h = [2-0.4j 1.5+1.8j 1];
h = [1 0 0 0 0 0];

mBits = 2; % 2 bit/sym => QPSK
repNum = 3; % ECC
nvar = 0.01; % noise var
SNR = 10*log10(1/nvar);
overSamp = 2;

Ld = length(h);
Lr = overSamp*length(h);

% create bits
num_bits = 2^13; % data bit len
num_train_bits = 2^10; % train bit len
bits = randi([0 1], 1, num_bits);
trBits = randi([0 1], 1, num_train_bits);

% create symbols
num_symb = num_bits*repNum/mBits; % num data symbls
num_symb_train = num_train_bits/mBits; % num train symbls
symb_in = EncryptorPath(bits, repNum); % symbles after TX chain (before modulation)
symb_in_train = symbMap(trBits);

%% simulation
% oversample
os_el = ones(overSamp,1);
symb_chan_pre = kron([symb_in_train;symb_in],os_el);

% channel
channel_out = ChannelPass(symb_chan_pre, h, nvar);
symb_chan_train = channel_out(1:overSamp*num_symb_train+(Ld-1)); % train symbs after channel with memory of L-1
symb_chan = channel_out(overSamp*num_symb_train+1:end);

%% Turbo

% set EQ
mu = 0.001; % update step size
EQ_turbo = AdaEQ(Lr, Ld,mu, overSamp); % set
hard = @(x) getHard(x);
maxiter = 10;

% train equlizer
symb_eq_train = EQ_turbo.turboEqualize_train(symb_chan_train,symb_in_train, num_symb_train);
err = calcError(symb_eq_train, symb_in_train, hard);
disp(['Turbo - train BER: ' num2str(err) '  MSE:' num2str(mean(abs(symb_eq_train - symb_in_train)))]);

% set params
MSE_eq = zeros(maxiter,1);
err_eq = zeros(maxiter,1);
SE_bit = zeros((num_symb)*maxiter,1);
pre_in_symb = symb_in_train(end-Ld+1:end);
pre_chan_symb = symb_chan_train(end-(Ld-1)-Lr+1:end-(Ld-1));
symb_chan_input = [pre_chan_symb;symb_chan]; % add L smples for time channel response 
inds = (1:num_symb);

% equlize data
figure;hold on;
for i = 1:maxiter
    if i == 1 % have no dn_ yet => run simple eq
        symb_eq = EQ_turbo.normalEqualize_run(symb_chan_input, num_symb);
    else
        symb_dn_input = [pre_in_symb;dn_]; % add L smples for time channel response 
        symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input, num_symb);
    end
    [dn_, Decoded] = DecoderPath(symb_eq, repNum); % decode and re-encode
    
    %print and plot
    err_eq(i) = calcError(Decoded, bits);
    MSE_eq(i) = mean(abs(symb_eq - symb_in).^2);
    SE_bit(inds) = abs(symb_eq - symb_in).^2;
    plot(inds,SE_bit(inds)); drawnow;
    fprintf(['Turbo - iteration ' num2str(i) ': BER: ' num2str(err_eq(i)) '  MSE:' num2str(MSE_eq(i))]); fprintf('\n');
    inds = inds + num_symb;
end

Stat = EQ_turbo.Statistics;
xlabel("Filters Update number"); ylabel("bit SE");
title("bit MSE on each filters update"); grid on;
hold off;

figure;semilogy(smooth(SE_bit,1000)); grid on;


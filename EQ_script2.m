close all;
clear;clc;

% channel
h = [2-0.4j 1.5+1.8j 1 1.2-1.3j 0.8+1.6j];
% h = [0.9 0.1];
L = length(h);

mBits = 2; % 2 bit/sym => QPSK
repNum = 1; % ECC
nvar = 0.05; % noise var

% create bits
num_bits = 2^15; % data bit len
num_train_bits = 2^13; % train bit len
bits = randi([0 1], 1, num_bits);
trBits = randi([0 1], 1, num_train_bits);

% create symbols
num_symb = num_bits*repNum/mBits; % num data symbls
num_symb_train = num_train_bits/mBits; % num train symbls
symb_in = EncryptorPath(bits, repNum); % symbles after TX chain (before modulation)
symb_in_train = symbMap(trBits);

%%
% channel
channel_out = ChannelPass([symb_in_train;symb_in], h, nvar);
symb_chan_train = channel_out(1:num_symb_train+(L-1)); % train symbs after channel with memory of L-1
symb_chan = channel_out(num_symb_train+1:end);

%% Turbo

% set EQ
mu = 0.001; % update step size
EQ_turbo = AdaEQ(L,mu); % set
hard = @(x) getHard(x);
maxiter = 10;

% train equlizer
symb_eq_train = EQ_turbo.turboEqualize_train(symb_chan_train,symb_in_train);
err = calcError(symb_eq_train, symb_in_train, hard);
disp(['Turbo - train BER: ' num2str(err) '  MSE:' num2str(mean(abs(symb_eq_train - symb_in_train)))]);

% set params
MSE_eq = zeros(maxiter,1);
err_eq = zeros(maxiter,1);
SE_bit = zeros((num_symb)*maxiter,1);
pre_in_symb = symb_in_train(end-L+1:end);
pre_chan_symb = symb_chan_train(end-(L-1)-L+1:end-(L-1));
symb_chan_input = [pre_chan_symb;symb_chan]; % add L smples for time channel response 
inds = (1:num_symb);

% equlize data
figure;hold on;
for i = 1:maxiter
    if i == 1 % have no dn_ yet => run simple eq
        symb_eq = EQ_turbo.normalEqualize_run(symb_chan_input);
    else
        symb_dn_input = [pre_in_symb;dn_]; % add L smples for time channel response 
        symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input);
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

figure;semilogy(smooth(SE_bit,1000))


%% DFE

% set EQ
mu = 0.001; % update step size
EQ_dfe = AdaEQ(L,mu); % set
hard = @(x) getHard(x);
maxiter = 10;

% train equlizer
symb_eq_train = EQ_dfe.turboEqualize_train(symb_chan_train,symb_in_train);
err = calcError(symb_eq_train, symb_in_train, hard);
disp(['Turbo - train BER: ' num2str(err) '  MSE:' num2str(mean(abs(symb_eq_train - symb_in_train)))]);

% set params
MSE_eq = zeros(maxiter,1);
err_eq = zeros(maxiter,1);
SE_bit = zeros((num_symb)*maxiter,1);
pre_in_symb = symb_in_train(end-L+1:end);
pre_chan_symb = symb_chan_train(end-(L-1)-L+1:end-(L-1));
symb_chan_input = [pre_chan_symb;symb_chan]; % add L smples for time channel response 
inds = (1:num_symb);
symb_eq = [pre_chan_symb;symb_chan];
% equlize data
figure;hold on;
for i = 1:maxiter
    symb_eq = EQ_dfe.dfeEqualize(symb_chan_input);
    [dn_, Decoded] = DecoderPath(symb_eq, repNum); % decode and re-encode
    
%     symb_eq = [pre_chan_symb;symb_eq];
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

figure;semilogy(smooth(SE_bit,1000))























%%

% figure;plot3(1:length(Stat.symbHard), real(Stat.symbHard), imag(Stat.symbHard),'.'); hold on;
% plot3(1:length(Stat.sn), real(Stat.sn), imag(Stat.sn),'.')
% ca = Stat.sn;
% ca(abs(Stat.sn)<1.05*sqrt(2)) = -inf;
% plot3(1:length(ca), real(ca), imag(ca),'.k')
% figure;plot3(1:length(Stat.E), 10*real(Stat.E), 10*imag(Stat.E),'.')

% figure;plot3(1:length(Stat.symbHard), real(Stat.symbHard), imag(Stat.symbHard),'.'); hold on;
% plot3(1:length(Stat.sn), real(Stat.sn), imag(Stat.sn),'.')

%% Normal EQ

% set EQ
EQ_dfe = AdaEQ(L,0.001); % set

% train equlizer
symb_eq_train = EQ_dfe.normalEqualize_train(symb_chan_train,symb_in_train, 0.001);
err = calcError(symb_eq_train, symb_in_train, hard);
disp(['normal - train error: ' num2str(err)]);

pre_in_symb = symb_in_train(end-L+1:end);
pre_chan_symb = symb_chan_train(end-(L-1)-L+1:end-(L-1));

symb_dn_input = [pre_in_symb;symb_chan];
inds = (1:num_symb);

% set params
maxiter = 10; 
MSE_nr = zeros(maxiter,1);
err_nr = zeros(maxiter,1);

% equlize data
for i = 1:maxiter
    symb_eq = EQ_dfe.normalEqualize_run(symb_dn_input, EQ_dfe.Mu);

    symb_dn_input = [pre_in_symb;symb_eq];
    
    [~, Decoded] = DecoderPath(symb_eq, repNum);
    
    err_nr(i) = calcError(Decoded, bits);
    MSE_nr(i) = mean(abs(symb_eq - symb_in).^2);
    
    disp(['normal - ieration ' num2str(i) ': BER: ' num2str(err_nr(i)) '  MSE:' num2str(MSE_nr(i))]);
    fprintf('\n')
%     if err == 0
%         break;
%     end
end
%%
% figure(2);semilogy(MSE_eq);hold on;
% semilogy(MSE_nr);hold off;
% xlabel("Number of Turbo iteration"); ylabel("MSE");
% title("MSE"); grid on;
% 
% figure(3);hold on;
% plot(err_eq*100)
% plot(err_nr*100);hold off;
% xlabel("Number of Turbo iteration"); ylabel("bit error rate [%]");
% title("bit error rate"); grid on;
% 






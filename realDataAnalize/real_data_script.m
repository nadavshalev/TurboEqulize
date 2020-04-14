close all;
clear;
clc;

%% cahnnels option (uncomment dsierd channel)

% load('./11-Mar-2019 15_25_09_c1.mat'); % not working
% load('./11-Mar-2019 15_25_09_c2.mat'); % one iter
% load('./11-Mar-2019 15_25_09_c3.mat'); 

% load('./11-Mar-2019 15_24_14_c1.mat'); % not working
load('./11-Mar-2019 15_24_14_c2.mat');
% load('./11-Mar-2019 15_24_14_c3.mat');

% load('./11-Mar-2019 15_22_22_c1.mat');  % not working
% load('./11-Mar-2019 15_22_22_c2.mat');
% load('./11-Mar-2019 15_22_22_c3.mat');

%% init

% turbo hyper-params
turbo.Ld = 40; % BF tap num
turbo.Lr = 10; % FF tap num
turbo.mu = 0.004; % step size
turbo.iter = 5; % number of turbo iterations

% set ECC
ecc.type = 'ldpc';
ecc.n = 64800;
ecc.r = 0.75;
ecc.k = ecc.r * ecc.n;
ecc.H = dvbs2ldpc(ecc.r);
ecc.encoder = comm.LDPCEncoder('ParityCheckMatrix',ecc.H);
ecc.decoder = comm.LDPCDecoder('ParityCheckMatrix',ecc.H, ...
                                       'DecisionMethod' , 'Soft decision', ...
                                       'OutputValue', 'Whole codeword');    
                                   
                              
% set Intrlv
intrlv.row = record.intrlv.row;
intrlv.col = record.intrlv.col;
intrlv.step = 3; % hstep is the slope of the diagonal

% set demodulation
pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4,...
                                 'Variance',1);
% LLR to probability
L2P = @(x) 1./(1+exp(x));


% set input params
inputs.num_msg_enc = ecc.n;

inputs.msg_bits = record.TX.info_msg_with_CRC(:);
inputs.num_msg_bits = ecc.k;

inputs.train_symb = record.RX.train_symb(:);
inputs.num_train_symb = length(inputs.train_symb);

inputs.msg_symb = record.TX.dataSymbls(:);
inputs.num_msg_symb = length(record.TX.dataSymbls);

inputs.padd_bits = record.TX.codeBits_zp(ecc.n+1:end)';

% set channel
chann.overSamp = 2;

chann.out = record.RX.chann_out(:);
chann.num_out = length(chann.out);

chann.train_symb = chann.out(1:chann.overSamp*inputs.num_train_symb+(turbo.Ld-1)); % train symbs after channel with memory of L-1
chann.msg_symb = chann.out(chann.overSamp*inputs.num_train_symb+1:end);
chann.num_train = length(chann.train_symb);
chann.num_msg = length(chann.msg_symb);

%% Turbo

% set EQ object
EQ_turbo = AdaEQ(turbo.Lr, turbo.Ld, turbo.mu, chann.overSamp, 200); % set
hard = @(x) getHard(x);

% train equlizer
eq.train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
% train results
eq.train_err = calcError(eq.train_symb,inputs.train_symb, hard);
eq.train_mse = mean(abs(eq.train_symb-inputs.train_symb).^2);
eq.train_papr = max(abs(eq.train_symb)) / mean(abs(eq.train_symb));
fprintf('Train: \t\t\tBER: %f \realDataChannelPlot(chann, inputs);tMSE:%f \tPAPR:%f -----\n', eq.train_err, eq.train_mse, eq.train_papr);
% plot train results
x = abs(eq.train_symb-inputs.train_symb).^2;
figure(1);plot(x);drawnow;
grid('on')
xlabel('Symbol Number')
ylabel('Square Error')
title('Equalization Train Error')
hold off

% set params
msg_pre_in_symb = inputs.train_symb(end-turbo.Ld+1:end);
msg_pre_chan_symb = chann.train_symb(end-(turbo.Ld-1)-turbo.Lr+1:end-(turbo.Ld-1));
symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
inds = (1:inputs.num_msg_symb);

% equlize data
figure(2);hold on;
for i = 1:turbo.iter+1
    if i == 1 % have no dn_ yet => run simple eq
        symb_eq = record.RX.res.dfe_out_soft(inputs.num_train_symb+1:end);
    else
        symb_dn_input = [msg_pre_in_symb;dn_]; % add L smples for time channel response 
        symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input, inputs.num_msg_symb);
    end
    
    [dn_, Decoded] = DecoderPath_real(symb_eq, ecc, intrlv, pskDemod, inputs.padd_bits);
    
    % print and plot
    eq.err_eq(i) = calcError(Decoded, inputs.msg_bits);
    eq.msg_mse(i) = mean(abs(symb_eq - inputs.msg_symb).^2);
    eq.papr(i) = max(abs(symb_eq)) / mean(abs(symb_eq));
    SE_bit(inds) = abs(symb_eq - inputs.msg_symb).^2;
    figure(2);plot(inds,SE_bit(inds)); drawnow;
	fprintf('Turbo - iteration %d: \tBER: %f \tMSE:%f \tPAPR:%f\n', i, eq.err_eq(i), eq.msg_mse(i), eq.papr(i));
    
    inds = inds + inputs.num_msg_symb;
    
end
grid('on')
xlabel('Symbol Number')
ylabel('Square Error')
legend('DFE', 'iter1', 'iter2', 'iter3', 'iter4', 'iter5')
title('Equalization Error')

figure(3);stem(eq.err_eq)
hold on; stem(eq.err_eq(1))
grid('on')
xlabel('Iteration')
ylabel('BER')
title('Equalization BER')
legend('Turbo', 'DFE')

%% chann est
realDataChannelPlot(chann, inputs);
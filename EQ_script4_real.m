close all;
clear;
clc;

%% settings

dfe.Lr = 14;
dfe.Ld = 4;
dfe.mu = 0.14;
dfe.type = 'RLS';
turbo.type = 'LIAT';

% hard channels
turbo.Ld = 40;
turbo.Lr = 10;
turbo.mu = 0.004;

% % % easy channels
% turbo.Ld = 1;
% turbo.Lr = 8;
% turbo.mu = 0.08;

%% 43
% load('./LiatModem/11-Mar-2019 10_28_43_c1.mat');
% turbo.Ld = 1;
% turbo.Lr = 7;
% turbo.mu = 0.03;
% turbo.type = 'LIAT';
% turbo.iter = 24;

% load('./LiatModem/11-Mar-2019 10_28_43_c2.mat');
% turbo.Ld = 1;
% turbo.Lr = 5;
% turbo.mu = 0.08;
% turbo.type = 'MATLAB';
% turbo.iter = 15;
% dfe.Lr = 11;
% dfe.Ld = 1;
% dfe.mu = 0.14;
% dfe.type = 'LMS';

% load('./LiatModem/11-Mar-2019 10_28_43_c3.mat');
% turbo.Ld = 1;
% turbo.Lr = 4;
% turbo.mu = 0.1;
% turbo.type = 'MATLAB';
% turbo.iter = 2;
% dfe.Lr = 5;
% dfe.Ld = 1;
% dfe.mu = 0.14;
% dfe.type = 'LMS';

%% 33
% load('./LiatModem/11-Mar-2019 09_11_33_c1.mat');
% turbo.Ld = 1;
% turbo.Lr = 8;
% turbo.mu = 0.08;
% turbo.type = 'LIAT';
% turbo.iter = 2;

% % NOT COMB
% load('./LiatModem/11-Mar-2019 09_11_33_c2.mat');
% turbo.Ld = 1;
% turbo.Lr = 6;
% turbo.mu = 0.11;
% turbo.type = 'LIAT';
% turbo.iter = 2;

% load('./LiatModem/11-Mar-2019 09_11_33_c3.mat');
% turbo.Ld = 1;
% turbo.Lr = 8;
% turbo.mu = 0.08;
% turbo.type = 'LIAT';
% turbo.iter = 2;

%% 39

% load('./LiatModem/11-Mar-2019 09_12_39_c1.mat');
% turbo.Ld = 1;
% turbo.Lr = 8;
% turbo.mu = 0.08;
% turbo.type = 'LIAT';

% % NOT WORKING
% load('./LiatModem/11-Mar-2019 09_12_39_c2.mat');
% turbo.Ld = 2;
% turbo.Lr = 15;
% turbo.mu = 0.11;
% turbo.type = 'LIAT';

% load('./LiatModem/11-Mar-2019 09_12_39_c3.mat');
% turbo.Ld = 2;
% turbo.Lr = 15;
% turbo.mu = 0.11;
% turbo.type = 'LIAT';
turbo.iter = 5;
%% HARD --------------------------
%% 09

% load('./LiatModem/11-Mar-2019 15_25_09_c1.mat');
% turbo.Ld = 30;
% turbo.Lr = 6;
% turbo.mu = 0.0072;
% turbo.type = 'LIAT';
% turbo.iter = 2;

% good example
% load('./LiatModem/11-Mar-2019 15_25_09_c2.mat');
% turbo.Ld = 30;
% turbo.Lr = 3;
% turbo.mu = 0.005;
% turbo.type = 'LIAT';
% turbo.iter = 2;

% % good example
% load('./LiatModem/11-Mar-2019 15_25_09_c3.mat');
% turbo.Ld = 30;
% turbo.Lr = 3;
% turbo.mu = 0.005;
% turbo.type = 'LIAT';
% turbo.iter = 2;
%% 14

% % NOT WORKING!
% load('./LiatModem/11-Mar-2019 15_24_14_c1.mat');
% turbo.Ld = 40;
% turbo.Lr = 10;
% turbo.mu = 0.004;
% turbo.type = 'LIAT';
% turbo.iter = 10;

% BEST example + graph
% BEST params graph
% load('./LiatModem/11-Mar-2019 15_24_14_c2.mat');
% turbo.Ld = 40;
% turbo.Lr = 10;
% turbo.mu = 0.004;
% turbo.type = 'LIAT';
% turbo.iter = 13;

% % BEST example + graph
% load('./LiatModem/11-Mar-2019 15_24_14_c3.mat');
% turbo.Ld = 40;
% turbo.Lr = 10;
% turbo.mu = 0.004;
% turbo.type = 'LIAT';
% turbo.iter = 10;
%% 22

% % NOT COMB 
% load('./LiatModem/11-Mar-2019 15_22_22_c1.mat');
% turbo.Ld = 13;
% turbo.Lr = 6;
% turbo.mu = 0.008;
% turbo.type = 'LIAT';
% turbo.iter = 2;

% % good example
% load('./LiatModem/11-Mar-2019 15_22_22_c2.mat');
% turbo.Ld = 40;
% turbo.Lr = 6;
% turbo.mu = 0.0065;
% turbo.type = 'LIAT';
% turbo.iter = 4;

% good example
% load('./LiatModem/11-Mar-2019 15_22_22_c3.mat');
% turbo.Ld = 40;
% turbo.Lr = 6;
% turbo.mu = 0.0065;
% turbo.type = 'LIAT';
% turbo.iter = 3;
%% init

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

pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4,...
                                 'Variance',1);
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


% set matlab DFE
Matlab_DFE = comm.DecisionFeedbackEqualizer('Algorithm',dfe.type, ...
    'NumForwardTaps',dfe.Lr,'NumFeedbackTaps',dfe.Ld,...
    'Constellation',[1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2),...
    'ReferenceTap',1,'InputSamplesPerSymbol',2,...
    'StepSize', dfe.mu);


% realDataChannelPlot(chann, inputs);
%%


% turbo.iter = 5;
% turbo.Ld = 17;
% turbo.Lr = 11;
% turbo.mu = 0.008;

% set EQ
EQ_turbo = AdaEQ(turbo.Lr, turbo.Ld, turbo.mu, chann.overSamp, 200); % set
hard = @(x) getHard(x);

% train equlizer
eq.train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
% train results
eq.train_err = calcError(eq.train_symb,inputs.train_symb, hard);
eq.train_mse = mean(abs(eq.train_symb-inputs.train_symb).^2);
eq.train_papr = max(abs(eq.train_symb)) / mean(abs(eq.train_symb));
fprintf('Train: \t\t\tBER: %f \realDataChannelPlot(chann, inputs);tMSE:%f \tPAPR:%f -----\n', eq.train_err, eq.train_mse, eq.train_papr);

x = abs(eq.train_symb-inputs.train_symb).^2;
figure;plot(x);
grid('on')
xlabel('Symbol Number')
ylabel('Square Error')
title('Equalization Train Error')

% set params
msg_pre_in_symb = inputs.train_symb(end-turbo.Ld+1:end);
msg_pre_chan_symb = chann.train_symb(end-(turbo.Ld-1)-turbo.Lr+1:end-(turbo.Ld-1));
symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
inds = (1:inputs.num_msg_symb);

% equlize data
figure;hold on;
for i = 1:turbo.iter
    if i == 1 % have no dn_ yet => run simple eq
        if strcmp(turbo.type,'LIAT')
            symb_eq = record.RX.res.dfe_out_soft(inputs.num_train_symb+1:end);
        elseif strcmp(turbo.type,'MATLAB')
            [eq_matlab,~,~] = Matlab_DFE(chann.out ,inputs.train_symb);
            symb_eq = eq_matlab(inputs.num_train_symb+1:end);
        end
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
    plot(inds,SE_bit(inds)); drawnow;
	fprintf('Turbo - iteration %d: \tBER: %f \tMSE:%f \tPAPR:%f\n', i, eq.err_eq(i), eq.msg_mse(i), eq.papr(i));
    
    inds = inds + inputs.num_msg_symb;
    
end
grid('on')
xlabel('Symbol Number')
ylabel('Square Error')
% legend('DFE', 'iter1', 'iter2', 'iter3', 'iter4', 'iter5', 'iter6')
title('Equalization Error')


% figure;semilogy(smooth(SE_bit,1000)); grid on;
figure;stem(eq.err_eq)
hold on; stem(eq.err_eq(1))
grid('on')
xlabel('Iteration')
ylabel('BER')
title('Equalization BER')
legend('Turbo', 'DFE')
%%
% LrArr = 6:16;
% LdArr = 1:5;
% muArr = 0.01:0.03:0.12;
% N = length(LrArr) * length(LdArr) * length(muArr)
% % LrArr = 5:7;
% % LdArr = 1:2;
% % muArr = 0.1:0.005:0.13;
% [resData, minData] = searchEQParams(dn_, chann, inputs, muArr, LrArr, LdArr);

% LrArr = 2:20;
% LdArr = 1:15;
% muArr = 0.002:0.001:0.02;
% N = LrArr * LdArr * muArr;
% % LrArr = 4:12;
% % LdArr = 1:4;
% % muArr = 0.09:0.01:0.1;
% [resData, minData] = searchEQParamsMATLAB(chann, inputs, LrArr, LdArr, muArr);
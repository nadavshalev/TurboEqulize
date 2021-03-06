 close all;
 clear;
 clc;

%% settings

% set params
% chanels

Porat_and_al = [2-0.4j,1.5+1.8j,1,1.2-1.3j,0.8+1.6j];
Proakis_A = [0.04,-0.05,0.07,-0.21,-0.5,0.72,0.36,0,0.21,0.03,0.07];
Proakis_B = [0.407,0.815,0.407];
Matlab_Channel = [0.227,0.460,0.688,0.460,0.227];
Matlab_Channel2 =[.986,.845,.237,.12345+.31j,0];
Matlab_Channel3 =[0.623,0.489+0.234j,0.398i, 0.21,0];
My_Channel = [1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256];

chann.h = Porat_and_al;

% L = 10;
% chann.h = randi([-2 2],1,L) + 1j*randi([-2 2],1,L);


chann.SNR = 10; % noise var

chann.nVar = 10 ^ (-chann.SNR/10);
chann.overSamp = 1;

Ld = length(chann.h);
Lr = chann.overSamp*length(chann.h);

ecc.type = 'ldpc';

inputs.num_train_bits = 2^11;

inputs.num_msg_enc = 64800;
switch ecc.type
    case 'ldpc'
        ecc.n = inputs.num_msg_enc;
        ecc.r = 0.75;
        ecc.k = ecc.r * ecc.n;
        inputs.num_msg_bits = ecc.k;
    case 'none'
        inputs.num_msg_bits = inputs.num_msg_enc;
end

disp('========= Turbo Simulation =========');
fprintf('SNR: %f, Noise Var: %f, chan L: %d\n', chann.SNR, chann.nVar, length(chann.h));
fprintf('train num: %d, msg num %d (%f)\n',inputs.num_train_bits,inputs.num_msg_bits, inputs.num_train_bits/inputs.num_msg_bits);

%% init

% create bits
inputs.msg_bits = randi([0 1], 1, inputs.num_msg_bits)';
inputs.train_bits = randi([0 1], 1, inputs.num_train_bits)';

% set ECC
switch ecc.type
    case 'ldpc'
        ecc.H = dvbs2ldpc(ecc.r);
        ecc.encoder = comm.LDPCEncoder('ParityCheckMatrix',ecc.H);
        ecc.decoder = comm.LDPCDecoder('ParityCheckMatrix',ecc.H, ...
                                       'DecisionMethod' , 'Soft decision', ...
                                       'OutputValue', 'Whole codeword');
        inputs.enc_bits = ecc.encoder(inputs.msg_bits);
    case 'none'
        inputs.enc_bits = inputs.msg_bits;
end

% set Intrlv
intrlv.row = 10;
intrlv.col = inputs.num_msg_enc/intrlv.row;
intrlv.step = 3; % hstep is the slope of the diagonal

% create symbs
inputs.msg_symb = EncryptorPath2(inputs.enc_bits, intrlv);
inputs.train_symb = symbMap(inputs.train_bits);
inputs.num_msg_symb = length(inputs.msg_symb);
inputs.num_train_symb = length(inputs.train_symb);

pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4);
chann.channel = comm.AWGNChannel('NoiseMethod',"Signal to noise ratio (SNR)", 'SNR', chann.SNR);

%% channel

% set oversample
chann.symb_oversamp = kron([inputs.train_symb;inputs.msg_symb],ones(chann.overSamp,1));

% simulate channel
chann.out = chann.channel(conv(chann.symb_oversamp,chann.h, "full"));

chann.train_symb = chann.out(1:chann.overSamp*inputs.num_train_symb+(Ld-1)); % train symbs after channel with memory of L-1
chann.msg_symb = chann.out(chann.overSamp*inputs.num_train_symb+1:end);
chann.num_train = length(chann.train_symb);
chann.num_msg = length(chann.msg_symb);

%% Turbo

% set EQ
mu = 0.001; % update step size
maxiter = 10;
EQ_turbo = AdaEQ(Lr, Ld,mu, chann.overSamp); % set
hard = @(x) getHard(x);

% train equlizer
eq.train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
eq.train_err = calcError(eq.train_symb,inputs.train_symb, hard);
eq.train_mse = mean(abs(eq.train_symb-inputs.train_symb));
disp(['Turbo - train BER: ' num2str(eq.train_err) '  MSE:' num2str(eq.train_mse)]);

% set params
msg_pre_in_symb = inputs.train_symb(end-Ld+1:end);
msg_pre_chan_symb = chann.train_symb(end-(Ld-1)-Lr+1:end-(Ld-1));
symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
inds = (1:inputs.num_msg_symb);

% equlize data
figure;hold on;
for i = 1:maxiter
    if i == 1 % have no dn_ yet => run simple eq
        symb_eq = EQ_turbo.normalEqualize_run(symb_chan_input, inputs.num_msg_symb);
    else
        symb_dn_input = [msg_pre_in_symb;dn_]; % add L smples for time channel response 
        symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input, inputs.num_msg_symb);
    end
    [dn_, Decoded] = DecoderPath2(symb_eq, ecc, intrlv, pskDemod); % decode and re-encode
    %print and plot
    eq.err_eq(i) = calcError(Decoded, inputs.msg_bits);
    eq.msg_mse(i) = mean(abs(symb_eq - inputs.msg_symb).^2);
    SE_bit(inds) = abs(symb_eq - inputs.msg_symb).^2;
    plot(inds,SE_bit(inds)); drawnow;
    fprintf(['Turbo - iteration ' num2str(i) ': BER: ' num2str(eq.err_eq(i)) '  MSE:' num2str(eq.msg_mse(i))]); fprintf('\n');
    inds = inds + inputs.num_msg_symb;
end
eq.msg_symb = symb_eq;

Stat = EQ_turbo.Statistics;
xlabel("Filters Update number"); ylabel("bit SE");
title("bit MSE on each filters update"); grid on;
hold off;

figure;semilogy(smooth(SE_bit,1000)); grid on;

%% Turbo Combine

% set EQ
mu = 0.001; % update step size
maxiter = 10;
EQ_turbo = AdaEQ(Lr, Ld,mu, chann.overSamp); % set
hard = @(x) getHard(x);

numTotal = inputs.num_msg_symb + inputs.num_train_symb;
isTrain = logical([ones(inputs.num_train_symb,1);zeros(inputs.num_msg_symb,1)]);
inds = (1:inputs.num_msg_symb);
% equlize data
figure;hold on;
for i = 1:maxiter
    if i == 1 % have no dn_ yet => run simple eq
        [symb_eq, train_eq] = EQ_turbo.normalEqualizeCombined_run(chann.out, inputs.train_symb, isTrain, numTotal);
        symb_eq = circshift(symb_eq,-1);
    else
        [symb_eq, train_eq] = EQ_turbo.turboEqualizeCombined(chann.out,dn_, inputs.train_symb, isTrain, numTotal);
    end
    [dn_, Decoded] = DecoderPath2(symb_eq, ecc, intrlv, pskDemod); % decode and re-encode
    %print and plot
    eq.err_eq(i) = calcError(Decoded, inputs.msg_bits);
    eq.msg_mse(i) = mean(abs(symb_eq - inputs.msg_symb).^2);
    SE_bit(inds) = abs(symb_eq - inputs.msg_symb).^2;
    plot(inds,SE_bit(inds)); drawnow;
    fprintf(['Turbo - iteration ' num2str(i) ': BER: ' num2str(eq.err_eq(i)) '  MSE:' num2str(eq.msg_mse(i))]); fprintf('\n');
    inds = inds + inputs.num_msg_symb;
    train_err = calcError(train_eq,inputs.train_symb, hard);
    train_mse = mean(abs(train_eq-inputs.train_symb));
    disp(['Turbo - train BER: ' num2str(train_err) '  MSE:' num2str(train_mse)]);

end
eq.msg_symb = symb_eq;

Stat = EQ_turbo.Statistics;
xlabel("Filters Update number"); ylabel("bit SE");
title("bit MSE on each filters update"); grid on;
hold off;

figure;semilogy(smooth(SE_bit,1000)); grid on;

%% Matlab DFE

DFE = comm.DecisionFeedbackEqualizer('Algorithm','LMS', ...
    'NumForwardTaps',2*Lr+1,'NumFeedbackTaps',2*Ld+1,'StepSize',mu,...
    'Constellation',[1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2),...
    'ReferenceTap',1,'InputSamplesPerSymbol',2);

[Matlab_Equalized,~,~] = DFE(chann.out ,inputs.train_symb);

Matlab_Equalized = Matlab_Equalized((length(inputs.train_symb)+1):end);

Matlab_error = (abs(inputs.msg_symb - Matlab_Equalized(1:(end-2)))).^2;

MSE_matlab = mean(Matlab_error);

Matlab_BER = calcError(inputs.msg_symb,Matlab_Equalized(1:(end-2)), hard);

figure();
semilogy(Matlab_error,'r');
grid on;
xlabel("bit number");
ylabel("SE");


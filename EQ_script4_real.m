close all;
clear;
clc;

%% settings

load('./LiatModem/11-Mar-2019 09_12_39.mat');
Ld = 2;
Lr = 12;

% load('./LiatModem/11-Mar-2019 09_11_33.mat');


% load('./LiatModem/11-Mar-2019 10_28_43.mat');
% Ld = 2;
% Lr = 8;



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
                                   
inputs.num_msg_enc = ecc.n;
inputs.num_msg_bits = ecc.k;
                                   
% set Intrlv
intrlv.row = record.intrlv.row;
intrlv.col = record.intrlv.col;
intrlv.step = 3; % hstep is the slope of the diagonal

pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4,...
                                 'Variance',1);
L2P = @(x) 1./(1+exp(x));



inputs.train_symb = record.RX.train_symb(:);
inputs.num_train_symb = length(inputs.train_symb);
inputs.msg_bits = record.TX.info_msg_with_CRC(:);
inputs.msg_symb = record.TX.dataSymbls(:);
inputs.num_msg_symb = length(record.TX.dataSymbls);
inputs.padd_bits = record.TX.codeBits_zp(ecc.n+1:end)';

chann.overSamp = 2;

chann.out = record.RX.chann_out(:);
chann.num_out = length(chann.out);

chann.train_symb = chann.out(1:chann.overSamp*inputs.num_train_symb+(Ld-1)); % train symbs after channel with memory of L-1
chann.msg_symb = chann.out(chann.overSamp*inputs.num_train_symb+1:end);
chann.num_train = length(chann.train_symb);
chann.num_msg = length(chann.msg_symb);

Matlab_DFE = comm.DecisionFeedbackEqualizer('Algorithm','RLS', ...
    'NumForwardTaps',12,'NumFeedbackTaps',4,...
    'Constellation',[1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2),...
    'ReferenceTap',1,'InputSamplesPerSymbol',2);


%%

% set EQ
mu = 0.05; % update step size
EQ_turbo = AdaEQ(Lr, Ld,mu, chann.overSamp); % set
hard = @(x) getHard(x);

% train equlizer
eq.train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
eq.train_err = calcError(eq.train_symb,inputs.train_symb, hard);
eq.train_mse = mean(abs(eq.train_symb-inputs.train_symb).^2);
disp(['Turbo - train BER: ' num2str(eq.train_err) '  MSE:' num2str(eq.train_mse)]);

% figure;plot(abs(eq.train_symb-inputs.train_symb).^2);


% set params
maxiter = 10;
EQ_turbo.Mu = 0.05;
msg_pre_in_symb = inputs.train_symb(end-Ld+1:end);
msg_pre_chan_symb = chann.train_symb(end-(Ld-1)-Lr+1:end-(Ld-1));
symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
inds = (1:inputs.num_msg_symb);

% equlize data
figure;hold on;
for i = 1:maxiter
    if i == 1 % have no dn_ yet => run simple eq
        symb_eq = record.RX.res.dfe_out_soft(inputs.num_train_symb+1:end);
%         [eq_matlab,~,~] = Matlab_DFE(chann.out ,inputs.train_symb);
%         symb_eq = eq_matlab(inputs.num_train_symb+1:end);
    else
        symb_dn_input = [msg_pre_in_symb;dn_]; % add L smples for time channel response 
        symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input, inputs.num_msg_symb);
    end
    
    [dn_, Decoded] = DecoderPath_real(symb_eq, ecc, intrlv, pskDemod, inputs.padd_bits);
    
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

figure;semilogy(smooth(SE_bit,1000)); grid on;


%%
% sz = 500;
% ids = 1:sz;
% x = [];
% idx = [];
% for i = 1:floor(inputs.num_msg_symb/sz)
%     x(i) = (symb_eq(ids)' * inputs.msg_symb(ids)) / sz;
%     ids = ids + sz;
% end
% 
% figure;plt3(x); title('xcor 3D');
% t = sz:sz:i*sz;
% 
% figure;plot(t,abs(x));
% figure;plot(t,unwrap(angle(x)));
% 
% ph = exp(-1j*angle(x));
% 
% x = [];
% idx = [];
% ids = 1:sz;
% for i = 1:floor(inputs.num_msg_symb/sz)
%     x(i) = ph(i) * (symb_eq(ids)' * inputs.msg_symb(ids)) / sz;
%     ids = ids + sz;
% end
% 
% figure;plt3(x); title('xcor 3D');
% t = sz:sz:i*sz;
% 
% figure;plot(t,abs(x));
% figure;plot(t,unwrap(angle(x)));
% 
% 
% %% MATLAB
% mu = 0.01; % update step size
% hard = @(x) getHard(x);
% 
% % rls(0.99,1)
% 
% Matlab_DFE = comm.DecisionFeedbackEqualizer('Algorithm','RLS', ...
%     'NumForwardTaps',12,'NumFeedbackTaps',4,'StepSize',mu,...
%     'Constellation',[1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2),...
%     'ReferenceTap',1,'InputSamplesPerSymbol',2);
% 
% % sa = smooth(unwrap(angle(chann.out)),100);
% % e = exp(1j*sa);
% 
% [eq_matlab,~,~] = Matlab_DFE(chann.out ,inputs.train_symb);
% 
% % Matlab_Equalized = Matlab_Equalized((length(inputs.train_symb)+1):end);
% % 
% % Matlab_Equalized = Matlab_Equalized(1:(end-((Ld-1)/2)));
% 
% train_eq_matlab = eq_matlab(1:inputs.num_train_symb);
% msg_eq_matlab = eq_matlab(inputs.num_train_symb+1:end);
% 
% Matlab_error = (abs(inputs.msg_symb - msg_eq_matlab)).^2;
% 
% MSE_matlab = mean(Matlab_error);
% 
% Matlab_BER = calcError(inputs.msg_symb,msg_eq_matlab, hard);
% fprintf(['MATLAB dfe - BER: ' num2str(Matlab_BER) '  MSE:' num2str(MSE_matlab)]); fprintf('\n');
%     
% figure; plot(Matlab_error,'r');
% grid on;
% xlabel("bit number");
% ylabel("SE");
% 


% close all;
% clear;
% clc;

%% settings

load('./LiatModem/11-Mar-2019 09_12_39.mat');
% load('./LiatModem/11-Mar-2019 09_11_33.mat');
% load('./LiatModem/11-Mar-2019 10_28_43.mat');




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


Ld = 5;
Lr = 2 * Ld;

inputs.num_train_symb = 1000;
inputs.train_symb = record.RX.train_symb(:);
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


%%

% set EQ
mu = 0.01; % update step size
maxiter = 3;
EQ_turbo = AdaEQ(Lr, Ld,mu, chann.overSamp); % set
hard = @(x) getHard(x);

% train equlizer
eq.train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
eq.train_err = calcError(eq.train_symb,inputs.train_symb, hard);
eq.train_mse = mean(abs(eq.train_symb-inputs.train_symb));
disp(['Turbo - train BER: ' num2str(eq.train_err) '  MSE:' num2str(eq.train_mse)]);


% set params
EQ_turbo.Mu = 0.01;
msg_pre_in_symb = inputs.train_symb(end-Ld+1:end);
msg_pre_chan_symb = chann.train_symb(end-(Ld-1)-Lr+1:end-(Ld-1));
symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
inds = (1:inputs.num_msg_symb);

% equlize data
figure;hold on;
for i = 1:maxiter
    if i == 1 % have no dn_ yet => run simple eq
        symb_eq = EQ_turbo.normalEqualize_run(symb_chan_input, inputs.num_msg_symb);
%         [dn_, ~] = DecoderPath_real(getHard(chann.msg_symb(1:chann.overSamp:end)), ecc, intrlv, pskDemod, inputs.padd_bits);
%         symb_dn_input = [msg_pre_in_symb;dn_];
%         symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input, inputs.num_msg_symb);
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
xlabel("Filters Update number"); ylabel("bit SE");
title("bit MSE on each filters update"); grid on;
hold off;

figure;semilogy(smooth(SE_bit,1000)); grid on;



%%
sz = 50;
ids = 1:sz;
x = [];
for i = 1:floor(inputs.num_msg_symb/sz)
    x(i) = (symb_eq(ids)' * inputs.msg_symb(ids)) / sz;
    ids = ids + sz;
end

figure;plt3(x); title('xcor 3D');
t = sz:sz:i*sz;
figure;plot(t,abs(x)); title('xcor');
%% comppadd_bits
% 
% % diff in size from tx to rx
% cut_factor = (length(record.RX.res.dfe_out_hard) - length(record.TX.dataSymbls));
% % cut the beggining of rx
% dfe_out = record.RX.res.dfe_out_hard(cut_factor+1:end);
% 
% dmd = pskDemod(dfe_out);
% 
% bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
% bits_dintrlv = bits_dintrlv(1:ecc.n);
% 
% decoded_llr = ecc.decoder(bits_dintrlv);
% decoded_bits = L2P(decoded_llr);
% Decoded = double(decoded_bits(1:ecc.k) > 0.5);
% 
% % res
% 
% demapErr = sum(abs(dmd/2 - record.RX.res.after_demap(:)));
% dintrlvErr = sum(abs(bits_dintrlv/2 - record.RX.res.after_deintlv(:)));
% bitErr = sum(abs(Decoded - record.TX.info_msg_with_CRC(:)));
% 
% fprintf('map error: %f, dinterlive error: %f, bit error: %d\n', demapErr, dintrlvErr, bitErr)
% 
% %% complete loop
% padd_bits = record.TX.codeBits_zp(ecc.n+1:end)';
% dfe_out_soft = record.RX.res.dfe_out_soft(cut_factor+1:end);
% 
% 
% 
% 
% 
% 
% [dn_, Decoded] = DecoderPath_real(dfe_out_soft, ecc, intrlv, pskDemod, padd_bits);
% 
% 
% % %     demodulate
% %     dmd = pskDemod(dfe_out_soft);
% %     
% % %     deinterlive
% %     bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
% %     bits_dintrlv = bits_dintrlv(1:ecc.n);
% %     
% % %     decode
% %     decoded_llr = ecc.decoder(bits_dintrlv);
% %     decoded_bits = L2P(decoded_llr);
% %     Decoded = double(decoded_bits(1:ecc.k) > 0.5);
% %     
% %     % re-encrypt
% %     decoded_bits_pd = [decoded_bits;padd_bits];
% %     bits_intrv = helscanintrlv(decoded_bits_pd,intrlv.row,intrlv.col,intrlv.step);
% %     dn_ = symbMap(bits_intrv);
% 
% 
% 
% 
% 
% 
% 
% 
% loopbitErr = sum(abs(Decoded - record.TX.info_msg_with_CRC(:)));
% loopSymbErr = sum(abs(dn_ - record.TX.dataSymbls(:)));
% 
% fprintf('CHAIN - BER: %d, SER: %d\n', loopbitErr, loopSymbErr)
% 
% figure;plt3(dn_)
% hold on;plt3(dfe_out)
% 
% 
% %% Iter
% 
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


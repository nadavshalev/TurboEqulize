%% cleaning worksapce 
close all;
clear;
clc;

%% settings

% set params
version = "new";                % new or old version of matlab
rng(43);           
EbN0s = 0:10;                 % Eb/N0 values to check
turbo_iterations = 5;           % number of turbo iterations
chann.overSamp = 1;             % number of samples per symbol
ecc.type = 'none';              % encoding ldpc/none
inputs.num_train_bits = 2^11;   % number of training bits
inputs.num_msg_enc = 64800;     % number of bits in the encoded massage
mu = 0.001;                     % step size for equalizers
code_rate = 0.75;               % code rate of the encoding

% different channels to choose from
Porat_and_al = [2-0.4j,1.5+1.8j,1,1.2-1.3j,0.8+1.6j];
Proakis_A = [0.04,-0.05,0.07,-0.21,-0.5,0.72,0.36,0,0.21,0.03,0.07];
Proakis_B = [0.407,0.815,0.407];
Proakis_C = [0.227,0.460,0.688,0.460,0.227];
Proakis_D = [0.1275,0.450,0.750,0.450,0.1275];
Matlab_Channel2 =[0.986,0.845,0.237,0.12345+0.31j,0];
Matlab_Channel3 =[0.623,0.489+0.234j,0.398j, 0.21,0];
My_Channel = [1,1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256];
new_channel = [1+j,0.9-j,0.8+j,0.7-j,0.6+j,0.5-j,0.4+j,0.3-j,0.2+j,0.1-j,0.05+j];
five_tap = [2-0.4j,1.5+1.8j,1,1.2-1.3j,0.8+1.6j];
C1 = [1,1,1,1,1,1,1,1,1,1,1];

% put the name of the channel you want here
chann.h = Proakis_A;

% normlize channel
chann.h = chann.h/norm(chann.h);

% set filters length with realtion to the channels
Ld = length(chann.h);
Lr = chann.overSamp*length(chann.h);

% saving channels name for printing
channel.name =  num2str(chann.h);
channel.name = regexprep(channel.name,'\s+',' ');
channel.name(channel.name=='i')='j';
channel.name = insertAfter(channel.name,' ',', ');
channel.name = ['[' channel.name ']'];

% settings coding parameters
switch ecc.type
    case 'ldpc'
        ecc.n = inputs.num_msg_enc;
        ecc.r = code_rate;
        ecc.k = ecc.r * ecc.n;
        inputs.num_msg_bits = ecc.k;
    case 'none'
        inputs.num_msg_bits = inputs.num_msg_enc;
end

% print general informatio about the simulation
disp('========= Turbo Simulation =========');
fprintf('train num: %d, msg num %d (%f)\n',inputs.num_train_bits,inputs.num_msg_bits, inputs.num_train_bits/inputs.num_msg_bits);

%% init

% create info bits and training bits
inputs.msg_bits = randi([0 1], 1, inputs.num_msg_bits)';
inputs.train_bits = randi([0 1], 1, inputs.num_train_bits)';

% ECC
switch ecc.type
    case 'ldpc'
        % set Encoding and decoding objects
        ecc.H = dvbs2ldpc(ecc.r);
        ecc.encoder = comm.LDPCEncoder('ParityCheckMatrix',ecc.H);
        ecc.decoder = comm.LDPCDecoder('ParityCheckMatrix',ecc.H, ...
                                       'DecisionMethod' , 'Soft decision', ...
                                       'OutputValue', 'Whole codeword');
        % encode info bits                           
        inputs.enc_bits = ecc.encoder(inputs.msg_bits);
    case 'none'
        inputs.enc_bits = inputs.msg_bits;
end

% set Interleaving paramerts
intrlv.row = 10;
intrlv.col = inputs.num_msg_enc/intrlv.row;
intrlv.step = 3; % hstep is the slope of the diagonal

% interleave and map encoded info bits
inputs.msg_symb = EncryptorPath2(inputs.enc_bits, intrlv);

% only map training bits
inputs.train_symb = symbMap(inputs.train_bits);

% save lengths of trainng and info symbols
inputs.num_msg_symb = length(inputs.msg_symb);
inputs.num_train_symb = length(inputs.train_symb);

% creating demapping object
pskDemod = comm.PSKDemodulator(4,'BitOutput',true,...
                                 'DecisionMethod','Approximate log-likelihood ratio',...
                                 'PhaseOffset',pi/4);
% fits to constelation [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2) 
constelation = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2) ;

% create matlab equalizers for first iteration
if(version == "new")
    
    % define new matlab DFE
    DFE = comm.DecisionFeedbackEqualizer('Algorithm','LMS', ...
       'NumForwardTaps',2*Lr+1,'NumFeedbackTaps',2*Ld+1,'StepSize',mu,...
       'Constellation',constelation,...
       'ReferenceTap',1,'InputSamplesPerSymbol',chann.overSamp);
    
    % define new matlab liniear equlaizer
    lin_eq = comm.LinearEqualizer('Algorithm','LMS', ...
    'NumTaps',2*Lr+1,'StepSize',mu,...
    'Constellation',constelation,...
    'ReferenceTap',1,'InputSamplesPerSymbol',chann.overSamp); 
else
    
    % define old matlab DFE
    DFE = dfe(2*Lr+1,2*Ld+1,lms(mu),constelation);
    DFE.nSampPerSym = chann.overSamp;   
    DFE.RefTap = 1;
    
    % define old matlab liniear equlaizer
    lin_eq = lineareq(2*Lr+1,lms(mu),constelation);
    lin_eq.nSampPerSym = chann.overSamp;
    lin_eq.RefTap = 1;
end                             
                             
% index of current Eb/N0                                                         
EbN0_index = 1;

% containers for storing bers of each equalizer
BER = zeros(turbo_iterations,length(EbN0s));    % Turbo
Matlab_BER = zeros(1,length(EbN0s));            % Matlab DFE
Matlab_BER_lin = zeros(1,length(EbN0s));        % Matlab linear
our_lin_BER = zeros(1,length(EbN0s));           % our liniear

%% simulating over different EbN0s
for EbNo = EbN0s
    % print current EbN0
    disp("EbNo = " +num2str(EbNo));
    
    % add AWGN acoording to current EbN0
    chann.channel = comm.AWGNChannel('NoiseMethod',"Signal to noise ratio (Eb/No)",'EbNo',EbNo,'BitsPerSymbol',2,'SamplesPerSymbol',chann.overSamp);

    % set oversample
    chann.symb_oversamp = kron([inputs.train_symb;inputs.msg_symb],ones(chann.overSamp,1));

    % simulate passing through the channel
    chann.out = chann.channel(conv(chann.symb_oversamp,chann.h, "full"));

    % split to train symbols and info symbols from channel output
    chann.train_symb = chann.out(1:chann.overSamp*inputs.num_train_symb+(Ld-1)); % train symbs after channel with memory of L-1
    chann.msg_symb = chann.out(chann.overSamp*inputs.num_train_symb+1:end);
    
    % save lengts of train and info symbols from channel output
    chann.num_train = length(chann.train_symb);
    chann.num_msg = length(chann.msg_symb);

%     %% Matlab Equalizer's
%     %   option to search for best parameters
%     %   LrArr = 2:40;
%     %   LdArr = 2:40;
%     %   muArr = 0.002:0.001:0.02;
%     %   [resData,minData]=searchEQParamsMATLAB(chann, inputs, LrArr, LdArr, muArr);
%     
%     % reset whegits from previous iteration
%     reset(DFE);
%     reset(lin_eq);
%     
%     % train and equalize through matlab equalizers
%     if(version == "new")
%         [DFE_Equalized,~,~] = DFE(chann.out ,inputs.train_symb);
%         [Linear_Equalizer,~,~] = lin_eq(chann.out ,inputs.train_symb);
%     else
%         DFE_Equalized = equalize(DFE,chann.out ,inputs.train_symb);
%         Linear_Equalizer = equalize(lin_eq,chann.out ,inputs.train_symb);
%     end
%     
%     % cutt training equalzied symbol fron the output
%     DFE_Equalized = DFE_Equalized((length(inputs.train_symb)+1):end);
%     
%     Linear_Equalizer = Linear_Equalizer((length(inputs.train_symb)+1):end);
% 
%     % cut extra bits fron the output
%     DFE_Equalized = DFE_Equalized(1:(end-((length(chann.h)-1)/chann.overSamp)));
%     
%     Linear_Equalizer = Linear_Equalizer(1:(end-((length(chann.h)-1)/chann.overSamp)));
%     
%     % passing the equalzied symbols in the decoding and rencoding scheme
%     [dn_DFE, Decoded_DFE] = DecoderPath2(DFE_Equalized, ecc, intrlv, pskDemod); 
% 
%     [dn_Liniear, Decoded_LIN] = DecoderPath2(Linear_Equalizer, ecc, intrlv, pskDemod); 
% 
%     % calculate BERs for Matlab's equalizers
%     Matlab_BER(EbN0_index) = calcError(inputs.msg_bits,Decoded_DFE);
%     
%     Matlab_BER_lin(EbN0_index) = calcError(inputs.msg_bits,Decoded_LIN);
%     

    %% Turbo

    % set EQ
    EQ_turbo = AdaEQ(Lr, Ld,mu, chann.overSamp); % set
    hard = @(x) getHard(x);
    
    % set linear equlazier for first iteration
    our_linear = AdaEQ(Lr, Ld,mu, chann.overSamp);

    % train turbo equlizer
    eq.train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);

    % train liniear equalizer for first iteration
    our_linear.liniearEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
        
    % add symbols from end of training to the start of the info massage for
    % the filters to work with
    msg_pre_chan_symb = chann.train_symb(end-(Ld-1)-Lr+1:end-(Ld-1));
    symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
    msg_pre_in_symb = inputs.train_symb(end-Ld+1:end);

    % just vector for plotting MSE
    inds = (1:inputs.num_msg_symb);
    
    % equalize symbols with the liniear equalizer
    our_liniear_equalized = our_linear.liniearEqualize(symb_chan_input, inputs.num_msg_symb);
    
    % passing the equalzied symbols in the decoding and rencoding scheme
    [dn_our_liniear, Decoded_our_liniear] = DecoderPath2(our_liniear_equalized, ecc, intrlv, pskDemod);
    
    % calculate BER of the liniear equalizer
    our_lin_BER(EbN0_index) = calcError(inputs.msg_bits,Decoded_our_liniear);
    
    % set first dn_ to the liniear equalizer output
    dn_ = dn_our_liniear;

    % equlize data with turbo equalzier
    % plot MSE downgrade on last EbNo iteration
    if (EbNo == EbN0s(end))
        figure;
        hold on;
        grid on;
    end    
    for i = 1:turbo_iterations
        
        symb_dn_input = [msg_pre_in_symb;dn_]; % add L smples for time channel response 
        
        % equalize symbols
        symb_eq = EQ_turbo.turboEqualize(symb_chan_input,symb_dn_input, inputs.num_msg_symb);
        
        % decode and rencode
        [dn_, Decoded] = DecoderPath2(symb_eq, ecc, intrlv, pskDemod); 
        
        % store data
        eq.err_eq(i) = calcError(Decoded, inputs.msg_bits);
        BER(i,EbN0_index) = calcError(Decoded, inputs.msg_bits);
        eq.msg_mse(i) = mean(abs(symb_eq - inputs.msg_symb).^2);
        SE_bit(inds) = abs(symb_eq - inputs.msg_symb).^2;
        
        % plot MSE downgrade on last EbNo iteration
        if (EbNo == EbN0s(end))
            plot(inds,SE_bit(inds)); drawnow;
            fprintf(['Turbo - iteration ' num2str(i) ': BER: ' num2str(eq.err_eq(i)) '  MSE:' num2str(eq.msg_mse(i))]); fprintf('\n');
        end   
        inds = inds + inputs.num_msg_symb;
    end
    eq.msg_symb = symb_eq;

    Stat = EQ_turbo.Statistics;

    xlabel("Filters Update number"); ylabel("bit SE");
    title("bit MSE on each filters update Eb/N0 = "+num2str(EbNo));
    grid on;
    hold off;

    EbN0_index = EbN0_index +1;
end

% printing BER vs. EbNo graphs
figure();

% print Matlab DFE BER
% semilogy(EbN0s,Matlab_BER,'-o');
% hold on;

% print our liniear equalizer BER
semilogy(EbN0s,our_lin_BER,'-x');
hold on;
%semilogy(EbN0s,Matlab_BER_lin,'-x');

grid on;
xlabel("Eb/N0 [dB]");
ylabel("Bit Error Rate");
title("\mu = "+num2str(mu)+" ,            training = " +num2str(inputs.num_train_bits)+" bits,"+ "          Turbo iterations = "+num2str(turbo_iterations)+",           channel = "+channel.name);

% printing turbo equalizer BERS for each turbo iteration
for i = 1:length(BER(:,1))
    % print only first 4 iterations or 10t'h iteration
    if (i<=5) || (i==10)
        semilogy(EbN0s,BER(i,:),'-*');
    end    
end
legend("Linear","Turbo 1","Turbo 2","Turbo 3","Turbo 4","Turbo 5");
hold on;



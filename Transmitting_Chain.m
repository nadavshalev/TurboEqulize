clear all;
clc;
close all;

config = 1;

LDPCWord = 1; % Number of LDPC words in one transmission packet

nCRC = 16;
kLDPC = 48600;
nLDPC = 64800;
R_code = kLDPC/nLDPC;

nInformationBits =  kLDPC-nCRC;

M = 4;

NullFct = 18;
PilotFctVec = 2;
PilotFct = PilotFctVec;

NCarriersVVec = [2048, 2622, 1024, 1536] + 1;

NCarriersV = NCarriersVVec(config);

Nfft_vec = [4096,4096,2048,2048];
Nfft = Nfft_vec(config);
    
NullIndx = [10:NullFct:(NCarriersV-1)/2, (NCarriersV+3)/2+9:18:NCarriersV, NCarriersV];
Nn = length(NullIndx); % Number of null symbols in each block
    
PilotIndx = [1:PilotFct:1+(NCarriersV-1)/2, (NCarriersV+3)/2:PilotFct:(NCarriersV-PilotFct+1)]; % The indexes of the pilot tones used for channel estimation (The indexes are with respect to the signal bandwidth and not the entire spectrum[-fs/2,fs/2] ). % SYSTEM
Np = length(PilotIndx); % Number of pilot symbols in each block
      
DataIndx = (1:NCarriersV);
DataIndx([PilotIndx,NullIndx]) = [];
    
Ns_OFDM_block = NCarriersV - Np - Nn; % number of data symbols in 1 OFDM block
     
Nb_OFDM = ceil(LDPCWord * nLDPC / (Ns_OFDM_block * log2(M))); % Number of blocks in long OFDM packet

nCodeBits = Nb_OFDM * Ns_OFDM_block * log2(M);

padding_size = nCodeBits-nLDPC;

filename = 'C:\Users\Asaf Gendler\OneDrive - Technion\Project_Turbo\Liat\MATLAB_code\files\txOutput_config_1.mat';
load(filename);


massage_length = 65448;

massage_with_crc = info_msg_with_CRC;

% encoding
ParityCheckMatrix = dvbs2ldpc(3/4);
Encoder = comm.LDPCEncoder(ParityCheckMatrix);
coded_massage = Encoder(massage_with_crc')';

% Padding

padding = codeBits_zp((length(codeBits)+1):end);
coded_massage_with_padding = [coded_massage,padding];

% Interleaving
% rows = 10;
rows = Ns_OFDM_block * log2(M);
% cols = ceil(length(coded_massage_with_padding)/rows);
cols = Nb_OFDM;
intrlvHstep = 3; 

interleaved_coded_massage = helscanintrlv(coded_massage_with_padding',rows,cols,intrlvHstep)';

% maaping
codeBits_mat = reshape(interleaved_coded_massage,log2(M),[]);
constellation = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2);

mapper = comm.QPSKModulator('PhaseOffset',pi/4);

codeBits_integer = codeBits_mat(1,:)*2+codeBits_mat(2,:);
massage_in_symbols = mapper(codeBits_integer').';

% [0,0] ----> 1  ----> (1+j)/sqrt(2)
% [0,1] ----> 2  ----> (-1+j)/sqrt(2)
% [1,0] ----> 3  ----> (1-j)/sqrt(2)
% [1,1] ----> 4  ----> (-1-j)/sqrt(2)

bits = interleaved_coded_massage(:);
rshape = reshape(bits,2,length(bits)/2).';
coplexBits = rshape(:,1)+1j*rshape(:,2);
massage_in_symbols = 2*(coplexBits - (1+1j)/2);
massage_in_symbols = massage_in_symbols(:) /sqrt(2);
massage_in_symbols = (conj(massage_in_symbols*exp(1j*pi/4))*exp(-1j*pi/4)).';

a=2;
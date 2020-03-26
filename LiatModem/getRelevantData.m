% must run liat's script first!
% file_name = '11-Mar-2019 09_12_39.mat';
file_name = '11-Mar-2019 09_11_33.mat';
% file_name = '11-Mar-2019 10_28_43.mat';

record.train_symb = params.rxParams.payloadParams.SC_params.training_symbols;
record.TX.info_msg_bits = info_msg_bits;
record.TX.info_msg_with_CRC = info_msg_with_CRC;
record.TX.crc = [1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  0  1];
record.TX.codeBits = codeBits;
record.TX.codeBits_zp = codeBits_zp;
record.TX.tx_bits = CodeBits_tx;
record.TX.dataSymbls = dataSymbols;
record.intrlv.row = interleaverParams.nIntrlvRows;
record.intrlv.col = interleaverParams.nIntrlvCols;

save(['../../code/LiatModem/' file_name],'record');
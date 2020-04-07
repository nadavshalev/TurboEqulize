% must run liat's script first!
% file_name = '11-Mar-2019 09_12_39_c3.mat';
% file_name = '11-Mar-2019 09_11_33_c3.mat';
% file_name = '11-Mar-2019 10_28_43_c3.mat';
% file_name = '11-Mar-2019 15_24_14_c3.mat';
% file_name = '11-Mar-2019 15_25_09_c3.mat';
file_name = '11-Mar-2019 15_22_22_c3.mat';

record.TX.info_msg_bits = info_msg_bits;
record.TX.info_msg_with_CRC = info_msg_with_CRC;
record.TX.crc = [1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  0  1];
record.TX.codeBits = codeBits;
record.TX.codeBits_zp = codeBits_zp;
record.TX.tx_bits = CodeBits_tx;
record.TX.dataSymbls = dataSymbols;

record.RX.train_symb = params.rxParams.payloadParams.SC_params.training_symbols;
record.RX.chann_out = rxOutput_SC.demodOutput.rxSignal_DFE_in;
record.RX.res.dfe_out_soft = rxOutput_SC.demodOutput.DFEOut.estSymbols_const;
record.RX.res.dfe_out_hard = rxOutput_SC.demodOutput.DFEOut.estSymbols;
record.RX.res.res_bits = rxOutput_SC.rx_msg_bits;
record.RX.res.after_demap = rxOutput_SC.codeBits_rx;
record.RX.res.after_deintlv = rxOutput_SC.LLRs;
record.RX.isTrain = rxOutput_SC.demodOutput.DFEOut.training_vec;


record.intrlv.row = interleaverParams.nIntrlvRows;
record.intrlv.col = interleaverParams.nIntrlvCols;


save(['../../code/LiatModem/' file_name],'record');
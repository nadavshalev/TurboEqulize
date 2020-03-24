    
    % General Information
    
    nCRC = 16;
    kLDPC = 48600;
    nLDPC = 64800;
    nCode = nLDPC;
    R_code = kLDPC/nLDPC;
    %nTrainingSymbols = 2000;
    
    % encoding and decoding

    ParityCheckMatrix = dvbs2ldpc(3/4);
    hEnc = comm.LDPCEncoder(ParityCheckMatrix);
    hDec = comm.LDPCDecoder(ParityCheckMatrix);
    
    % Mapping
    M = 4; % Modulation order
    nBitsPerSymbol = log2(M);
    constellation = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2);
    
    codeBits_map2symbols_sc = codeBits_mat_sc(1,:)*2+codeBits_mat_sc(2,:)+1;
    dataSymbols.sc = symbolParams.constellation(codeBits_map2symbols_sc);
        
    % [0,0] ----> 1  ----> (1+j)/sqrt(2)
    % [0,1] ----> 2  ----> (-1+j)/sqrt(2)
    % [1,0] ----> 3  ----> (1-j)/sqrt(2)
    % [1,1] ----> 4  ----> (-1-j)/sqrt(2)

    
    % Interleaving
    nCodeBits = nLDPC;
    nIntrlvRows.sc = 10;
    nIntrlvCols.sc = nCodeBits/nIntrlvRows.sc;
    intrlvHstep = 3; % hstep is the slope of the diagonal, that is, the amount by which
    %the row index increases as the column index increases by one

    CodeBits_tx.sc = helscanintrlv(codeBits,nRows.sc,nCols.sc,h_step);
    
    
    % Demmaping
    Var = var((demodOutput.rx_est_DataSymbols));
    x = real(DecodedSymbols);
    y = imag(DecodedSymbols);
    xy = reshape(sqrt(2) * [y; x] / Var / 1, 1, []);
    
    % deinterleaving
    rx_codeBits_zp = helscandeintrlv(rx_codeBits,nRows,nCols,h_step);



    
    
    
    
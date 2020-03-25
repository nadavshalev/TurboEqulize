function [dn_, Decoded] = DecoderPath2(sn, ldpc, intrlv, demod)
    L2P = @(x) 1./(1+exp(x));
    
    dmd = demod(sn);
    bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
    decoded_llr = ldpc.decoder(bits_dintrlv);
    decoded_bits = L2P(decoded_llr);
    
    Decoded = double(decoded_bits(1:ldpc.k) > 0.5);
    
    % re-encrypt
    dn_ = EncryptorPath2(decoded_bits, intrlv);
end


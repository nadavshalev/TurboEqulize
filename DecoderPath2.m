function [dn_, Decoded] = DecoderPath2(sn, ecc, intrlv, demod)
    L2P = @(x) 1./(1+exp(x));
    
    dmd = demod(sn);
    bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
    
    switch ecc.type
        case 'ldpc'
            decoded_llr = ecc.decoder(bits_dintrlv);
            decoded_bits = L2P(decoded_llr);
            Decoded = double(decoded_bits(1:ecc.k) > 0.5);
        case 'none'
            decoded_bits = -bits_dintrlv;
            Decoded = double(decoded_bits > 0);
    end
    
    % re-encrypt
    dn_ = EncryptorPath2(decoded_bits, intrlv);
end


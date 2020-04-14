function [dn_, Decoded] = DecoderPath_real(sn, ecc, intrlv, demod, pd)
    L2P = @(x) 1./(1+exp(x));
    
%     demodulate
    dmd = demod(sn);
    
%     deinterlive
    bits_dintrlv = helscandeintrlv(dmd,intrlv.row,intrlv.col,intrlv.step);
    bits_dintrlv = bits_dintrlv(1:ecc.n);
    
%     decode
    decoded_llr = ecc.decoder(bits_dintrlv);
    decoded_bits = L2P(decoded_llr);
    Decoded = double(decoded_bits(1:ecc.k) > 0.5);
    
    % re-encrypt
    decoded_bits_pd = [decoded_bits;pd];
    bits_intrv = helscanintrlv(decoded_bits_pd,intrlv.row,intrlv.col,intrlv.step);
    dn_ = symbMap(bits_intrv);
end


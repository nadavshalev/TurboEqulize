function [dn] = EncryptorPath2(codedBits, intrlv)
    % using LDPC and helscanintrlv    
    bits_intrv = helscanintrlv(codedBits,intrlv.row,intrlv.col,intrlv.step);
    dn = symbMap(bits_intrv);
end


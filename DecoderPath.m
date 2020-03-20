function [dn_, Decoded] = DecoderPath(sn, repNum)
    K = length(sn)/repNum;
    % deinterleaved
    idk = interleave(sn,true);
    % demapping
    ick = symbDemap(idk);
    % decode
    iaktmp = reshape(ick, [repNum, K*2]);
    iak = sum(iaktmp,1)/repNum;
    Decoded = iak > 0.5;
    
    % re-encrypt
    dn_ = EncryptorPath(iak, repNum);
end


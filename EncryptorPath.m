function [dn] = EncryptorPath(bits,repNum)
    % coded
    ck = kron(bits, ones(1,repNum));
    % mapping
    symbs = symbMap(ck);
    % interleaved
    dn = interleave(symbs);
end


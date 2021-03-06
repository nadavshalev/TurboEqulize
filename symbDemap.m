function bits = symbDemap(symbs)
        symbs = symbs(:) * sqrt(2);
        symbs = conj(symbs*exp(1j*pi/4))*exp(-1j*pi/4); % for liat's modem (see symbMap.m)
        coplexBits = symbs / 2 + (1+1j)/2;
        realBits = real(coplexBits);
        imagBits = imag(coplexBits);
        bits = zeros(2*length(symbs),1);
        bits(1:2:end) = realBits;
        bits(2:2:end) = imagBits;
end


function symbs = symbMap(bits)
        bits = bits(:);
        rshape = reshape(bits,2,length(bits)/2).';
        coplexBits = rshape(:,1)+1j*rshape(:,2);
        symbs = 2*(coplexBits - (1+1j)/2);
%          = tanh(real())+1j*tanh(imag()));
end


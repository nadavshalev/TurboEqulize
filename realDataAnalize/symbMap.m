function symbs = symbMap(bits)
        bits = bits(:);
        rshape = reshape(bits,2,length(bits)/2).';
        coplexBits = rshape(:,1)+1j*rshape(:,2);
        symbs = 2*(coplexBits - (1+1j)/2) / sqrt(2);
        
        % set constelation acording liat's modem:
        % [0,0] ----> 1  ----> (1+j)/sqrt(2)
        % [0,1] ----> 2  ----> (-1+j)/sqrt(2)
        % [1,0] ----> 3  ----> (1-j)/sqrt(2)
        % [1,1] ----> 4  ----> (-1-j)/sqrt(2)
        symbs = conj(symbs*exp(1j*pi/4))*exp(-1j*pi/4);
end


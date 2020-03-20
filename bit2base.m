function nums = bit2base(bits,m)
    rshape = reshape(bits,log2(m),length(bits)/log2(m))';
    nums = bin2dec(num2str(rshape));
end


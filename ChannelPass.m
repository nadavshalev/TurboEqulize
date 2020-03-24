function chSymbs = ChannelPass(symbs, h, nvar)
    chSymbs = conv(symbs,h, "full");
    chSymbs = chSymbs + sqrt(nvar/2) * (randn(size(chSymbs)) + 1j * randn(size(chSymbs)));
end


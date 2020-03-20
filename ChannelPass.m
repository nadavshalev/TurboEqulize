function chSymbs = ChannelPass(symbs, h, nvar)
    chSymbs = conv(symbs,h, "full");
    chSymbs = chSymbs + nvar * randn(size(chSymbs));
end


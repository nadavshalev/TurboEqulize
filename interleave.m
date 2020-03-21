function resSymbs = interleave(symbs,toDeinterleave)
    if nargin < 2
        toDeinterleave = false;
    end
    fact = 2;
    if toDeinterleave
        tmp = reshape(symbs,[fact,length(symbs)/fact])';
    else
        tmp = reshape(symbs,[length(symbs)/fact,fact])';
    end
    resSymbs = tmp(:);
end


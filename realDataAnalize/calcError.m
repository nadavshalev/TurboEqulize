function ErrNum = calcError(a, b, hard)
    if nargin < 3
        ErrNum = sum(abs(a - b)) / length(a);
    else
        ErrNum = sum(abs(hard(a) ~= hard(b))) / length(a);
    end
end


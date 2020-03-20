function [resVal] = softXorMod2(softVec)
hardVec = mod(round(softVec),2);
resHard = mod(sum(hardVec),2);

distVal = abs(softVec - round(softVec));
if isempty(distVal(distVal>0))
    resVal = resHard;
else
    resSoft = mean(distVal(distVal>0)); 
    resVal = abs(resHard-resSoft); % norm hard output with soft data
end
end

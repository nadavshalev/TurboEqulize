function chann = channEst(chann_msg, prior_msg, h_len)
osamp = 2;
estArr = zeros(osamp,h_len);

for i = 1:osamp
    [x,l] = xcorr(chann_msg(i:osamp:end), prior_msg(i:osamp:end));
    start_ind = find(l==0);

    estArr(i,:) = round(x(start_ind: start_ind+h_len-1)/ (length(x)/2),2);
end

est(1) = estArr(1,1);
for i=2:h_len
    est(i) = estArr(mod((i-1),osamp)+1, floor((i-1)/osamp)+1) - est(i-1);
end

chann = est(:);
end


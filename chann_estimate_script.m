% clear;
h = [2-0.4j 1.5+1.8j 1 1.2-1.3j 0.8+1.6j];
% h = [1 -1];
lh = length(h);

msg_len = 1e6;
osamp = 2;

bits = randi([0 1], 1, msg_len)';
symb = 2*bits-1;
% symb = symbMap(bits);

oversamp = kron(symb,ones(osamp,1));

c_out = conv(oversamp, h,'same');

[x,l] = xcorr(c_out, oversamp);
start_ind= find(l==0);

est_all = round(x(start_ind-1: start_ind+lh-1)/ (length(x)/2),2);

[x1,l1] = xcorr(c_out(1:2:end), oversamp(1:2:end));
start_ind1 = find(l1==0);

est1 = round(x1(start_ind1: start_ind1+lh-1)/ (length(x1)/2),2);

[x2,l2] = xcorr(c_out(2:2:end), oversamp(2:2:end));
start_ind2 = find(l2==0);

est2 = round(x2(start_ind2: start_ind2+lh-1)/ (length(x2)/2),2);

est(1) = est1(1);
for i=2:lh
    if mod(i,2) == 0
        est(i) = est2(i/2)-est(i-1);
    else
        est(i) = est1((i+1)/2)-est(i-1);
    end
end

est(:)
%%

estArr = zeros(osamp,lh);
for i = 1:osamp
    [x,l] = xcorr(c_out(i:osamp:end), oversamp(i:osamp:end));
    start_ind= find(l==0);

    estArr(i,:) = round(x(start_ind: start_ind+lh-1)/ (length(x)/2),2);
end

est(1) = estArr(1,1);
for i=2:lh-1
    est(i) = estArr(mod((i-1),osamp)+1, floor((i-1)/osamp)+1) - est(i-1);
end
est(:)
%%
ch = channEst(c_out, oversamp, 5);

figure;plot(abs(fft(ch,1e4)));


%% REAL DATA

part = 30;
padd = 2^7;
chLen = 20;

inputs.msg_symb = inputs.msg_symb(:);
chann.msg_symb = chann.msg_symb(:);

prior = kron(inputs.msg_symb,ones(2,1));
data = chann.msg_symb;

cut_part = mod(length(prior),part)-1;

if cut_part >= 0
    prior(end-cut_part:end) = [];
    data(end-cut_part:end) = [];
end


inArr = reshape(prior,[length(prior)/part,part]);
chanArr = reshape(data,[length(prior)/part,part]);
ch = zeros(part, chLen);
% chF = zeros(part,padd);
% figure;hold on;
for i = 1:part
    ch(i,:) = channEst(chanArr(:,i), inArr(:,i), chLen);
    ch(i,:) = ch(i,:) / std(ch(i,:));
%     chF(i,:) = abs(fft(ch(i,:),padd));
%     plot(abs(ch(i,:)),':o');
end

% figure;surf(chF);
% ylabel('time');
% xlabel('freq');
% zlabel('Power');

figure;surf(abs(ch));
xlabel('Coefs');
ylabel('Time Slice');
zlabel('Power');
title('Channel Estimation');















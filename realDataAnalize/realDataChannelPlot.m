function realDataChannelPlot(chann, inputs, n_chann, slice)
if nargin < 3
    n_chann = 20;
end
if nargin == 3
    slice = n_chann;
elseif nargin < 4
    slice = 30;
end

% set prior to oversample 2
prior = kron(inputs.msg_symb(:),ones(2,1));
data = chann.msg_symb(:);

% cut data to round multiply of slice size
cut_part = mod(length(prior),slice)-1;
if cut_part >= 0
    prior(end-cut_part:end) = [];
    data(end-cut_part:end) = [];
end

% reshape into slices
inArr = reshape(prior,[length(prior)/slice,slice]);
chanArr = reshape(data,[length(prior)/slice,slice]);

% calc channel for each slice
ch = zeros(slice, n_chann);
for i = 1:slice
    ch(i,:) = channEst(chanArr(:,i), inArr(:,i), n_chann);
    ch(i,:) = ch(i,:) / std(ch(i,:)); % norm energy
end

% plot
figure;surf(abs(ch));
xlabel('Coefs');
ylabel('Time Slice');
zlabel('Power');
title('Channel Estimation');
end


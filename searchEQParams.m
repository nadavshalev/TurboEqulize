function [resData, minData] = searchEQParams(dn_, chann, inputs, muArr, LrArr, LdArr)
N = length(muArr)*length(LrArr)*length(LdArr);
disp(N);

resData = [];
i = 1;
for mu = muArr
    for Lr = LrArr
        for Ld = LdArr
            
            msg_pre_in_symb = inputs.train_symb(end-Ld+1:end);
            msg_pre_chan_symb = chann.train_symb(end-(Ld-1)-Lr+1:end-(Ld-1));
            symb_chan_input = [msg_pre_chan_symb;chann.msg_symb]; % add L smples for time channel response 
            symb_dn_input = [msg_pre_in_symb;dn_];
            
            EQ_turbo = AdaEQ(Lr, Ld,mu, chann.overSamp, 200); % set
            
            train_symb = EQ_turbo.turboEqualize_train(chann.train_symb,inputs.train_symb, inputs.num_train_symb);
            
            symb_eq = EQ_turbo.turboEqualize( symb_chan_input, symb_dn_input, inputs.num_msg_symb);
            msg_mse = mean(abs(symb_eq - inputs.msg_symb).^2);
            
            rs.mu = mu;
            rs.Lr = Lr;
            rs.Ld = Ld;
            rs.mse = msg_mse;
            resData = [resData;rs];
            
            fprintf('(%d,%d) \tMSE: %f  \t(Lr: %d, \tLd: %d, \tmu: %f)\n', i, N, rs.mse, rs.Lr, rs.Ld, rs.mu);
            i = i+1;
        end
    end
end

[~,i] = min([resData.mse]);
minData = resData(i);

fprintf('best MSE: %f (Lr: %d, Ld: %d, mu: %f)\n', minData.mse, minData.Lr, minData.Ld, minData.mu);

figure;scatter3([resData.Lr], [resData.Ld], [resData.mse], 60, [resData.mu])
ax = gca;
xlabel('Lr')
ylabel('Ld')
zlabel('mse')
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'mu';

end


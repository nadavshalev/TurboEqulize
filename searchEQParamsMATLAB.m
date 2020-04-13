  
function [resData, minData] = searchEQParamsMATLAB(chann, inputs, LrArr, LdArr, muArr)
N = length(LrArr)*length(LdArr)*length(muArr);
disp(N);

resData = [];
i = 1;
for mu = muArr
    for Lr = LrArr
        for Ld = LdArr
            Matlab_DFE = comm.DecisionFeedbackEqualizer('Algorithm','LMS', ...
                'NumForwardTaps',Lr,'NumFeedbackTaps',Ld,...
                'Constellation',[1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i] / sqrt(2),...
                'ReferenceTap',1,'InputSamplesPerSymbol',2,...
                'StepSize', mu);
            [eq_matlab,~,~] = Matlab_DFE(chann.out ,inputs.train_symb);
            symb_eq = eq_matlab(inputs.num_train_symb+1:end);
            symb_eq = symb_eq(1:(end-((length(chann.h)-1)/2)));
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
% ax = gca;
% xlabel('Lr')
% ylabel('Ld')
% zlabel('mse')
% cb = colorbar;                                     % create and label the colorbar
% cb.Label.String = 'mu';

end

% function [beta_est,RMSE,var_beta,SNR,p_value] = GeneraModel(X_opt,Y)

%% INPUTS %%

X_opt = [1,-1,-1,-1,1,-1,-1,1;...
         1,-1,1,-1,-1,1,-1,-1;...
         1,-1,1,-1,1,1,1,-1;...
         1,1,1,1,-1,-1,-1,1;...
         1,-1,-1,1,-1,-1,1,-1;...
         1,-1,1,1,1,-1,-1,1;...
         1,1,1,-1,-1,-1,1,-1;...
         1,1,-1,-1,-1,1,-1,1;...
         1,1,-1,1,1,1,-1,-1;...
         1,-1,-1,1,-1,1,1,-1;...
         1,1,1,1,1,1,1,1;...
         1,1,-1,-1,1,-1,1,1];
     
Y = transpose([10.94,15.79,25.96,35.92,22.92,23.54,47.44,19.8,29.48,17.13,43.75,40.86]);


%% OUTPUTS %%



%% PROCESO %%

n = size(X_opt,1);
p = size(X_opt,2);

beta_est = ((inv(transpose(X_opt)*X_opt))*transpose(X_opt))*Y;
RMSE = sqrt((1/(n - p))*transpose(Y - X_opt*beta_est)*(Y - X_opt*beta_est));

var_beta = diag((RMSE^2)*(inv(transpose(X_opt)*X_opt)));
SNR = beta_est./(sqrt(var_beta));
p_value = (1 - tcdf(abs(SNR),n - 1)) + tcdf(- abs(SNR),n - 1);   %REVISAR

% end
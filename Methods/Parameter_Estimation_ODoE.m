function [p, err_per_round] = Parameter_Estimation_ODoE(D, Y, model)

    k = size(D,2);
    n_exp = size(D,1);
    n_rounds = size(Y,2);
    
    % Getting the X matrix by the design D
    
    X = GeneraX_by_D_ODoE(D,model);
    
    % Estimating the parameters by least squares, for each experimental round
      
    p = ((inv(transpose(X)*X))*transpose(X))*Y;
    
    % Calculate of relative error in the estimation
    
    err_per_round = zeros(1,n_rounds);
    
    Y_model = X*p;
    
    for i_r = 1:n_rounds
        
        err_per_round(1,i_r) = (100/n_exp)*sum( sqrt( ( (Y(:,i_r) - Y_model(:,i_r))./Y(:,i_r) ).^2 ) );
        
    end
        
end
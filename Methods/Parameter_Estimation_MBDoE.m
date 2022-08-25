function [p, err_per_round] = Parameter_Estimation_MBDoE(D, Y, model, n_p, Ap, bp, Aeqp, beqp, lbp, ubp, nonlconp)

    k = size(D,2);
    n_exp = size(D,1);
    n_rounds = size(Y,2);
    
    options = optimoptions('ga','Display','None');
    
    p = zeros(n_p, n_rounds);
    err_per_round = zeros(1, n_rounds);
    
    for i_r = 1:n_rounds
        
        % Function objetive definition
        
        fun = @(p)fobj(p, D, Y(:,i_r), model);
        
        % Optimization
        
        p_opt = ga(fun, n_p, Ap, bp, Aeqp, beqp, lbp, ubp, nonlconp, options);
        p(:, i_r) = transpose(p_opt);
        
        % Calculate of relative error in the estimation, for each round
        
        error = 0;
        
        for i_exp = 1:n_exp
            
            error = error + (100/n_exp)*sqrt(((Y(i_exp,i_r) - model(D(i_exp,:), p_opt))/Y(i_exp,i_r))^2);
            
        end
        
        err_per_round(1,i_r) = error;
        
    end
    
    function val = fobj(p, D, y, model_resp)
        
        val = 0;
        
        for i = 1:size(D,1)
        
            val = val + ((y(i) - model_resp(D(i,:),p))^2);
            
        end
        
    end
end
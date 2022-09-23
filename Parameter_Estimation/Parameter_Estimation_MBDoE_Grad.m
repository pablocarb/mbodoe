function [p, err_per_round] = Parameter_Estimation_MBDoE_Grad(D, Y, model, n_p, seeds, Ap, bp, Aeqp, beqp, lbp, ubp, nonlconp)

    n_exp = size(D,1);
    n_rounds = size(Y,2);
    
    options = optimoptions('fmincon','Display','None');
    
    p = zeros(n_p, n_rounds);
    err_per_round = zeros(1, n_rounds);
    
    % Generating the seeds random initial points, parameters
    
    p_0 = zeros(seeds, n_p);
    
    for i_p = 1:n_p
        
        p_0(:,i_p) = lbp(i_p)*ones(seeds,1) + (ubp(i_p) - lbp(i_p))*rand(seeds,1);
        
    end    
    
    for i_r = 1:n_rounds
        
        i_r
        
        % Function objetive definition
        
        fun = @(p)fobj(p, D, Y(:,i_r), model);
        
        % Multilocal optimization
        
        % - Optimization of each initial point case
        
        P_opt = zeros(seeds, n_p);
        F_opt = zeros(seeds,1);
        
        for i_seed = 1:seeds
        
            [p_opt, f_opt] = fmincon(fun, p_0(i_seed, :), Ap, bp, Aeqp, beqp, lbp, ubp, nonlconp, options);
            
            P_opt(i_seed,:) = p_opt;
            F_opt(i_seed) = f_opt;
        
        end
        
        % - Conserve the best
        
        [~,i_best] = min(F_opt);
        p_best = P_opt(i_best,:);
        
        p(:, i_r) = transpose(p_best);
        
        % Calculate of relative error in the estimation, for each round
        
        error = 0;
        y_span = max(Y(:,i_r)) - min(Y(:,i_r));
        
        for i_exp = 1:n_exp
            
            error = error + (100/n_exp)*sqrt(((Y(i_exp,i_r) - model(D(i_exp,:), p_opt))/y_span)^2);
            
        end
        
        err_per_round(1,i_r) = error;
        
    end
    
    function val = fobj(p, D, y, model_resp)
        
        val = 0;
        
        for i = 1:size(D,1)
        
            val = val + ((y(i) - model_resp(D(i,:),p))^2);  % Sum of quadratic errors
            
        end
        
    end
end
function results = Experiment_Simulation2(fun_exp, D, max_error_perc, n_rounds)

    % Add a noise signal of constant amplitude based on the response span
    
    n_exp = size(D,1);
    results = zeros(n_exp, n_rounds);
    
    y_ideal = zeros(n_exp,1);
    for j = 1:n_exp        
            y_ideal(j) = fun_exp(D(j,:));
    end
    
    y_span = max(y_ideal) - min(y_ideal);
    
    
    
    for i = 1:n_rounds
        error_percentage = - max_error_perc + 2*max_error_perc*rand(n_exp,1);            
        results(:,i) = y_ideal + error_percentage/100*y_span;
    end
    
end
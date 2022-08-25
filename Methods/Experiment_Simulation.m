function results = Experiment_Simulation(fun_exp, D, max_error_perc, n_rounds)

    n_exp = size(D,1);
    results = zeros(n_exp, n_rounds);
    

    
    for i = 1:n_rounds
        for j = 1:n_exp
        
            y_ideal = fun_exp(D(j,:));
            error_percentage = - max_error_perc + 2*max_error_perc*rand;
            results(j,i) = y_ideal + ((error_percentage/100)*y_ideal); 
            %PC:
%            results(j,i) = y_ideal + ((error_percentage/100)*100000); 
        end
    end
    
end
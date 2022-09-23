function [D_opt, j_D_opt] = Detmax_Grad(fun_y, p, n_vars, n_exp, sigma, tol, lambda, delta, n_comp, n_init, A, b, Aeq, beq, lb, ub, nonlcon)

% DETMAX ALGORITHM
% Calculates the optimal desiogn of experiments to estimate the parameters 
% of the model function of response, fun_y, taking the local value of the
% parameters, p, the value of the standard deviation of the Gaussian noise
% of the support points, sigma, and a non-degenerate initial design of
% experiments.

%   INPUTS:
%   fun_y: Model function of response
%   p: Local value of parameters
%   sigma: Standard deviation of Gaussian noise.
%   lambda: longitud de la "excursion"
%
%   OUTPUTS:
%   D_opt: Optimal design experiments
%   j_D_opt: D-criteria value

% PARAMETERS

options = optimoptions('fmincon','Display','None');

% PROCESS

% Find a valid good condition design

D_1 = zeros(n_exp,n_vars);

tries = 1;
max_tries = 100;

for i_k = 1:n_vars
    
    D_1(:,i_k) = lb(i_k) + (ub(i_k) - lb(i_k))*rand(n_exp,1);
end

while det(Fisher(fun_y, D_1, p, sigma)) < tol || isnan(det(Fisher(fun_y, D_1, p, sigma)))

    if tries > max_tries
        
        break
        
    end

    for i_k = 1:n_vars

        D_1(:,i_k) = lb(i_k) + (ub(i_k) - lb(i_k))*rand(n_exp,1);
    end
    
    tries = tries + 1;

end

if tries > max_tries   % Not able to find a good initial design
    
    D_opt = [];
    j_D_opt = 0;
    
else
    
    % Initial design (First iteration)

    k = 1;
    D = D_1;
    F = Fisher(fun_y, D, p, sigma);
    F_inv = inv(F);
    
    % Save the D-criteria values for stop criteria
    
    Dcrit_saved = zeros(1,n_comp);
    
    % Solve the optimization problem (adding "lambda" support points, excursion)
    
    for i_save = 1:n_comp

        for i = 1:lambda
            
            % Generating the initial points for the optimization
            
            X_0 = zeros(n_init,n_vars);
            X_opt = zeros(n_init,n_vars);
            F_opt = zeros(1,n_init);
            
            for i_k = 1:n_vars

                X_0(:,i_k) = lb(i_k) + (ub(i_k) - lb(i_k))*rand(n_init,1);
                
            end
            
            for i_0 = 1:n_init     % Multilocal optimization

                fun = @(x)fobj(fun_y ,x, p, F_inv, sigma);                          % Function to optimize definition
                [x_opt_i, f_opt_i] = fmincon(fun, X_0(i_0,:), A ,b, Aeq, beq, lb, ub, nonlcon,options);   % Optimization
                X_opt(i_0,:) = x_opt_i;
                F_opt(i_0) = f_opt_i;
                
            end
            
            % Selecting the best case
            
            [~,i_opt] = min(F_opt);
            x_opt = X_opt(i_opt,:);

            % Inverse actualization
            sens_opt = Sensivity(fun_y, x_opt, p);
            F_inv = F_inv - ((F_inv*transpose(sens_opt)*sens_opt*F_inv)/(sigma^2 + sens_opt*F_inv*transpose(sens_opt)));

            % Adding point
            D = [D; x_opt];

        end
        
        % Search and sustitute the worst support points

        for i = 1:lambda

            n_supp = size(D,1);
            Delta_val = zeros(1, n_supp);

            for j = 1:n_supp

                s_j = Sensivity(fun_y, D(j,:), p);
                Delta_val(j) = (1/(sigma^2))*(s_j*F_inv*transpose(s_j));  % Cost Function to minimize

            end

            % Select the worst point support
            [~, j_min] = min(Delta_val);
            x_worst = D(j_min,:);

            % Inverse actualization
            sens_worst = Sensivity(fun_y, x_worst, p);
            F_inv = F_inv + ((F_inv*transpose(sens_worst)*sens_worst*F_inv)/(sigma^2 - sens_worst*F_inv*transpose(sens_worst)));

            % Remove the point
            D(j_min, :) = [];

        end
        
        % Save the case
        
        Dcrit_saved(i_save) = det(Fisher(fun_y, D, p, sigma));
        
    end

    % WHILE LOOP FOR STOP CONDITION

    while Dcrit_saved(n_comp) > Dcrit_saved(n_comp - 1) && ((Dcrit_saved(n_comp) - Dcrit_saved(1))/Dcrit_saved(n_comp))*100 > delta

        k = k + 1;

        % Solve the optimization problem (adding "lambda" support points, excursion)

        for i = 1:lambda
            
            % Generating the initial points for the optimization
            
            X_0 = zeros(n_init,n_vars);
            X_opt = zeros(n_init,n_vars);
            F_opt = zeros(1,n_init);
            
            for i_k = 1:n_vars

                X_0(:,i_k) = lb(i_k) + (ub(i_k) - lb(i_k))*rand(n_init,1);
                
            end
            
            for i_0 = 1:n_init     % Multilocal optimization

                fun = @(x)fobj(fun_y ,x, p, F_inv, sigma);                          % Function to optimize definition
                [x_opt_i, f_opt_i] = fmincon(fun, X_0(i_0,:), A ,b, Aeq, beq, lb, ub, nonlcon,options);   % Optimization
                X_opt(i_0,:) = x_opt_i;
                F_opt(i_0) = f_opt_i;
                
            end
            
            % Selecting the best case
            
            [~,i_opt] = min(F_opt);
            x_opt = X_opt(i_opt,:);

            % Inverse actualization
            sens_opt = Sensivity(fun_y, x_opt, p);
            F_inv = F_inv - ((F_inv*transpose(sens_opt)*sens_opt*F_inv)/(sigma^2 + sens_opt*F_inv*transpose(sens_opt)));

            % Adding point
            D = [D; x_opt];

        end

        % Search and sustitute the worst support points

        for i = 1:lambda

            n_supp = size(D,1);
            Delta_val = zeros(1, n_supp);

            for j = 1:n_supp

                s_j = Sensivity(fun_y, D(j,:), p);
                Delta_val(j) = (1/(sigma^2))*(s_j*F_inv*transpose(s_j));  % Cost Function to minimize

            end

            % Select the worst point support
            [~, j_min] = min(Delta_val);
            x_worst = D(j_min,:);

            % Inverse actualization
            sens_worst = Sensivity(fun_y, x_worst, p);
            F_inv = F_inv + ((F_inv*transpose(sens_worst)*sens_worst*F_inv)/(sigma^2 - sens_worst*F_inv*transpose(sens_worst)));

            % Remove the point
            D(j_min, :) = [];

        end
        
        % Save the case as the last, remove the first one
        
        Dcrit_saved(1) = [];
        Dcrit_saved = [Dcrit_saved det(Fisher(fun_y, D, p, sigma))];

    end

    % Outputs

    D_opt = D;
    j_D_opt = det(Fisher(fun_y, D_opt, p, sigma));

end

% ADDITIONAL FUNCTIONS

function val = fobj(fun_y ,x, p, F_inv, sig)
    
    sens = Sensivity(fun_y, x, p);
    val = - (1/(sig^2))*(sens*F_inv*transpose(sens));
    
end

end

function [D_opt, j_D_opt] = Detmax_disc_ga_norm(fun_y, p, n_levels, n_exp, sigma, tol, fp_dot_min, fp_dot_max, lambda, A, b, Aeq, beq, lb, ub, nonlcon)

% DETMAX ALGORITHM
% Calculates the optimal design of experiments to estimate the parameters 
% of the model function of response, fun_y, taking the local value of the
% parameters, p, the value of the standard deviation of the Gaussian noise
% of the support points, sigma, and a non-degenerate initial design of
% experiments.

%   INPUTS:
%   fun_y: Model function of response
%   p: Local value of parameters
%   n_levels: Number of levels on each factor (coordinate of support point)
%   n_exp: Number of experiments
%   sigma: Standard deviation of Gaussian noise.
%   tol: Tolerance for the  Fisher determinant of the initial non-degenerate experiment design matrix.
%   lambda: longitud de la "excursion"
%
%   OUTPUTS:
%   D_opt: Optimal design experiments
%   j_D_opt: D-criteria value

% PARAMETERS

n_vars = length(n_levels);  % Dimension of the support points/ Number of factors
intcon = 1:n_vars;

options = optimoptions('ga','Display','None');

% PROCESS

% Find a valid good condition design

D_1 = zeros(n_exp,n_vars);

tries = 1;
max_tries = 100;

for i_k = 1:n_vars
    
    D_1(:,i_k) = lb(i_k) + ((ub(i_k) - lb(i_k))/(n_levels(i_k) - 1))*(randi(n_levels(i_k), n_exp, 1) - 1);

end

while det(Fisher_norm(fun_y, D_1, p, sigma, fp_dot_min, fp_dot_max)) < tol || isnan(det(Fisher_norm(fun_y, D_1, p, sigma, fp_dot_min, fp_dot_max)))

    if tries > max_tries
        
        break
        
    end

    for i_k = 1:n_vars

        D_1(:,i_k) = lb(i_k) + ((ub(i_k) - lb(i_k))/(n_levels(i_k) - 1))*(randi(n_levels(i_k), n_exp, 1) - 1);

    end
    
    tries = tries + 1;

end

if tries > max_tries
    
    D_opt = [];
    j_D_opt = [];
    
else
    % Initial design (First iteration)

    k = 1;
    D = D_1;
    F = Fisher_norm(fun_y, D, p, sigma, fp_dot_min, fp_dot_max);
    F_inv = inv(F);

    % Solve the optimization problem (adding "lambda" support points, excursion)

    for i = 1:lambda

        fun = @(level)fobj(fun_y ,level, n_levels, lb, ub, p, F_inv, sigma, fp_dot_min, fp_dot_max);   % Function to optimize definition
        level_opt = ga(fun, n_vars, A ,b, Aeq, beq, ones(1,n_vars), n_levels, nonlcon, intcon, options);      % Optimization
        x_opt = level2point(level_opt, n_levels, lb, ub);

        % Inverse actualization

        sens_opt = Sensivity_norm(fun_y, x_opt, p, fp_dot_min, fp_dot_max);
        F_inv = F_inv - ((F_inv*transpose(sens_opt)*sens_opt*F_inv)/(sigma^2 + sens_opt*F_inv*transpose(sens_opt)));

        % Adding point
        D = [D; x_opt];

    end

    % Search and sustitute the worst support points

    for i = 1:lambda

        n_supp = size(D,1);
        Delta_val = zeros(1, n_supp);

        for j = 1:n_supp

            s_j = Sensivity_norm(fun_y, D(j,:), p, fp_dot_min, fp_dot_max);
            Delta_val(j) = (1/(sigma^2))*(s_j*F_inv*transpose(s_j));  % Cost Function to minimize

        end

        % Select the worst point support
        [~, j_min] = min(Delta_val);
        x_worst = D(j_min,:);

        % Inverse actualization
        sens_worst = Sensivity_norm(fun_y, x_worst, p, fp_dot_min, fp_dot_max);
        F_inv = F_inv + ((F_inv*transpose(sens_worst)*sens_worst*F_inv)/(sigma^2 - sens_worst*F_inv*transpose(sens_worst)));

        % Remove the point
        D(j_min, :) = [];

    end

    % WHILE LOOP FOR STOP CONDITION

    while det(Fisher_norm(fun_y, D, p, sigma, fp_dot_min, fp_dot_max)) > det(Fisher_norm(fun_y, D_1, p, sigma, fp_dot_min, fp_dot_max))

        k = k + 1;
        D_1 = D;

        % Solve the optimization problem (adding "lambda" support points, excursion)

        for i = 1:lambda

            fun = @(level)fobj(fun_y ,level, n_levels, lb, ub, p, F_inv, sigma, fp_dot_min, fp_dot_max);   % Function to optimize definition
            level_opt = ga(fun, n_vars, A ,b, Aeq, beq, ones(1,n_vars), n_levels, nonlcon, intcon, options);      % Optimization
            x_opt = level2point(level_opt, n_levels, lb, ub);

            % Inverse actualization

            sens_opt = Sensivity_norm(fun_y, x_opt, p, fp_dot_min, fp_dot_max);
            F_inv = F_inv - ((F_inv*transpose(sens_opt)*sens_opt*F_inv)/(sigma^2 + sens_opt*F_inv*transpose(sens_opt)));

            % Adding point
            D = [D; x_opt];

        end

        % Search and sustitute the worst support points

        for i = 1:lambda

            n_supp = size(D,1);
            Delta_val = zeros(1, n_supp);

            for j = 1:n_supp

                s_j = Sensivity_norm(fun_y, D(j,:), p, fp_dot_min, fp_dot_max);
                Delta_val(j) = (1/(sigma^2))*(s_j*F_inv*transpose(s_j));  % Cost Function to minimize

            end

            % Select the worst point support
            [~, j_min] = min(Delta_val);
            x_worst = D(j_min,:);

            % Inverse actualization
            sens_worst = Sensivity_norm(fun_y, x_worst, p, fp_dot_min, fp_dot_max);
            F_inv = F_inv + ((F_inv*transpose(sens_worst)*sens_worst*F_inv)/(sigma^2 - sens_worst*F_inv*transpose(sens_worst)));

            % Remove the point
            D(j_min, :) = [];

        end

    end



    % Outputs

    D_opt = D;
    j_D_opt = det(Fisher_norm(fun_y, D_opt, p, sigma, fp_dot_min, fp_dot_max));
    
end

% ADDITIONAL FUNCTIONS

    function val = fobj(fun_y ,level, n_lev, Lb, Ub, p, F_inv, sig, fp_dot_min, fp_dot_max)

        % Calculate the x value selected the level

        x = level2point(level, n_lev, Lb, Ub);

        % Calculate the objetive function

        sens = Sensivity_norm(fun_y, x, p, fp_dot_min, fp_dot_max);
        val = - (1/(sig^2))*(sens*F_inv*transpose(sens));

    end

    function x = level2point(level, n_lev, Lb, Ub)
        
        x = Lb + ((Ub - Lb)./(n_lev - 1)).*(level - 1);
        
    end
end

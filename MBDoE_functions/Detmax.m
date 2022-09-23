function [D_opt, j_D_opt] = Detmax(fun_y, p, sigma, D_1, lambda, A, b, Aeq, beq, lb, ub, nonlcon)

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
%   D_1: Initial non-degenerate xperiment design matrix.
%   lambda: longitud de la "excursion"
%
%   OUTPUTS:
%   D_opt: Optimal design experiments
%   j_D_opt: D-criteria value

% PARAMETERS

n_vars = size(D_1, 2);  % Dimension of the support points
options = optimoptions('ga','Display','None');

% PROCESS

% Initial design (First iteration)

k = 1;
D = D_1;
F = Fisher(fun_y, D, p, sigma);
F_inv = inv(F);

% Solve the optimization problem (adding "lambda" support points, excursion)

for i = 1:lambda
    
    fun = @(x)fobj(fun_y ,x, p, F_inv, sigma);              % Function to optimize definition
    x_opt = ga(fun, n_vars, A ,b, Aeq, beq, lb, ub, nonlcon,options);   % Optimization
    
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

% WHILE LOOP FOR STOP CONDITION

while det(Fisher(fun_y, D, p, sigma)) > det(Fisher(fun_y, D_1, p, sigma))
    
    k = k + 1;
    D_1 = D;

    % Solve the optimization problem (adding "lambda" support points, excursion)

    for i = 1:lambda

        fun = @(x)fobj(fun_y ,x, p, F_inv, sigma);              % Function to optimize definition
        x_opt = ga(fun, n_vars, A ,b, Aeq, beq, lb, ub, nonlcon,options);   % Optimization

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
    
end



% Outputs

D_opt = D;
j_D_opt = det(Fisher(fun_y, D_opt, p, sigma));

% ADDITIONAL FUNCTIONS

function val = fobj(fun_y ,x, p, F_inv, sig)
    
    sens = Sensivity(fun_y, x, p);
    val = - (1/(sig^2))*(sens*F_inv*transpose(sens));
    
end

end

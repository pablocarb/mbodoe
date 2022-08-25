function [D_opt, j_D_opt] = Fedorov(fun_y, p, sigma, D_1, A, b, Aeq, beq, lb, ub, nonlcon)

% FEDOROV ALGORITHM
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
%
%   OUTPUTS:
%   D_opt: Optimal design experiments
%   j_D: D-criteria value

% PARAMETERS

delta = 1e-6;         % Stop Criteria
n_l = size(D_1, 1);     % Number of support points
n_vars = size(D_1, 2);  % Dimension of the support points
n_p = length(p);      % Number of parameters

% PROCESS

% Initial design (First iteration)

tol = 1e-6;

if det(Fisher(fun_y, D_1, p, sigma)) < tol    % Checking the first design
    
    error("Bad conditioning first desing")
    
end

k = 1
D = D_1
j_D = det(Fisher(fun_y, D, p, sigma))
x_add = zeros(n_l, n_vars);     % Posible substitutes, solutions of the first optimization problem
del_max = zeros(1, n_l);        % Value of cost function at the end of the optimization
F_inv = inv(Fisher(fun_y, D_1, p, sigma));

for i = 1:n_l     % Solve the optimization problem for each support point
    
    fun = @(x)fobj(fun_y, D(i,:), x, p, F_inv, sigma);              % Function to optimize definition
    [x_opt, fval] = ga(fun, n_vars, A ,b, Aeq, beq, lb, ub, nonlcon); % Optimization
    x_add(i,:) = x_opt;
    del_max(i) = - fval;
        
end

[cost,i_max] = max(del_max);  % What is the best case?
cost
x_best = x_add(i_max,:);     % The point to put into design
x_i_delete = D(i_max,:);   % The point to eliminate


% Iteration: Solving the optimization problem recursively

while (-1*fobj(fun_y, x_i_delete, x_best, p, F_inv, sigma)) > delta   % STOP Condition
    
    k = k + 1
    
    % Actualization of the design matrix
    
    D(i_max,:) = x_best;
    D
    j_D = det(Fisher(fun_y, D, p, sigma))
    
    % Actualization of inverse Fisher matrix
    
    M = zeros(n_p, 2);
    M(:,1) = 1i*(transpose(Sensivity(fun_y, x_i_delete, p))/(sqrt(sigma^2)));
    M(:,2) = transpose(Sensivity(fun_y, x_best, p))/(sqrt(sigma^2));
    F_inv = (eye(n_p) - F_inv*M*(eye(2) + transpose(M)*F_inv*M)\transpose(M))*F_inv;
    
    % Preallocating the vectors used in optimization
    
    x_add = zeros(n_l, n_vars);     % Posible substitutes, solutions of the first optimization problem
    del_max = zeros(1, n_l);        % Value of cost function at the end of the optimization

    for i = 1:n_l     % Solve the optimization problem for each support point

        fun = @(x)fobj(fun_y, D(i,:), x, p, F_inv, sigma);              % Function to optimize definition
        [x_opt, fval] = ga(fun, n_vars, A ,b, Aeq, beq, lb, ub, nonlcon); % Optimization
        x_add(i,:) = x_opt;
        del_max(i) = - fval;

    end
    
    [cost,i_max] = max(del_max);  % What is the best case?
    cost
    x_best = x_add(i_max,:);     % The point to put into design
    x_i_delete = D(i_max,:);   % The point to eliminate

end

D_opt = D;
j_D_opt = det(Fisher(fun_y, D_opt, p, sigma));

% ADDITIONAL FUNCTIONS

function val = fobj(fun_y ,x_i, x, p, F_inv, sig)
    
    sens_x = Sensivity(fun_y, x, p);
    sens_xi = Sensivity(fun_y, x_i, p);
    
    d1_x = (1/(sig^2))*sens_x*F_inv*transpose(sens_x);
    d1_xi = (1/(sig^2))*sens_xi*F_inv*transpose(sens_xi);
    d2 = (1/(sig^2))*sens_xi*F_inv*transpose(sens_x); 
    
    val = - d1_x + d1_xi - (d2^2) + d1_x*d1_xi;
    
end


end



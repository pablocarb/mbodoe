clc;clear all;close all;

%% INPUTS %%

%COMMON INPUTS

% - Max number of experiments

n_max = 20;

% - Number of initial valid design seeds

seeds = 15;

% - Design Space

Lb = 0;
Ub = 20;

% - Number of factors/dimension of support point

k = 1;

%ODoE INPUTS

% - Min number of experiments

n_min_ODoE = (2*k + 1) + k*(k - 1)/2 + 1;

% - Levels considered per factor

levels_ODoE = 200;

% - Max iterations in CE algorithm

Max_iter = 100;

% - Kind of LP-structure model response

model_ODoE = 'quadratic';

%MBDoE INPUTS

% - Number of parameters in model and minimum number of experiments
% required

%n_p = 4;
n_p = 3; %PC
n_min_MBDoE = n_p + 1;

% - Number of levels to consider

levels_MBDoE = 200;

% - Function with the model response, depends on experimental conditions
% and parameters value

model_MBDoE = @(x,p)Ex_model(x,p);

% - Standard deviation of experimental Gaussian noise (Fisher weight)

sigma = 1;

% - Tolerance of bad condition design

tol = 1e-15;

% - Parameters space and number of initial values

np_init_MBDoE = 15;

Lb_p_MBDoE = [0 0 0 0];
Ub_p_MBDoE = [1 1 1 1];

% - Excursions

lambda = 6;

% - Unsconstrained problem considered, only design space

A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];

% EXPERIMENT SIMULATION INPUTS

% - Experiment parameters

%p_exp = [1 1 1 0.1];
p_exp = [1e8 1e-6 1e-2];

% - Experiment model

%model_experiment = @(x)Ex_model(x,p_exp);
model_experiment = @(x)Model_case_1(x,p_exp); %PC

% - Experiment rounds

n_rounds = 100;

% - Max.Relative Error added (+/- x%)

max_rel_error = 10;

% MODEL VALIDATION

% - Number of points to validate

n_exp_val = 100;

%% PROCESS %%

% - Initializing parallel computing

%parpool

% - Preallocating the final plots results

max_error_ODoE = zeros();
max_error_MBDoE = zeros();

% - Cell structure for results

results_ODoE = cell((n_max - n_min_ODoE) + 1, 7);
results_MBDoE = cell((n_max - n_min_MBDoE) + 1, 7);

%Getting the optimal design of experiments with ODoE 

disp("Optimal design of experiments with ODoE")

tic

index = 1;

for n = n_min_ODoE:n_max
    
    n
    [D_opt,~,Deff_max,~,~,~,~] = Dopt_CE(n,k,levels_ODoE,model_ODoE,Lb,Ub,Max_iter,seeds);
    results_ODoE(index,1) = {n};
    results_ODoE(index,2) = {D_opt};
    results_ODoE(index,3) = {Deff_max};
    index = index + 1;
    
end

t_Design_ODoE = toc;

%Getting the optimal design of experiments with MBDoE

%% PC: MBDoE Section

disp("Optimal design of experiments with MBDoE")

tic

index = 1;

%Getting the maximum/minimum value of the sensitive vector components (to scale)

% disp("Calculating the sensitivy vector limits")
% 
% [fp_dot_min, fp_dot_max] = Sensivity_max_min(model_MBDoE, k, n_p, A, b, Aeq, beq, Lb, Ub, [], [], [], [], Lb_p_MBDoE, Ub_p_MBDoE);

for n = n_min_MBDoE:n_max
    
    n
    % Generate the random initial parameter values
    
    p_init = zeros(np_init_MBDoE,n_p);
    
    for i = 1:n_p
        
        p_init(:,i) = Lb_p_MBDoE(i)*ones(np_init_MBDoE,1) + (Ub_p_MBDoE(i) - Lb_p_MBDoE(i))*rand(np_init_MBDoE,1);
        
    end
    
    % Generate a p_structure to save the experiments designs
    
    p_optimal_designs = cell(np_init_MBDoE, seeds);
    j_d_optimal_p = zeros(np_init_MBDoE, seeds);
    
    for j = 1:np_init_MBDoE
        j
        for s = 1:seeds
            
            % Find the optimal experiments design

            [D_opt_MBDoE,j_D_opt] = Detmax_disc_ga(model_MBDoE, p_init(j,:), levels_MBDoE, n, sigma, tol, lambda, A, b, Aeq, beq, Lb, Ub, nonlcon);

            p_optimal_designs(j,s) = {D_opt_MBDoE};
            j_d_optimal_p(j,s) = j_D_opt;

        end
        
    end
    
    % Save the best design of experiments
    
    [best_per_col,i_rows] = max(j_d_optimal_p,[],1);
    [best_j_d, i_col] = max(best_per_col);
    i_row = i_rows(i_col);
    
    results_MBDoE(index,1) = {n};
    results_MBDoE(index,2) = {p_optimal_designs{i_row, i_col}};
    results_MBDoE(index,3) = {((best_j_d)^(1/n_p))/n};
    
    index = index + 1;
    
end

t_Design_MBDoE = toc;

%Simulating the experiments
%% PC: Simulating the experiments for ODoE designs

disp("Simulating the experiments for ODoE designs")

tic

n_cases_ODoE = n_max - n_min_ODoE + 1;
n_cases_MBDoE = n_max - n_min_MBDoE + 1;

for index = 1:n_cases_ODoE
    
    index
    
    results = Experiment_Simulation(model_experiment, results_ODoE{index,2}, max_rel_error, n_rounds);
    results_ODoE(index,4) = {results};
    
end

t_experiments_ODoE = toc;

%% PC: Simulating the experiments for MBDoE designs

disp("Simulating the experiments for MBDoE designs")

tic

for index = 1:n_cases_MBDoE
    
    index
    
    results = Experiment_Simulation(model_experiment, results_MBDoE{index,2}, max_rel_error, n_rounds);
    results_MBDoE(index,4) = {results};
    
end

t_experiments_MBDoE = toc;

%Parameters estimation for ODoE: Least squares
%% PC: Estimating the parameters for ODoE Models

disp("Estimating the parameters for ODoE Models")


tic

for index = 1:n_cases_ODoE
    
    index
    
    D = results_ODoE{index, 2};
    Y = results_ODoE{index, 4};
    
    % Estimate the parameters of the model, one by each experimental round
    
    [p_est_ODoE, err_per_round_ODoE] = Parameter_Estimation_ODoE(D, Y, model_ODoE);
    
    results_ODoE{index, 5} = p_est_ODoE;
    results_ODoE{index, 6} = err_per_round_ODoE;
    
end

t_model_ODoE = toc;

%Parameters estimation for MBDoE: Numerical optimization process
%% PC: Estimating the parameters for MBDoE Models

disp("Estimating the parameters for MBDoE Models")

tic

for index = 1:n_cases_MBDoE
    
    index
    
    D = results_MBDoE{index, 2};
    Y = results_MBDoE{index, 4};
    
    % Estimate the parameters of the model, one by each experimental round
    
    [p_est_MBDoE, err_per_round_MBDoE] = Parameter_Estimation_MBDoE(D, Y, model_MBDoE, n_p, [], [], [], [], Lb_p_MBDoE, Ub_p_MBDoE, []);
    
    results_MBDoE{index, 5} = p_est_MBDoE;
    results_MBDoE{index, 6} = err_per_round_MBDoE;
    
end

t_model_MBDoE = toc;


%% PC: Error estimation

%Error estimation: Validation of models

% - Generating the random points to validate

Points_val = zeros(n_exp_val,k);

for i = 1:k
    
    Points_val(:,i) = Lb(i) + (Ub(i) - Lb(i))*rand(n_exp_val,1);
    
end
% - Experiments

Y_val = Experiment_Simulation(model_experiment, Points_val, max_rel_error, 1);

% - Calculate the error in these points for each model
%% PC: ODoE validation

disp("Validating the ODoE models")

X = GeneraX_by_D_ODoE(Points_val,model_ODoE);

for index = 1:n_cases_ODoE
    
    parameters_val_ODoE = results_ODoE{index,5};
    n_rounds = size(parameters_val_ODoE,2);
    errors_val = zeros(1,n_rounds);
    n = results_ODoE{index,1};
    
    % Model responses
    
    Y_model = X*parameters_val_ODoE;
    
    % Error
    
    for i_r = 1:n_rounds
        
        errors_val(i_r) = (100/(n_exp_val))* sum(abs((Y_val - Y_model(:,i_r))./Y_val));
        
    end
    
    results_ODoE{index, 7} = errors_val;
    
end

%% PC: MBDoE validation

disp("Validating the MBDoE models")

for index = 1:n_cases_MBDoE
    
    parameters_val_MBDoE = results_MBDoE{index,5};
    n_rounds = size(parameters_val_ODoE,2);
    errors_val = zeros(1,n_rounds);
    n = results_ODoE{index,1};
    
    % Model responses & error
   
    for i_r = 1:n_rounds
        
        Y_model = zeros(n_exp_val,1);
        
        for i_exp = 1:n_exp_val
            
            Y_model(i_exp) = model_MBDoE(Points_val(i_exp), transpose(parameters_val_MBDoE(:,i_r)));
            
        end
        
        errors_val(i_r) = (100/n_exp_val)* sum(abs((Y_val - Y_model)./Y_val));
        
    end
    
    results_MBDoE{index, 7} = errors_val;
    
end

% delete(gcp)

%% PLOTS %%

% Error estimation/validation

n_experimentos_ODoE = n_min_ODoE:n_max;
n_experimentos_MBDoE = n_min_MBDoE:n_max;

err_est_mean_ODoE = zeros(1, n_max - n_min_ODoE);
err_est_max_ODoE = zeros(1, n_max - n_min_ODoE);
err_val_mean_ODoE = zeros(1, n_max - n_min_ODoE);
err_val_max_ODoE = zeros(1, n_max - n_min_ODoE);

err_est_mean_MBDoE = zeros(1, n_max - n_min_ODoE);
err_est_max_MBDoE = zeros(1, n_max - n_min_ODoE);
err_val_mean_MBDoE = zeros(1, n_max - n_min_ODoE);
err_val_max_MBDoE = zeros(1, n_max - n_min_ODoE);

for index = 1:n_cases_ODoE
    
    resu_error_est = results_ODoE{index,6};
    resu_error_val = results_ODoE{index,7};
    
    err_est_mean_ODoE(index) = mean(resu_error_est);
    err_est_max_ODoE(index) = max(resu_error_est);
    err_val_mean_ODoE(index) = mean(resu_error_val);
    err_val_max_ODoE(index) = max(resu_error_val);
    
end

for index = 1:n_cases_MBDoE
    
    resu_error_est = results_MBDoE{index,6};
    resu_error_val = results_MBDoE{index,7};
    
    err_est_mean_MBDoE(index) = mean(resu_error_est);
    err_est_max_MBDoE(index) = max(resu_error_est);
    err_val_mean_MBDoE(index) = mean(resu_error_val);
    err_val_max_MBDoE(index) = max(resu_error_val);
    
end

figure
plot(n_experimentos_ODoE,err_est_mean_ODoE,n_experimentos_MBDoE,err_est_mean_MBDoE)

figure
plot(n_experimentos_ODoE,err_est_max_ODoE,n_experimentos_MBDoE,err_est_max_MBDoE)

% INTERESTING PLOTS

figure
title('Mean Relative Error Experiments vs Models')
plot(n_experimentos_ODoE,err_val_mean_ODoE,'ro',n_experimentos_MBDoE,err_val_mean_MBDoE,'b*')
ylim([0 100])
xlabel("Number of experiments")
ylabel("Mean Rel. Error(%)")
legend("ODoE", "MBDoE")

figure
title("Max Relative Error Experiments vs Models")
plot(n_experimentos_ODoE,err_val_max_ODoE,'ro',n_experimentos_MBDoE,err_val_max_MBDoE,'b*')
ylim([0 100])
xlabel("Number of experiments")
ylabel("Max. Rel. Error(%)")
legend("ODoE", "MBDoE")


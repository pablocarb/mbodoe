function F = Fisher_norm(fun_y, D, p, sigma, fp_dot_min, fp_dot_max)

% Calculates the Fisher matrix of the response model, fun_y,
% experiment design, D, and local value of parameters, p. Assume the same
% Gaussian Noise in all dimensions.
%
%   INPUTS:
%   fun_y: Model function of response
%   D: Experiment design matrix, support points in rows
%   p: Local value of parameters
%   sigma: Standard deviation of Gaussian noise.
%   fp_dot_min: Global minimun of all the parameter derivatives of model
%   response
%   fp_dot_max: Global maximun of all the parameter derivatives of model
%   response
%
%   OUTPUTS:
%   F: Fisher Matrix

% PARAMETERS

n_l = size(D, 1);   % Number of support points
n_p = length(p);    % Número de parámetro

% PROCESS

F = zeros(n_p);

for i = 1:n_l
    
    F = F + (transpose(Sensivity_norm(fun_y, D(i, :), p, fp_dot_min, fp_dot_max))*Sensivity_norm(fun_y, D(i, :), p, fp_dot_min, fp_dot_max));

end

F = (1/(sigma^2))*F;

end
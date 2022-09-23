function F = Fisher(fun_y, D, p, sigma)

% Calculates the Fisher matrix of the response model, fun_y,
% experiment design, D, and local value of parameters, p. Assume the same
% Gaussian Noise in all dimensions.
%
%   INPUTS:
%   fun_y: Model function of response
%   D: Experiment design matrix, support points in rows
%   p: Local value of parameters
%   sigma: Standard deviation of Gaussian noise.
%
%   OUTPUTS:
%   F: Fisher Matrix

% PARAMETERS

n_l = size(D, 1);   % Number of support points
n_p = length(p);    % Número de parámetro

% PROCESS

F = zeros(n_p);

for i = 1:n_l
    
    F = F + (transpose(Sensivity(fun_y, D(i, :), p))*Sensivity(fun_y, D(i, :), p));

end

F = (1/(sigma^2))*F;

end
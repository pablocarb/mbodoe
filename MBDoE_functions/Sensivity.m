function s = Sensivity(fun_y, x, p)

% Calculates the sensivity vector of the response model, fun_y,
% with the parameters at the support point x, and local value of
% parameters, p.
%
%   INPUTS:
%   fun_y: Model function of response
%   x: Support point
%   p: Local value of parameters
%
%   OUTPUTS:
%   s: Sensivity vector

% PARAMETROS

n_p = length(p);     % Número de parámetros

% CÁLCULO

h_p = 1e-6*p;      % Paso relativo al orden de los parámetros

y_h = zeros(1, n_p);
y_mh = zeros(1, n_p);

for i = 1:n_p
    
    p_h = p;  
    p_mh = p; 
    
    p_h(i) = p(i) + h_p(i);
    p_mh(i) = p(i) - h_p(i);

    y_h(i) = fun_y(x, p_h);
    y_mh(i) = fun_y(x, p_mh);

end

s = (y_h - y_mh)./(2*h_p);

end


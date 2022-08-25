function y = Exercise6_2(x, p)

%EXERCISE6_2 of Walter & Pronzato
%   Example of function of response.
%   INPUTS:
%   x(dim = 1): Support point
%   p(dim = 2): Parameters
%   OUTPUTS:
%   y: Model response

y = (p(1)*exp(-p(2)*x));

end



function y = Exercise6_3(x, p)

%EXERCISE6_3 of Walter & Pronzato
%   Example of function of response.
%   INPUTS:
%   x(dim = 1): Support point
%   p(dim = 4): Parameters
%   OUTPUTS:
%   y: Model response

y = (p(1)*exp(-p(2)*x)) + p(3)*exp(-p(4)*x);

end


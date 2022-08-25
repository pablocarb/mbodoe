function y = Model_case_1a(x,p)

% x: Observation time
% p: Parameters--> 1: K_p*K_I
%                  2: d_m
%                  3: d_p

%y = (p(1)/(p(2)*p(3)))*(1 - exp(-p(2)*x) - exp(-p(3)*x));
y = p(1) + p(2)*x + p(3)*x^2;

end
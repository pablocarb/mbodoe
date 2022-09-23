function y = Model_case_1(x,p)

% x: Observation time
% p: Parameters--> 1: K_p*K_I
%                  2: d_m
%                  3: d_p

y = (p(1)/(p(2)*p(3)))*(1 - (1/(1-p(2)/p(3)))*exp(-p(2)*x) - (1/(1-p(3)/p(2)))*exp(-p(3)*x)) ;

end
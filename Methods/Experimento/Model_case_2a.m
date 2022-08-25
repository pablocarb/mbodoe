function y = Model_case_2a(x,p)

%% INPUTS %%
% x: Observation time
% p: Parameters--> 1: cn*Kp
%                  2: KM

%% OUTPUTS %%
% y: Value of product concentration at observation time

%% PROCESS %%

% Defined Parameters

Kcat = 1;
Kd = 1e-4;
K1 = 10;
K_1 = 1e-2;
n = 2;
I = 1e-2;
Kleak = 1e-12;
Kdeg = 1e-6;
PrT = 1;

% Assignation of parameters

cn_Kp = p(1);
KM = p(2);

% Integration

y0 = [1000; 0; 0; 0];
tf = x;

[~,y_values] = ode15s(@(t,y)odefun(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT),[0, tf],y0);
y_final = y_values(end,:);
y = y_final(4);

function dydt = odefun(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT)

    dydt = zeros(4,1);
    
    for i = 1:4
        
        if y(i)<0
            
            y(i) = 0;
            
        end
    end
    
    dydt(1) = - Kd*y(1) - ((Kcat*y(3)*y(1))/(KM + y(1)));  % Sustrato
    dydt(2) = K1*(I^n) - (K_1*y(2));   % Promotor Activado
    dydt(3) = (Kleak + cn_Kp)*PrT - (Kleak + cn_Kp)*y(2) - (Kdeg*y(3));  % Enzima Total
    dydt(4) = ((Kcat*y(3)*y(1))/(KM + y(1)));    % Producto

end

end
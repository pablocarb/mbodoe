clc; clear all; close all;

% function [t,P] = Experimento(Kcat,KM,Kp,Cn)

Kcat = [10,20,100];
KM = [1e3,500,200];
Kp = [1e-9,1e-8,1.5e-8];
Cn = [5,10,100];

% Calcula la producción en una simulación de experimento

%% INPUTS %%

%% OUTPUTS %%

%% PROCESO %%

% Parámetros

Kd = 1e-4;
K1 = 10;
K_1 = 1e-2;
n = 1.85;
I = 1e-2;
Kleak = 1e-12;
Kdeg = 1e-6;
PrT = 10;

% Integración + Detección de evento

y0 = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0];
tf = 100000000;
options = odeset('Events', @(t,y)FinalReact_Event(t,y,Kcat,KM,Kp,Cn,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT),'NonNegative',1:10);

[t,y,te,ye,ie] = ode15s(@(t,y)odefun(t,y,Kcat,KM,Kp,Cn,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT),[0, tf],y0,options);
test = odefun(te,ye,Kcat,KM,Kp,Cn,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT);


function [value,isterminal,direction] = FinalReact_Event(t,y,Kcat,KM,Kp,Cn,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT)

    epsilon = 1e-5;
    Dy_Dt = odefun(t,y,Kcat,KM,Kp,Cn,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT);
    n_Dy_Dt = norm(Dy_Dt,2);
    
    if n_Dy_Dt < epsilon
        
        value = 0;
        
    else
        
        value = 1;
        
    end
    
    isterminal = 1;
    direction = 0;

end

function dydt = odefun(t,y,Kcat,KM,Kp,Cn,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT)

    n_reac = length(Kcat);
    n_ecs = 4 + 3*(n_reac - 1);
    dydt = zeros(n_ecs,1);
    
    for i = 1:n_ecs
        
        if y(i)<0
            
            y(i) = 0;
            
        end
    end
    
    dydt(1) = - Kd*y(1) - ((Kcat(1)*y(3)*y(1))/(KM(1) + y(1)));  % Sustrato de la primera reacción
    
    for i = 1:n_reac
        
        dydt(2 + (i - 1)*3) = K1*(I^n) - (K_1*y(2 + (i - 1)*3)); % Promotor activado 
        dydt(3 + (i - 1)*3) = (Kleak + Cn(i)*Kp(i))*PrT - (Kleak + Cn(i)*Kp(i))*y((2 + (i - 1)*3)) - (Kdeg*y(3 + (i - 1)*3)); % Encima Total
        
        if i == n_reac
            
            dydt(4 + (i - 1)*3) = ((Kcat(i)*y(3 + (i - 1)*3)*y(1 + (i - 1)*3))/(KM(i) + y(1 + (i - 1)*3)));    % Producto final
            
        else
            
            dydt(4 + (i - 1)*3) = ((Kcat(i)*y(3 + (i - 1)*3)*y(1 + (i - 1)*3))/(KM(i) + y(1 + (i - 1)*3))) -...
                - Kd*y(4 + (i - 1)*3) - ((Kcat(i + 1)*y(3 + i*3)*y(4 + (i - 1)*3))/(KM(i + 1) + y(4 + (i - 1)*3)));    % Productos - Sustratos
            
        end
        
    end

end



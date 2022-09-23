%function y = Model_case_2b(x,p)

clc; clear all; close all;

% x

x = [2 -4 6 -4 7]; % Values of amplitude steps

% p

KM = 1e3;
cn_Kp = 1e-8;

%% INPUTS %%
% x: Signal Design, vector of input values in perturbation equal intervals
%    of times
% p: Parameters--> 1: cn*Kp
%                  2: KM
%% OUTPUTS %%
% y: Integral of product concentration in certain time, starting in the
%    initial input signal instant 

%% PARAMETERS %%

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

% Mass Flow and Extraction Rates

S_in = 1;
F_Pr = 0.5;
F_ET = 0.5;
F_P = 0.5;

% Input Signal

Delta_time = 5000;   % Duration of input signal time

%% FIRST STEP: Integration until stationary state

% Integration + Event Detection

y0 = [0; 0; 0; 0];
tf = 2000000;
options = odeset('Events', @(t,y)FinalReact_Event(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P));

[t,y,te,ye,ie] = ode15s(@(t,y)odefun_1(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P),[0, tf],y0, options);

%% SECOND STEP: Introduction of signal

% Definition of time limits of signal

t_init_sign = te;
t_end_sign = te + Delta_time;

% Integration

y0 = transpose(ye);

[t_signal,y_signal] = ode15s(@(t,y)odefun_2(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P,x,t_init_sign,t_end_sign),[t_init_sign, t_end_sign],y0);

% Adding Results

t = [t;t_signal(2:end)];
y = [y;y_signal((2:end),:)];

%% THIRD STEP: Time after signal introduction

Delta_tf = 24*60*60;  % One day after the input signal finish

% Integration

y0 = transpose(y(end,:));
t0 = t(end);
tf = t0 + Delta_tf;

[t_oneday,y_oneday] = ode15s(@(t,y)odefun_1(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P),[t0, tf],y0);

% Adding Results

t = [t;t_oneday(2:end)];
y = [y;y_oneday((2:end),:)];

%% PLOT %%

figure 

% Sustrato

subplot(2,2,1)
plot(t,y(:,1),'b', 'LineWidth', 2)
title('Sustrato')
xlabel('t(s)','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('$S(\mu Mol)$','Interpreter','latex','FontWeight','bold','FontSize',14)
grid on

% Promotor Activado

subplot(2,2,2)
plot(t,y(:,2),'r', 'LineWidth', 2)
title('Prom. Activo')
xlabel('t(s)','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('$P_{r}^{*}(\mu Mol)$','Interpreter','latex','FontWeight','bold','FontSize',14)
grid on

% Enzima

subplot(2,2,3)
plot(t,y(:,3),'g', 'LineWidth', 2)
title('Enzimas')
xlabel('t(s)','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('$E_{T}(\mu Mol)$','Interpreter','latex','FontWeight','bold','FontSize',14)
grid on

% Producto

subplot(2,2,4)
plot(t,y(:,4),'m', 'LineWidth', 2)
title('Producto')
xlabel('t(s)','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('$P(\mu Mol)$','Interpreter','latex','FontWeight','bold','FontSize',14)
grid on

function [value,isterminal,direction] = FinalReact_Event(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P)

    epsilon = 1e-7;
    Dy_Dt = odefun_1(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P);
    n_Dy_Dt = norm(Dy_Dt,2);
    
    if n_Dy_Dt < epsilon
        
        value = 0;
        
    else
        
        value = 1;
        
    end
    
    isterminal = 1;
    direction = 0;

end


function dydt = odefun_1(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P)

    dydt = zeros(4,1);
    
    for i = 1:4
        
        if y(i)<0
            
            y(i) = 0;
            
        end
    end
    
    dydt(1) = - Kd*y(1) - ((Kcat*y(3)*y(1))/(KM + y(1))) + S_in;  % Sustrato
    dydt(2) = K1*(I^n) - (K_1*y(2)) - F_Pr*y(2);   % Promotor Activado
    dydt(3) = (Kleak + cn_Kp)*PrT - (Kleak + cn_Kp)*y(2) - (Kdeg*y(3)) - F_ET*y(3);  % Encima Total
    dydt(4) = ((Kcat*y(3)*y(1))/(KM + y(1))) - F_P*y(4);     %Producto

end

function dydt = odefun_2(t,y,cn_Kp,KM,Kcat,Kd,K1,K_1,n,I,Kleak,Kdeg,PrT,S_in,F_Pr,F_ET,F_P,Amp_u,t_i_sign,t_f_sign)

    dydt = zeros(4,1);
    
    for i = 1:4
        
        if y(i)<0
            
            y(i) = 0;
            
        end
    end
    
    dydt(1) = - Kd*y(1) - ((Kcat*y(3)*y(1))/(KM + y(1))) + S_in + Signal(t, Amp_u, t_i_sign, t_f_sign);  % Sustrato
    dydt(2) = K1*(I^n) - (K_1*y(2)) - F_Pr*y(2);   % Promotor Activado
    dydt(3) = (Kleak + cn_Kp)*PrT - (Kleak + cn_Kp)*y(2) - (Kdeg*y(3)) - F_ET*y(3);  % Encima Total
    dydt(4) = ((Kcat*y(3)*y(1))/(KM + y(1))) - F_P*y(4);     %Producto

end

function signal_value = Signal(t, amp, t_init, t_end)
    
    n_steps = length(amp);
    step_time = (t_end - t_init)/n_steps;
    
    if t > t_init && t < t_end
        
        index = ceil((t - t_init)/step_time);
        signal_value = amp(index);
        
    else
        
        signal_value = 0;
        
    end
    
end
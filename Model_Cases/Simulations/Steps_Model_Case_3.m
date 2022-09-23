clc; clear all; close all;

amp = [2 5 3 4 1]; % Values of amplitude steps

t_init = 1;
t_end = 5000;

t = t_init:t_end;
y = zeros(1,length(t));

for i = 1:length(t)
    
    y(i) = Signal(t(i), amp, t_init, t_end);
    
end

figure
plot(t,y,'b', 'LineWidth', 2)
title('Step Perturbation')
xlabel('t(s)','Interpreter','latex','FontWeight','bold','FontSize',14)
ylabel('$\dot S (\mu M)$','Interpreter','latex','FontWeight','bold','FontSize',14)
grid on

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
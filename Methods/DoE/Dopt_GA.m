% function [D_opt,Dnorm_opt,Deff_max,X_opt,L_opt,rvar_beta,VIF] = Dopt_CE(n,k,levels,model,Lb,Ub,Max_iter,seeds)

% Calculo del diseño óptimo de experimentos (sin Dummy variables)

%% INPUTS %%

clear all; clc; 

n = 100;
k = 6;
levels = [3,3,3,3,3,3];
model = 'quadratic';
Lb = [0,0,0,0,0,0];
Ub = [100,100,100,100,100,100];
seeds = 50;

%% OUTPUTS %%

%% PROCESO %%

% ¿Son los experimentos suficientes?

if strcmp(model,'lineal')
    
    p = k + 1;
       
elseif strcmp(model,'interactions')
    
    p = (k + 1) + k*(k - 1)/2;
    
elseif strcmp(model,'quadratic') 
    
    p = (2*k + 1) + k*(k - 1)/2;
    
else
    
    error('Unknown model')
    
end

if n <= p
    
    msg = ['Número de experimentos insuficiente. Se requieren como mínimo ',num2str(p + 1),' experimentos para el modelo ',model];
    error(msg);
    
end


% GENETIC ALGORITHM

DetX_opt = 0;
Dseed = zeros(1,seeds);
Dbest = zeros(1,seeds);
Tiempos = zeros(1,seeds);

tic

% Opciones de ga

nvars = n*k;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
intcon = 1:(n*k);
lb = ones((n*k),1);
ub = zeros((n*k),1);

for i_bound = 1:k

    ub((1:n) + n*(i_bound - 1)) = levels(i_bound);
    
end

% Cálculo

for i_seed = 1:seeds
    
    i_seed
    
    [x_opt,fval] = ga(@(x)Fobj(x,n,k,levels,model),nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon);
    
    DetX = - fval;
    L = reshape(x_opt,[n,k]);
    Dseed(i_seed) = DetX;
    
    if DetX > DetX_opt
        
        DetX_opt = DetX;
        L_opt = L;
        X_opt = GeneraX(L_opt,levels,model);
        
    end
    
    Dbest(i_seed) = DetX_opt;
    Tiempos(i_seed) = toc;
    
    if DetX_opt == (n^p)  % Diseño completamente ortogonal
        
        disp('DISEÑO ORTOGONAL ALCANZADO')
        break
        
    end
end

[D_opt,Dnorm_opt] = GeneraD(L_opt,levels,Lb,Ub);
Deff_seed = ((Dseed).^(1/p))/n;
Deff_opt = ((Dbest).^(1/p))/n;
Deff_max = Deff_opt(seeds)

rvar_beta = diag(inv(transpose(X_opt)*X_opt));
VIF = n*rvar_beta;






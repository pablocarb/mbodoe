function [D_opt,Dnorm_opt,Deff_max,X_opt,L_opt,rvar_beta,VIF] = Dopt_GA(n,k,levels,model,Lb,Ub,seeds)

% Calculo del diseño óptimo de experimentos (sin Dummy variables), con
% algoritmo genético

%% INPUTS %%

%% OUTPUTS %%

%% PROCESO %%

% ¿Son los experimentos suficientes?

if strcmp(model,'linear')
    
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
options = optimoptions('ga','Display','None');

for i_bound = 1:k

    ub((1:n) + n*(i_bound - 1)) = levels(i_bound);
    
end

% Cálculo

for i_seed = 1:seeds
    
    [x_opt,fval] = ga(@(x)Fobj(x,n,k,levels,model),nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
    
    DetX = - fval;
    L = reshape(x_opt,[n,k]);
    Dseed(i_seed) = DetX;
    
    if DetX > DetX_opt
        
        DetX_opt = DetX;
        L_opt = L;
        X_opt = GeneraX(L_opt,levels,model);
        
    end
    
    Dbest(i_seed) = DetX_opt;
    
    if DetX_opt == (n^p)  % Diseño completamente ortogonal
        
        disp('DISEÑO ORTOGONAL ALCANZADO')
        break
        
    end
end

[D_opt,Dnorm_opt] = GeneraD(L_opt,levels,Lb,Ub);
Deff_seed = ((Dseed).^(1/p))/n;
Deff_opt = ((Dbest).^(1/p))/n;
Deff_max =  max(Deff_opt);

rvar_beta = diag(inv(transpose(X_opt)*X_opt));
VIF = n*rvar_beta;

%% FUNCIONES %%

function fun = Fobj(x,n,k,levels,model)

    L_design = reshape(x,[n,k]);
    X = GeneraX(L_design,levels,model);
    fun = -det(transpose(X)*X);

end

end


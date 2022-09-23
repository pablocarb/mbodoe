function [D_opt,Dnorm_opt,Deff_max,X_opt,L_opt,rvar_beta,VIF] = Dopt_CE(n,k,levels,model,Lb,Ub,Max_iter,seeds)

% Calculo del diseño óptimo de experimentos (sin Dummy variables)

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


% COORDINATE EXCHANGE

DetX_opt = 0;
Dseed = zeros(1,seeds);
Dbest = zeros(1,seeds);
Tiempos = zeros(1,seeds);

tic

for i_seed = 1:seeds
    
% - Generación aleatoria de primer diseño + evaluación

    L = zeros(n,k);

    for i = 1:k

        L(:,i) = randi([1,levels(i)],n,1);

    end

    X = GeneraX(L,levels,model);
    DetX = evalX(X);

    % - Aplicación de algoritmo

    for iter = 1:Max_iter

        L_prev = L;

        for j_fac = 1:k     % Se recorre por columnas

            n_levels = levels(j_fac);
            evaluaciones = zeros(1,n_levels);

            for i_exp = 1:n

                ind_ant = L(i_exp,j_fac);
                evaluaciones(ind_ant) = DetX;   %Solución retenida

                for i_level = 1:n_levels

                    if i_level == ind_ant

                        % Ya calculado, no hace falta volver a calcular

                    else

                        L(i_exp,j_fac) = i_level;
                        X = GeneraX(L,levels,model);
                        evaluaciones(i_level) = evalX(X);

                    end

                end

                [DetX,level_opt] = max(evaluaciones);
                L(i_exp,j_fac) = level_opt;

            end

        end

        if L == L_prev

            break

        end
    end

    % Matrices de diseño
    
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
Deff_max = max(Deff_opt);

rvar_beta = diag(inv(transpose(X_opt)*X_opt));
VIF = n*rvar_beta;

%% FUNCIONES %%

function DetX = evalX(X)

    DetX = det(transpose(X)*X);

end

end


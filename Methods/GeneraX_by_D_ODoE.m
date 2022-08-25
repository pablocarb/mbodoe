function X = GeneraX_by_D_ODoE(D,model)

% Función que genera la matriz X a partir de la matriz de diseño de
% experimentos D, según el modelo seleccionado

n = size(D,1);   % Number of experiments
k = size(D,2);   % Number of factors

if strcmp(model,'lineal')
    
    X = [ones(n,1) D];
       
elseif strcmp(model,'interactions')
    
    X = [ones(n, 1) D zeros(n, k*(k - 1)/2)];
    
    term = k + 1;
    
    for i = 1:k
        for j = (1 + i):k

            term = term + 1;
            X(:,term) = D(:,i).*D(:,j); 

        end
    end
    
elseif strcmp(model,'quadratic') 
    
    X = [ones(n, 1) D zeros(n, k*(k - 1)/2) (D.^2)];
    
    term = k + 1;
    
    for i = 1:k
        for j = (1 + i):k

            term = term + 1;
            X(:,term) = D(:,i).*D(:,j); 

        end
    end
    
else
    
    error('Unknown model')
    
end

end

function X = GeneraX(L,levels,model)

% Función que genera la matriz X a partir de la asignación discreta de
% niveles, L

    n = size(L,1);
    k = size(L,2);

    if strcmp(model,'lineal')
        
        p = k + 1;
        
        X = zeros(n,p);
        X(:,1) = ones(n,1);

        for i = 1:k

            X(:,i + 1) = -1 + (2/(levels(i) - 1))*(L(:,i) - 1);

        end
        
    elseif strcmp(model,'interactions')
        
        p = (k + 1) + k*(k - 1)/2;
        
        X = zeros(n,p);
        X(:,1) = ones(n,1);

        for i = 1:k

            X(:,i + 1) = -1 + (2/(levels(i) - 1))*(L(:,i) - 1);

        end
        
        term = k + 1;
    
        for i = 1:k
            for j = (1 + i):k
            
                term = term + 1;
                X(:,term) = X(:,i + 1).*X(:,j + 1); 
            
            end
        end
        
    elseif strcmp(model,'quadratic')
        
        p = (2*k + 1) + k*(k - 1)/2;
        
        X = zeros(n,p);
        X(:,1) = ones(n,1);

        for i = 1:k

            X(:,i + 1) = -1 + (2/(levels(i) - 1))*(L(:,i) - 1);
            X(:,(k + 1) + (k*(k - 1)/2) + i) = X(:,i + 1).^2;

        end
           
        term = k + 1;
    
        for i = 1:k
            for j = (1 + i):k
            
                term = term + 1;
                X(:,term) = X(:,i + 1).*X(:,j + 1); 
            
            end
        end 
        
    end

end

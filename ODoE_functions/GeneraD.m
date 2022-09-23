function [D,Dnorm] = GeneraD(L,levels,Lb,Ub)

% Función que genera la matriz de diseño D a partir de la asignación 
% discreta de niveles, L
    
    n = size(L,1);
    k = size(L,2);
    
    Dnorm = zeros(n,k);
    D = zeros(n,k);
    
    M = (Lb + Ub)/2;
    Delta = (Ub - Lb)/2;
    
    for i = 1:k

        Dnorm(:,i) = -1 + (2/(levels(i) - 1))*(L(:,i) - 1);
        D(:,i) = (Dnorm(:,i)*Delta(i)) + M(i);

    end

end
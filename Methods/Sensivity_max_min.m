function [fp_dot_min, fp_dot_max] = Sensivity_max_min(fun_y, dim_x, dim_p, Ax, bx, Aeqx, beqx, Lbx, Ubx, Ap, bp, Aeqp, beqp, Lbp, Ubp)
    
    fp_dot_min = zeros(1, dim_p);
    fp_dot_max = zeros(1, dim_p);
    
    options = optimoptions('ga','Display','None','UseParallel', true, 'UseVectorized', false);
    
    for i = 1:dim_p
        
        % Previous definitions
        
        nvars = dim_x + dim_p;
        A = [Ax, zeros(size(Ax,1),dim_p); zeros(size(Ap,1,dim_x)), Ap];
        b = [bx;bp];
        Aeq = [Aeqx, zeros(size(Aeqx,1),dim_p); zeros(size(Aeqp,1,dim_x)), Aeqp];
        beq = [beqx;beqp];
        lb = [Lbx, Lbp];
        ub = [Ubx, Ubp];
        
        % Get the minimum value
        
        fun = @(x)fobj_min(x, fun_y, dim_x, dim_p, i);
        
        % Optimization
        
        x_min = ga(fun,nvars,A,b,Aeq,beq,lb,ub,[],options);
        fval_min = fun(x_min);
        fp_dot_min(i) = fval_min;
        
        % Get the maximum value
        
        fun = @(x)fobj_max(x, fun_y, dim_x, dim_p, i);
        
        % Optimization
        
        x_max = ga(fun,nvars,A,b,Aeq,beq,lb,ub,[],options);
        fval_max = fun(x_max);
        fp_dot_max(i) = -fval_max;
        
    end
    

    function val = fobj_min(x, fun_y, dim_x, dim_p, i_s)
        
        point = x(1:dim_x);
        p = x(dim_x + 1: dim_x + dim_p);
        
        h_p = 1e-6*p(i_s);      % Paso relativo al orden de los parámetros
        
        p_h = p;
        p_mh = p;
        
        p_h(i_s) =  p(i_s) + h_p;
        p_mh(i_s) =  p(i_s) - h_p;
        
        y_h = fun_y(point, p_h);
        y_mh = fun_y(point, p_mh);
        
        val = (y_h - y_mh)/(2*h_p);
        
    end

    function val = fobj_max(x, fun_y, dim_x, dim_p, i_s)
        
        point = x(1:dim_x);
        p = x(dim_x + 1: dim_x + dim_p);
        
        h_p = 1e-6*p(i_s);      % Paso relativo al orden de los parámetros
        
        p_h = p;
        p_mh = p;
        
        p_h(i_s) =  p(i_s) + h_p;
        p_mh(i_s) =  p(i_s) - h_p;
        
        y_h = fun_y(point, p_h);
        y_mh = fun_y(point, p_mh);
        
        val = -((y_h - y_mh)/(2*h_p));
        
    end


end


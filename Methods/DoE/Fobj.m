function fun = Fobj(x,n,k,levels,model)

    L = reshape(x,[n,k]);
    X = GeneraX(L,levels,model);
    fun = -det(transpose(X)*X);

end
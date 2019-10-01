function ModJError = FunctionModJError(N,i)
    for n = 1:N-1
        J(n) = (2*sqrt(n*(N-n)))/(sqrt((N^2)-1));
    end
    r = normrnd(0,i,[1,N-1]);
    ModJError = J + r; %At the moment the errors are just being added..
    
end

 
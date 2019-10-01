function ErrorHamiltonian = FunctionErrorHamiltonian3(N,JError,BError)
    ModJError = FunctionModJError(N,JError);
    ModBError = FunctionModBError(N,BError);
    C = diag((ModJError),-1); %Create the segments of repeated elements required
    D = diag((ModJError),+1);
    E = diag(ModBError);
    h = C + D + E;
    [V,D] = eig(h)
    I = eye(N);
    H = -kron(h,I)+kron(I,h);
    ErrorHamiltonian = blkdiag(0,h,-h,H);
    
end



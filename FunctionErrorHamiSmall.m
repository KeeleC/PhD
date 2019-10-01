function ErrorHamiSmall = FunctionErrorHamiSmall(N,JError,BError)
    ModJError = FunctionModJError(N,JError);
    ModBError = FunctionModBError(N,BError);
    C = diag((ModJError),-1); %Create the segments of repeated elements required
    D = diag((ModJError),+1);
    E = diag(ModBError);
    ErrorHamiSmall = C + D + E;
end
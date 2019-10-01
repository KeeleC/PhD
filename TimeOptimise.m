clear
%THIS CODE WORKS%
%For N = 3, 4 iterations completed per second
%For N = 5, 2 iterations per second
%For N = 7, ~1 iterations per second
%For N = 9, 0.5 iterations per second



%% INITIAL PARAMETERS
N = 3; %Number of sites
InitialSite = 1; %Position of initial excitation
FinalSite = N; %Position to extract excitation
JMax = 1
TStart = 2.; %Time gap for highest fidelity
TEnd = 2.5;
Gamma1 = 0.5 %Errors for dephasing/decoherence
Gamma2 = 0.0
Runs = 400 %Iterations



%% CREATION OF INITIAL STATE, FINAL STATE AND STANDARD HAMILTONIAN
%Q = FunctionQMatrix(N,Gamma2); %Noisy matrix amplitude damping
%AllNoise = M + Q;
DM = FunctionDensityMatrix(N,InitialSite); %Initial density matrix
FDM = FunctionFinalDensityMatrix(N,FinalSite); %Final desired state
H = FunctionHamiltonian3(N); %Standard Hamiltonian with no errors




%% HAMILTONIAN WITH ERRORS BUILT ACCORDING TO ERROR
for i = 1:Runs %run miltiple times to average
    i
    for B = 1:51 %maximum error on site
        BError = (B-1)*0.01;
        B1(B) = BError;
        for J = 1:21 %maximum error for couplings
            JError = (J-1)*0.01;
            J1(J) = JError;
            EH = FunctionErrorHamiltonian(N,JError,BError); %Create H with errors
            
            
            
            %% ADD NOISE
            %M = FunctionMatrixM(N,Gamma1,EH);
            %assert(ishermitian(M))
            %Q = FunctionQMatrix(N,Gamma2,ErrorHamiltonian);
            %assert(ishermitian(Q))
            %MQ = M + Q;
            %Fidelity = FunctionFindFidelity(T1,T2,DM,FDM,H);
            
            
            
            
            %% CONVERGING ON HIGHEST FIDELITY
            psi1 = expm(-1i*EH*TStart)*DM; %Create system at each time
            psi2 = expm(-1i*EH*TEnd)*DM;
            F1 = (abs(FDM'*psi1))^2; %Fidelity at the 2 times
            F2 = (abs(FDM'*psi2))^2;
            while abs(F1-F2) > 0.000005 %Tolerance
                if F1 > F2 
                    TEnd = TEnd-(TEnd-TStart)/3;
                    psi2 = expm(-1i*EH*TEnd)*DM;
                    F2 = (abs(FDM'*psi2))^2;
                elseif F2 > F1
                    TStart = TStart+(TEnd-TStart)/3;
                    psi1 = expm(-1i*EH*TStart)*DM;
                    F1 = (abs(FDM'*psi1))^2;
                end
            end
            Fidelity = F1;
            A(J,B,i) = Fidelity;
            
            
            %%ASSERTIONS
           assert(1-1e-20 <= A(1,1,i) <= 1+1e-20);
           assert(A(J,B,i) <= 1);
            
        end
    end
end




%% PLOTTING
levels = [0.8,0.85,0.90,0.95,0.96,0.97,0.98,0.99]  
Y = quantile(A,0.75,3);
% compute 75th percentile (third quartile)
A1 = mean(A,3);
surf(Y)
%[C,h] = contour(Y,20);
%w = h.LineWidth;
%h.LineWidth = 2;
%clabel(C,h)
colorbar
xlabel('BError')
ylabel('JError')
Y = -1:2:J;
X = -1:2:B;
xtickangle(-90)
yticks([Y]);
xticks([X]);
yticklabels({(Y./100)-0.01});
xticklabels({(X./100)-0.01});



             

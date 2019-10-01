
%For N = 3, ~1,2 iterations/sec
%For N = 5, ~0.86 iterations/sec
%THIS CODE WORKS%
clear 
N = 35; %Number of sites
p = 1; %Position of initial excitation
q = N; %Position we want to extract excitation
JMax = 1
BError = 0.5;
JError = 0.5;  
Gamma1 = 0.1;
TStart = 1.52
TEnd = 1.6

DM = [1;0;0;0]
FDM = [0,0,0,1]
Hami1 = [0,3*sqrt(11),0,0;3*sqrt(11),0,98,0;0,98,0,3*sqrt(11);0,0,3*sqrt(11),0]
DensityMatrix = FunctionDensityMatrix(N,p); %Initial density matrix
FinalDensityMatrix = FunctionFinalDensityMatrix(N,q); %Final desired state
%H = FunctionHamiltonian3(N); %Standard Hamiltonian with no errors   
ErrorHamiltonian = FunctionErrorHamiltonian(N,JError,BError)%Create H with errors
%HamiExactError = FunctionHamiExactError(N,JError,BError);
%M = FunctionMatrixM(N,Gamma1,ErrorHamiltonian);
for t = 1:3000
    T = (t-1)*0.001
    psi = expm(-1i*Hami1.*T)*[1;0;0;0];
    F = [0,0,0,1]*psi;%Fidelity calculation %%This takes a while
    Fidelity = (abs(F))^2;
    plot(T,Fidelity,'.'); hold on
end

 psi1 = expm(-1i*Hami1*TStart)*DM; %Create system at each time
 psi2 = expm(-1i*Hami1*TEnd)*DM;
            F1 = (abs(FDM*psi1))^2; %Fidelity at the 2 times
            F2 = (abs(FDM*psi2))^2;
            while abs(F1-F2) > 0.00000005 %Tolerance
                if F1 > F2 
                    TEnd = TEnd-(TEnd-TStart)/3;
                    psi2 = expm(-1i*Hami1*TEnd)*DM;
                    F2 = (abs(FDM*psi2))^2;
                elseif F2 > F1
                    TStart = TStart+(TEnd-TStart)/3;
                    psi1 = expm(-1i*Hami1*TStart)*DM;
                    F1 = (abs(FDM*psi1))^2;
                end
           end
            Fide= F1
            TR = TStart
            TR2 = TEnd
         
            

  


clear
%THIS CODE WORKS%

%% INITIAL PARAMETERS
N = 5; %Number of sites
nE = 2; %Encoding and decoding sites
nD = 2;
InitialSite = 1; %Position of initial excitation
FinalSite = N; %Position to extract excitation
JMax = 1
Gamma1 = 0.5 %Errors for dephasing/decoherence
Gamma2 = 0.0
Runs = 2000 %Iterations


OptimalTime = (pi/4)*(sqrt(N^(2)-1)); %Optimal time for PST chain

%% HAMILTONIAN WITH ERRORS BUILT ACCORDING TO ERROR
for p = 1:Runs %run miltiple times to average
    p
    for B = 1:21 %maximum error on site
        BError = (B-1)*0.1;
        B1(B) = BError;
        for J = 1:21 %maximum error for couplings
            JError = (J-1)*0.1;
            J1(J) = JError;
            ErrorHamiltonian = FunctionErrorHamiSmall(N,JError,BError); %Create H with errors
            
            
            
             %% CONVERGING ON HIGHEST FIDELITY
            [Eigenvecs,Eigenvals]= eig(ErrorHamiltonian);
            for i = 1:N
                EvolvedState(:,i) = Eigenvecs(1,i)*Eigenvecs(:,i)*(exp(Eigenvals(i,i)*(-1i)*OptimalTime)); %Calculate evolution at optimal time
            end 
            for i = 1:N
                DesiredState(:,i) = Eigenvecs(N,i)*Eigenvecs(:,i); %Desired state
            end
            Fidelity1 = sum(EvolvedState,2);
            Desired1 = sum(DesiredState,2);
            Fidelity1 = abs(dot(Fidelity1,conj(Desired1)))^2;
            
         
            AverageFidelity = (1/6)*(3+2*sqrt(Fidelity1)+Fidelity1);
            A(J,B,p) = AverageFidelity;
            MaxSingularValue(J,B,p) = FunctionEncode(ErrorHamiltonian,nE,nD,N,OptimalTime);
          
        end
        
            
        %%ASSERTIONS
        assert(1-1e-20 <= A(1,1,p) <= 1+1e-20);
        assert(A(J,B,p) <= 1);
        assert(1-1e-20 <= MaxSingularValue(1,1,p) <= 1+1e-20);
        assert(MaxSingularValue(J,B,p) <= 1);
        
    end
end

    
%% PLOTTING
%SURFACE PLOT 3D
figure(1)
Y = quantile(A,0.75,3);
surf(Y)

figure(2)
YE = quantile(MaxSingularValue,0.75,3);
surf(YE)



%CONTOUR PLOT 2D
figure(3)
levels = [0.4,0.5,0.505,0.51,0.52,0.53,0.54,0.55,0.6,0.7,0.8,0.85,0.90,0.95,0.96,0.97,0.98,0.99,1]  
Y = quantile(A,0.75,3);
%compute 75th percentile (third quartile)
A1 = mean(A,3);
[C,h] = contour(Y,25);
s = h.LevelList
h.LevelList = levels
w = h.LineWidth;
h.LineWidth = 2;
clabel(C,h);
colorbar;
LineColor = 'flat'
xlabel('BError')
ylabel('JError')
Y = -1:2:J;
X = -1:2:B;
xtickangle(-90)
yticks([Y]);
xticks([X]);
yticklabels({(Y./100)-0.01});
xticklabels({(X./100)-0.01});
grid

figure(4)
levels = [0.4,0.5,0.505,0.51,0.52,0.53,0.54,0.55,0.6,0.7,0.8,0.85,0.90,0.95,0.96,0.97,0.98,0.99,1]  
YE = quantile(MaxSingularValue,0.75,3);
%compute 75th percentile (third quartile)
A1E = mean(MaxSingularValue,3);
[C,h] = contour(YE,25);
s = h.LevelList
h.LevelList = levels
w = h.LineWidth;
h.LineWidth = 2;
clabel(C,h);
colorbar;
LineColor = 'flat'
xlabel('BError')
ylabel('JError')
Y = -1:2:J;
X = -1:2:B;
xtickangle(-90)
yticks([Y]);
xticks([X]);
yticklabels({(Y./100)-0.01});
xticklabels({(X./100)-0.01});
grid
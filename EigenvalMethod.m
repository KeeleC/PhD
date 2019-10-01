clear
%THIS CODE WORKS%

%% INITIAL PARAMETERS
N = 21; %Number of sites
InitialSite = 1; %Position of initial excitation
FinalSite = N; %Position to extract excitation
JMax = 1
TStart = 16.; %Time gap for highest fidelity
TEnd = 17.;
Gamma1 = 0.5 %Errors for dephasing/decoherence
Gamma2 = 0.0
Runs = 1000 %Iterations
Hami = [0,1,0;1,0,1;0,1,0];


%% HAMILTONIAN WITH ERRORS BUILT ACCORDING TO ERROR
for p = 1:Runs %run miltiple times to average
    p
    for B = 1:51 %maximum error on site
        BError = (B-1)*0.01;
        B1(B) = BError;
        for J = 1:21 %maximum error for couplings
            JError = (J-1)*0.01;
            J1(J) = JError;
            EH = FunctionErrorHamiSmall(N,JError,BError); %Create H with errors
            
            
            
         
            
            %% CONVERGING ON HIGHEST FIDELITY
            [Eigenvecs,Eigenvals]= eig(EH);
            for i = 1:N
                F1(i) = Eigenvecs(1,i)*Eigenvecs(N,i)*(exp(Eigenvals(i,i)*(-1i)*TStart));
                F2(i) = Eigenvecs(1,i)*Eigenvecs(N,i)*(exp(Eigenvals(i,i)*(-1i)*TEnd));
            end
            Fidel1 = sum(F1);
            Fidelity1 = Fidel1*conj(Fidel1);
            Fidel2 = sum(F2);
            Fidelity2 = Fidel2*conj(Fidel2);
            while abs(Fidelity1-Fidelity2) > 0.000005 %Tolerance
                if Fidelity1 > Fidelity2 
                    TEnd = TEnd-(TEnd-TStart)/3;
                    for i = 1:N
                        F2(i) = Eigenvecs(1,i)*Eigenvecs(N,i)*(exp(Eigenvals(i,i)*(-1i)*TEnd));
                    end
                    Fidel2 = sum(F2);
                    Fidelity2 = Fidel2*conj(Fidel2);
                elseif Fidelity2 > Fidelity1
                    TStart = TStart+(TEnd-TStart)/3;
                    for i = 1:N
                        F1(i) = Eigenvecs(1,i)*Eigenvecs(N,i)*(exp(Eigenvals(i,i)*(-1i)*TStart));
                    end
                    Fidel1 = sum(F1);
                    Fidelity1 = Fidel1*conj(Fidel1);
                end
            end
            Fidelity = Fidelity1;
            A(J,B,p) = Fidelity;
        end
            
        %%ASSERTIONS
        assert(1-1e-20 <= A(1,1,p) <= 1+1e-20);
        assert(A(J,B,p) <= 1);
    end
end

    
%% PLOTTING
%SURFACE PLOT 3D
%Y = quantile(A,0.75,3);
%surf(Y)

%CONTOUR PLOT 2D
levels = [0.8,0.85,0.90,0.95,0.96,0.97,0.98,0.99]  
Y = quantile(A,0.75,3);
% compute 75th percentile (third quartile)
A1 = mean(A,3);
[C,h] = contour(Y,20);
w = h.LineWidth;
h.LineWidth = 2;
clabel(C,h)
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

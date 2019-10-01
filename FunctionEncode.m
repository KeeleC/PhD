function MaxS = FunctionEncode(Hamiltonian,nE,nD,N,Time);

%% Parameters

    N = 301;
    JError = 0;
    BError = 0;
    Hamiltonian = FunctionErrorHamiSmall(N,JError,BError); 
    UniHami = FunctionUniHami(N);
    %   
    nE = 20;
    nD = 20;
    Time = (pi/4)*(sqrt(N^(2)-1));
% 
% %% Finding highest singular value
%     %Take an average over these for longer time (maybe 20 runs for 5000
%     %time steps?
    IdenE = eye(nE);
    IdenD  = eye(nD);
    Enc = padarray(IdenE,[0 N-nE],0,'post');
    Dec = padarray(IdenD,[N-nD 0],0,'pre');
% 
%     
%     EncHami = expm((-1i)*Hamiltonian*Time);
%     ProjectEnc = EncHami*Enc';
%     ProjectEnc2 = Dec'*ProjectEnc;
%     [~,SEnc,~] = svd(ProjectEnc2);
%     SEnc;
%     S = SEnc(1,1);
%     assert(1-1e-20 <= S <= 1+1e-20);
%     MaxS = (1/6)*(3+2*sqrt((S^2))+(S^2));
%     assert(1-1e-20 <= MaxS <= 1+1e-20);


%% For plotting against time

%     T1 = 0:0.01:50;
%     B = bsxfun(@times,UniHami,reshape(T1,1,1,[]));
%     NewM = -1i*B;
%     for i = 1:length(T1);
%         NewMM(:,:,i) = expm(NewM(:,:,i));
%     end
%     NewM2 = NewMM((N-nE)+1:N,1:nD,:);
%     for i = 1:length(T1);
%         i;
%         [U1(:,:,i),S1(:,:,i),V1(:,:,i)] = svd(NewM2(:,:,i));
%         Max(i) = S1(1,1,i);
%         MaxS(i) = (1/6)*(3+2*sqrt((Max(i)^2))+(Max(i)^2));
%         I = (i-1)*0.01;
%         plot(I,MaxS(i),'.b'); hold on
%     end

%% Another way to plot against time

    for t = 1:6000
            T = (t-1)*0.1;
            Unitary(:,:,t) = expm(-1i*UniHami*T); %Needs to be multiplied by t!
            Project1(:,:,t)  = Unitary(:,:,t)*Enc';
            Project(:,:,t)  = Dec'*Project1(:,:,t);
            [U(:,:,t) ,S(:,:,t) ,V(:,:,t) ] = svd(Project(:,:,t));
            S(:,:,t);
            %Encoding = padarray(V(1,:),[0 N-nE],0,'post'); %remove bracks
            %Decoding = padarray(U(:,1),[N-nD 0],0,'pre')';
            Max(t) = S(4,4,t);
            MaxS(t) = (1/6)*(3+2*sqrt((Max(t)^2))+(Max(t)^2));
            plot(T,MaxS(t),'.c'); hold on
    end
    
    [M,I]= max(MaxS);
    MaxS2 = M;
    Time = (I-1)*0.1;
end

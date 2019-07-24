function [Synch, Vm, S, D, instISI, spklist, LamRatio, Lambda2, LambdaVar] = ML_Net_Sweep(T, N, k, rho, I, SynStren, coefs, Syns);

rand('seed',sum(100*clock));
randn('seed',sum(100*clock));

load('ber_3000_0.010_4.0_1.4_1.3_1.2.dat');
Syns=cell(3000,1);
for i=1:1:3000
    Syns{i}=find(ber_3000_0_010_4_0_1_4_1_3_1_2(i,:));
end

if(nargin<1)
    T = 1000;
end;
if(nargin<2)
    N=3000;
end;
if(nargin<3)
    k=30;
end;
if(nargin<4)
    rho=0.01;
end;

if(nargin<5)
    I=0;
end;

if(nargin<6)
    SynStren=0.01;
end;
dt = 0.02;
NumSamples = T/dt;
NoiseScalar = (sqrt(dt)*0.2);

% 
% if(nargin<8)
%     Neighbors = ([-k/2:-1,1:k/2]);
%     CellNum = 1:N;
%     Syns = repmat(Neighbors,N,1) + repmat(CellNum',1,k);
%     ind = find(Syns<1);
%     Syns(ind)=Syns(ind)+N;
%     ind = find(Syns>N);
%     Syns(ind)=Syns(ind)-N;
%     
%     Num2Rewire = ceil(N*k*rho);
%     
%     randselect = randperm(N*k);
%     randsyns = randselect(1:Num2Rewire);
%     randtargets = ceil(rand(Num2Rewire,1)*N);
%     Syns(randsyns)=randtargets;
% end;


lambda2=0;
LambdaN=0;
LambdaVar=0;
if(0) %set to 1 to calculate eigen value spectrum
    if(iscell(Syns))
        K=length(Syns{1});
        CellNum = [];
        PostSyns=[];
        for i=1:N
            K=length(Syns{i});
            CellNum=[CellNum; i*ones(K,1)];
            PostSyns = [PostSyns; Syns{i}'];
        end;
         SparseSyns=sparse(CellNum, PostSyns,ones(length(CellNum),1));
         

    else
        
        CellNum=repmat(1:N,1,k);
        SparseSyns=sparse([1:N, CellNum], [1:N, Syns(:)'],[-30*ones(N,1); ones(k*N,1)], N, N);
        SparseSyns=sparse(CellNum, Syns(:), ones(k*N,1), N, N);
        SparseSyns(find(SparseSyns>1))=1;
    end;
    SparseSyns=SparseSyns-sparse(diag(sum(SparseSyns,1)));
  
    %[V,D] = eigs(SparseSyns,2,'SM');
    % DS(2,2)
    %tic
    %[V,D]=eig(FullSyn);
    %toc
%     subplot 212
%     plot(diag(D),'.')
%     diag(D)
%     lambda2=D(2,2);
    
    D = eig(full(SparseSyns));
    figure
    plot(D,'.');
    xlabel('Real(\lambda)');
    ylabel('Imag(\lambda)');
    LambdaVar = sum(abs(D(2:end)-mean(D(2:end))).^2)/(mean(D(2:end))^2*(N-1));
    %LambdaN = D(1);
    %[V,LambdaN]=eigs(SparseSyns,1);
    LamRatio=min(real(D(abs(D)>1e-12)))/max(real(D(abs(D)>1e-12)))
    Lambda2=max(real(D(abs(D)>1e-12)));
end;

if(nargin<7)
    coefs.C=(1);
    coefs.gL=(8); %8
    %coefs.gL_i=8; %8; %1.0;
    coefs.EL=(-53.23878);%-79.5 %-78.0; Tease this parameter down(up) to slow(speed) activation of excitatory neurons
    %EL_i=-78;
    coefs.gNa= (18.22315);%20;
    %gNa_i=15;%20; %20; %4.0;
    coefs.ENa=(60); %60.0;
    coefs.gK=(4); %4.0;
    %gK_i=10; %10; %4.0;
    coefs.EK=(-95.52116); %-90.0;
    coefs.Vhalfm=(-7.37257); %-30.0;
    coefs.km=(11.97239); %7.0;
    coefs.Vhalfn=(-16.34877); %-45.0;
    %Vhalfn_i=-42.2;
    coefs.kn=(4.21133); %5.0;
    coefs.tau=(1); %1.0;
    %MaxRand = 5 %0.57;%3.9;
    coefs.E_EPSP = (0);
    coefs.tauEPSPr = (0.25); %2.63;
    coefs.tauEPSPf  = (.5); %6.21;
    coefs.gEPSP=SynStren/( coefs.tauEPSPf - coefs.tauEPSPr);
end;


[t,vinit] = ode23(@ML_derivs,[0:0.01:100],[-60, 0, 0, 0],[],I, coefs);
[ind, spkind] = findpeaks(vinit(:,1),'minpeakheight',-20);
PredNumSpikes = length(ind)*N*1.2*T/100;
ISI = 100/length(ind);
ExpSynCurrent = 1/ISI*k*SynStren*(-mean(vinit(:,1)));

t = t(spkind(end-1):spkind(end));
t=t-t(1);
vinit=vinit(spkind(end-1):spkind(end),:);

ind = ceil((length(t))*rand(N,1));
V=(vinit(ind,1));
n=(vinit(ind,2));
lastSpkTime = -t(ind);
Ser = zeros(N,1);
Sef = zeros(N,1);
%y = [-60, 0, 0, 0];
%y = repmat(y,N,1);
%y(:,1) = y(:,1) + randn(N,1)*2;
%y = y';


if(nargout>5)
    LengthSpkList = PredNumSpikes;
    spklist = zeros(LengthSpkList,2);
end

Vm = zeros(NumSamples,1);
S = zeros(NumSamples,1);
Synch = zeros(NumSamples,1);

spknum=0;

t=0;
dVdt = zeros(N,1);
dndt = zeros(N,1);
dssdt = zeros(N,1);
dsfdt = zeros(N,1);
lastV = zeros(N,1);

ISI = ISI*ones(N,1);
D = zeros(T/dt,1);
instISI = zeros(T/dt,1);

%coefs = (coefs);

thresh = -20;
tempV=linspace(coefs.EK, coefs.ENa, 1024);
tempV = tempV';
NInfLUT = 1.0./(1.0+exp((coefs.Vhalfn - tempV)/coefs.kn ));
NTauLUT = exp(-0.07*tempV-3);
MInfLUT = (1.0./(1.0 + exp((coefs.Vhalfm - tempV)/coefs.km )));
LUTLength = length(MInfLUT);

PhaseLUT = exp(1i*linspace(0,2*pi,1024));
PhaseLUTLen=length(PhaseLUT);

%tic;
SynStrenSweep=linspace(0.01,0.01,NumSamples);
Tdrive=ones(1,ceil(NumSamples/3))*8;
I=linspace(8,-4,floor(2*NumSamples/3));
I=horzcat(Tdrive,I);

for i=1:NumSamples;
    %dydt = ML_derivs(t, y, I, coefs);
%     V = y(1:4:end);
%     n = y(2:4:end);
%     Ser = y(3:4:end);
%     Sef = y(4:4:end);
    coefs.gEPSP=SynStren/( coefs.tauEPSPf - coefs.tauEPSPr);
    Vind = ceil(LUTLength*(V-coefs.EK)/(coefs.ENa-coefs.EK));
    dVdt = 1.0./coefs.C*(I(i)-ExpSynCurrent-coefs.gL.*(V-coefs.EL) - ...
    coefs.gNa.*MInfLUT(Vind).*(V-coefs.ENa) - ...
    coefs.gK.*n.*(V-coefs.EK) + ...
    coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V));%(coefs.E_EPSP-V))
    dndt = (NInfLUT(Vind)-n)./NTauLUT(Vind);
    dssdt = -Ser/coefs.tauEPSPr;
    dsfdt = -Sef/coefs.tauEPSPf;

    %lasty = y;
    %y = y + dydt*dt;
    lastV = V;
    V = V + dVdt*dt + (randn(N,1))*NoiseScalar;
    n = n + dndt*dt;
    Ser = Ser + dssdt*dt;
    Sef = Sef + dsfdt*dt;
    
    spks = find((V>thresh).*(lastV<thresh));
    if(~isempty(spks))
        ISI(spks)=t-lastSpkTime(spks);
        lastSpkTime(spks) = t;
       
        for spk=1:length(spks);
            if(iscell(Syns))
                Ser(Syns{spks(spk),:})=Ser(Syns{spks(spk),:})+1;
                Sef(Syns{spks(spk),:})=Sef(Syns{spks(spk),:})+1;
            else
                Ser(Syns(spks(spk),:))=Ser(Syns(spks(spk),:))+1;
                Sef(Syns(spks(spk),:))=Sef(Syns(spks(spk),:))+1;
            end;
        end;
        if(nargout>5)
            spklist(spknum+1:spknum+1+length(spks),1)=t;
            spklist(spknum+1:spknum+length(spks),2)=spks;
            spknum=spknum+length(spks);
        end;
    end;
    Vm(i)=V(1);
    S(i)=(coefs.gEPSP.*(Sef(1)-Ser(1)).*(coefs.E_EPSP-V(1)));
    %Synch(i) = abs(mean(exp(j*2*pi*(t-lastSpkTime)/mean(ISI))));
    Phase = mod(ceil((PhaseLUTLen-1)*(t-lastSpkTime)/mean(ISI)),PhaseLUTLen)+1;
    Synch(i) = abs(mean(PhaseLUT(Phase)));
    D(i)=length(spks);
    instISI(i)=mean(ISI);
    t=t+dt;

end;
%toc
if(nargout>5)
    if spknum<length(spklist);
        spklist=spklist(1:spknum,:);
    end;
end;

figure;
% subplot 311;
% plot(Vm); axis tight;
% subplot 311;
% plot(S);
% ylabel('Synatpic drive');
subplot 311
%plot(spklist(:,1),spklist(:,2),'.')
plot(dt:dt:T,D,'k');
ylabel('Firing Density');

subplot 312
%plot(spklist(:,1),spklist(:,2),'.')
plot(dt:dt:T,instISI,'k');
ylabel('|ISI|');


% Calculate synchrony 
% FiducialNeuron = ceil(N/2);
% spks = find(spklist(:,2)==FiducialNeuron);
% ISI = diff(spklist(spks,1));
% Synch2=zeros(length(ISI)-1,2);
% for i=2:length(ISI);
%     ThisISI = ISI(i);
%     Synch2(i-1,1) = spklist(spks(i),1);
%     NormSpikeDiff = (spklist(spks(i-1)+1:spks(i)-1,1)-spklist(spks(i-1),1))/ISI(i);
%     Synch2(i-1,2) = abs(sum(exp(j*2*pi*NormSpikeDiff)))/length(NormSpikeDiff);
% end;
subplot 313
%plot(Synch2(:,1), Synch2(:,2),dt:dt:T, Synch);
plot(dt:dt:T, Synch,'k');
ylabel('Synch');

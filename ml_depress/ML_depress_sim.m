%function [Synch, Vm, S, D, instISI, spklist, LamRatio, Lambda2, LambdaVar] = ML_depress_sim(T, N, k, rho, I, SynStren, coefs);
function [Synch, T, dt, SyncTimeStep] = ML_depress_sim(T, N, k, rho, I, SynStren, coefs);

rand('seed',sum(100*clock));
randn('seed',sum(100*clock));

load('ber_3000_0.010_4.0_1.4_1.3_1.2.dat');
Syns=cell(3000,1);

for i=1:1:3000
    Syns{i}=find(ber_3000_0_010_4_0_1_4_1_3_1_2(i,:));
end

if(nargin<1)
    T = 200;
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
    I=-2;
end;

if(nargin<6)
    SynStren=0.02;
end;

Imag=I;

dt = 0.02;

NumSamples = T/dt;

SyncTimeStep=100;

% DepTau=1:10:100;
% FacTau=1:10:100;
% Kuramoto=zeros(length(DepTau),length(FacTau));
% Entropy=zeros(length(DepTau),length(FacTau));

% ElecStimAmp=-.9:0.1:.9;%0:0.01:.9;
% StimPeriodVect=1000./linspace(1,100,50);

Stimon=0;
Kindleon=0;
Ploton=1;
Saveon=0;
ThreeDon=0;
CurrentExternalon=0;

ElecStimAmp= 40;%linspace(0,50,50);%0:0.01:.9;
StimPeriodVect= linspace(2,4,20);%linspace(1,15,100);
CurrentVect= linspace(-5,5,10);

ElecStimAmp=40;
StimPeriodVect=5.5;%6.4,8.7
% d1=0.816;
% d2=0.575;
% f=0.917;

d1=0.7;%.7
d2=1;%1
f=0;
    %coefs.tauD1 = (380);60
    coefs.tauD1 = (30);
    %coefs.tauD2 = (9200);
    coefs.tauD2 = (920);
    %coefs.tauF = (94);30
    coefs.tauF = (90);
    
if(ThreeDon)
    Kuramoto=zeros(length(StimPeriodVect),length(ElecStimAmp),length(CurrentVect));
    Entropy=zeros(length(StimPeriodVect),length(ElecStimAmp),length(CurrentVect));
    Firing=zeros(length(StimPeriodVect),length(ElecStimAmp),length(CurrentVect));
else
    Kuramoto=zeros(length(StimPeriodVect),length(ElecStimAmp));
    Entropy=zeros(length(StimPeriodVect),length(ElecStimAmp));
    Firing=zeros(length(StimPeriodVect),length(ElecStimAmp));
    
end

if(CurrentExternalon)
    CurrentExternal=linspace(0,-6,NumSamples);
else
    CurrentExternal=zeros(NumSamples,1);
end

NoiseScalar = (sqrt(dt)*0.2);
Kindle=zeros(NumSamples,1);
Kindle_mag=4;
Kindle_dur=T/6;
Kindle_dur=Kindle_dur/dt;

if(Kindleon)
    for i=1:Kindle_dur
        Kindle(i)=Kindle_mag;
    end
end

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



% for FacTauInd=1:length(FacTau);
%     FacTauInd
%     for DepTauInd=1:length(DepTau);
%         DepTauInd

% for CurrentInd=1:length(CurrentVect);
%     CurrentInd
for FreqInd=1:length(StimPeriodVect);
    FreqInd
    
    for StrInd=1:length(ElecStimAmp);
        StrInd
        tic
        %         StimAmp=50;
        % StimFrequency=50;
        StimAmp=ElecStimAmp(StrInd);
        
        Stimvect=zeros(NumSamples,1);
        StimPeriod=StimPeriodVect(FreqInd)/dt;
        StimCount=StimPeriod;
        
        StimAmpLow=dt*StimAmp/2;
        StimDurationLow = 0; %2/dt;
        
        if(Stimon)
            for i=1:NumSamples
                if(i>StimCount)
                    Stimvect(i:i+5)=StimAmp;
                    StimCount=StimCount+StimPeriod;
                    %             for ja=i+1:i+StimDurationLow
                    %                 Stimvect(ja)=-StimAmpLow;
                    %             end
                    
                end
            end
        end
        %I=2;
        
        % [t,vinit] = ode23(@ML_derivs,[0:0.01:100],[-60, 0, 0, 0],[],CurrentVect(CurrentInd)+Kindle_mag, coefs);
        [t,vinit] = ode23(@ML_derivs,[0:0.01:100],[-60, 0, 0, 0],[],I + Kindle_mag, coefs);
        
        [ind, spkind] = findpeaks(vinit(:,1),'minpeakheight',-20);
        PredNumSpikes = length(ind)*N*1.2*T/100;
        ISI = 100/length(ind);
        ExpSynCurrent = 0; %1/ISI*k*SynStren*(-mean(vinit(:,1)));
        
        
        t = t(spkind(end-1):spkind(end));
        t=t-t(1);
        vinit=vinit(spkind(end-1):spkind(end),:);
        
        ind = ceil((length(t))*rand(N,1));
        V=(vinit(ind,1));
        n=(vinit(ind,2));
        lastSpkTime = -t(ind);
        Ser = zeros(N,1);
        Sef = zeros(N,1);
        AA = ones(N,1);
        D1 = ones(N,1);
        D2 = ones(N,1);
        F = ones(N,1);
        entropy = zeros(T/dt/SyncTimeStep,1);
        Synch = zeros(T/dt/SyncTimeStep,1);
        
        %y = [-60, 0, 0, 0];
        %y = repmat(y,N,1);
        %y(:,1) = y(:,1) + randn(N,1)*2;
        %y = y';
        
        
        
            LengthSpkList = PredNumSpikes;
            spklist = zeros(LengthSpkList,2);
        
        
        Vm = zeros(NumSamples,1);
        S = zeros(NumSamples,1);
        
        
        spknum=0;
        
        t=0;
        dVdt = zeros(N,1);
        dndt = zeros(N,1);
        dssdt = zeros(N,1);
        dsfdt = zeros(N,1);
        lastV = zeros(N,1);
        dd1dt = zeros(N,1);
        dd2dt = zeros(N,1);
        dfdt = zeros(N,1);
        
        
        
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
        
        %I=ones(NumSamples,1)*Imag;
        synchTime=0;
        scount=0;
        Term=0;
        
        for i=1:NumSamples;
            %dydt = ML_derivs(t, y, I, coefs);
            %     V = y(1:4:end);
            %     n = y(2:4:end);
            %     Ser = y(3:4:end);
            %     Sef = y(4:4:end);
            % coefs.gEPSP=SynStren/( coefs.tauEPSPf - coefs.tauEPSPr);
            Vind = ceil(LUTLength*(V-coefs.EK)/(coefs.ENa-coefs.EK));
            
            if(ThreeDon)
                dVdt = 1.0./coefs.C*(CurrentVect(CurrentInd)-ExpSynCurrent-coefs.gL.*(V-coefs.EL) - ...
                    coefs.gNa.*MInfLUT(Vind).*(V-coefs.ENa) - ...
                    coefs.gK.*n.*(V-coefs.EK) + ...
                    AA.*coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V) + ...
                    Stimvect(i)+Kindle(i)+CurrentExternal(i));%(coefs.E_EPSP-V))
            else
                dVdt = 1.0./coefs.C*(I-ExpSynCurrent-coefs.gL.*(V-coefs.EL) - ...
                    coefs.gNa.*MInfLUT(Vind).*(V-coefs.ENa) - ...
                    coefs.gK.*n.*(V-coefs.EK) + ...
                    AA.*coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V) + ...
                    Stimvect(i)+Kindle(i)+CurrentExternal(i));%(coefs.E_EPSP-V))
            end
            
            dndt = (NInfLUT(Vind)-n)./NTauLUT(Vind);
            dssdt = -Ser/coefs.tauEPSPr;
            dsfdt = -Sef/coefs.tauEPSPf;
            dd1dt = (1-D1)/coefs.tauD1;
            dd2dt = (1-D2)/coefs.tauD2;
            dfdt = (1-F)/coefs.tauF;
            %             dd1dt = (1-D1)/DepTau(DepTauInd);
            %             dd2dt = (1-D2)/coefs.tauD2;
            %             dfdt = (1-F)/FacTau(FacTauInd);
            
            %lasty = y;
            %y = y + dydt*dt;
            lastV = V;
            V = V + dVdt*dt + (randn(N,1))*NoiseScalar;
            n = n + dndt*dt;
            Ser = Ser + dssdt*dt;
            Sef = Sef + dsfdt*dt;
            D1 = D1 + dd1dt*dt;
            D2 = D2 + dd2dt*dt;
            F = F + dfdt*dt;
            
            spks = find((V>thresh).*(lastV<thresh));
            if(~isempty(spks))
                ISI(spks)=t-lastSpkTime(spks);
                lastSpkTime(spks) = t;
                
%                 for spk=1:length(spks);
%                     if(iscell(Syns))
%                         Ser(Syns{spks(spk),:})=Ser(Syns{spks(spk),:})+D1(spks(spk)).*D2(spks(spk)).*F(spks(spk));
%                         Sef(Syns{spks(spk),:})=Sef(Syns{spks(spk),:})+D1(spks(spk)).*D2(spks(spk)).*F(spks(spk));
%                         D1(spks(spk))=D1(spks(spk))*d1;
%                         D2(spks(spk))=D2(spks(spk))*d2;
%                         F(spks(spk))=F(spks(spk))+f;
%                     else
%                         Ser(Syns(spks(spk),:))=Ser(Syns(spks(spk),:))+D1(spks(spk)).*D2(spks(spk)).*F(spks(spk));
%                         Sef(Syns(spks(spk),:))=Sef(Syns(spks(spk),:))+D1(spks(spk)).*D2(spks(spk)).*F(spks(spk));
%                         D1(spks(spk))=D1(spks(spk))*d1;
%                         D2(spks(spk))=D2(spks(spk))*d2;
%                         F(spks(spk))=F(spks(spk))+f;
%                     end;
%                 end;
                
                for spk=1:length(spks);
                    if(iscell(Syns))
                        Ser(Syns{spks(spk),:})=Ser(Syns{spks(spk),:})+1;
                        Sef(Syns{spks(spk),:})=Sef(Syns{spks(spk),:})+1;
                        %AA(spks(spk))=D1(spks(spk)).*D2(spks(spk)).*F(spks(spk));
                        AA(spks(spk))=D1(spks(spk));
                        D1(spks(spk))=D1(spks(spk))*d1;
                        D2(spks(spk))=D2(spks(spk))*d2;
                        F(spks(spk))=F(spks(spk))+f;
                    else
                Ser(Syns(spks(spk),:))=Ser(Syns(spks(spk),:))+1;
                Sef(Syns(spks(spk),:))=Sef(Syns(spks(spk),:))+1;
                %AA(spks(spk))=D1(spks(spk)).*D2(spks(spk)).*F(spks(spk));
                AA(spks(spk))=D1(spks(spk));
                D1(spks(spk))=D1(spks(spk))*d1;
                D2(spks(spk))=D2(spks(spk))*d2;
                F(spks(spk))=F(spks(spk))+f;
                    end
                end
                
                    spklist(spknum+1:spknum+1+length(spks),1)=t;
                    spklist(spknum+1:spknum+length(spks),2)=spks;
                    spknum=spknum+length(spks);
               
            end;
            Vm(i)=V(1); %put stim vector here
            %S(i)=(coefs.gEPSP.*(Sef(1)-Ser(1)).*(coefs.E_EPSP-V(1)));
            %Synch(i) = abs(mean(exp(j*2*pi*(t-lastSpkTime)/mean(ISI))));
            % S(i)=mean(D1.*D2.*F); %D1(1)*D2(1)*F(1);
            S(i)=mean(AA);
            
            if(i>=synchTime)
                synchTime=synchTime+SyncTimeStep;
                scount = scount + 1;
                Phase = mod(ceil((PhaseLUTLen-1)*(t-lastSpkTime)/mean(ISI)),PhaseLUTLen)+1;
                Synch(scount) = abs(mean(PhaseLUT(Phase)));
                h = hist(angle(PhaseLUT(Phase)),linspace(-pi,pi,100))/N;
                ind=find(h);
                entropy(scount)=-sum(h(ind).*log(h(ind)));
            end
            
            D(i)=length(spks);
            instISI(i)=mean(ISI);
            t=t+dt;
            
        end;
        %toc
        
       
            if spknum<length(spklist);
                spklist=spklist(1:spknum,:);
            end;
    
        %        Kuramoto(DepTauInd, FacTauInd)=mean(Synch(floor(scount*4/5):end));
        %
        %        Entropy(DepTauInd, FacTauInd)=mean(entropy(floor(scount*4/5):end));
        
        %         Synch=Synch(1:T/SyncTimeStep/dt);
        %         entropy=entropy(1:T/SyncTimeStep/dt);
        
        if(ThreeDon)
            Kuramoto(FreqInd,StrInd,CurrentInd)=mean(Synch(floor(scount*4/5):end));
            Firing(FreqInd,StrInd,CurrentInd)=mean(D(floor(scount*4/5):end));
            Entropy(FreqInd,StrInd,CurrentInd)=mean(entropy(floor(scount*4/5):end));
        else
            Kuramoto(FreqInd,StrInd)=mean(Synch(floor(scount*4/5):end));
            Firing(FreqInd,StrInd)=mean(D(floor(scount*4/5):end));
            Entropy(FreqInd,StrInd)=mean(entropy(floor(scount*4/5):end));
        end
        
        for i=1:length(Stimvect)
            if(Stimvect(i)==0)
                Stimvect(i)=NaN;
            end
        end
         
        if(Ploton)
            figure(667);
            plot(spklist(:,1),spklist(:,2),'.k')
            ylabel('Neuron');
            xlabel('Time (ms)');
           xlim([100 200]);
            ylim([0 1000]);
            
            figure(666);
            clf
            % subplot 311;
            % plot(Vm); axis tight;
            
            subplot 411
            %plot(spklist(:,1),spklist(:,2),'.')
            plot(dt:dt:T,D,'k',dt:dt:T,Stimvect,'.r');
            ylabel('Spikes');
            set(gca,'XTick',[0;T/4;T/2;T*3/4;T]);
            set(gca,'XTickLabel',[]);
            set(gca,'Position',[0.1 0.78 0.87 0.2]);           
            hold on
            
            subplot 412
            %plot(spklist(:,1),spklist(:,2),'.')
            plot(dt:dt:T,instISI,'k','LineWidth',2);
            ylabel('ISI');
             set(gca,'XTick',[0;T/4;T/2;T*3/4;T]);
            set(gca,'XTickLabel',[]);
            set(gca,'Position',[0.1 0.54 0.87 0.2]);
            hold on
            
            subplot 413
            plot(dt:dt:T,S,'k','LineWidth',2);
            ylabel('Syn Drive');
             set(gca,'XTick',[0;T/4;T/2;T*3/4;T]);
            set(gca,'XTickLabel',[]);
            set(gca,'Position',[0.1 0.31 0.87 0.2]);
            hold on
            
           subplot 414
            %plot(Synch2(:,1), Synch2(:,2),dt:dt:T, Synch);
            %plot(dt:dt:T, Synch);
           % plot(1:T/dt/SyncTimeStep,entropy(1:T/dt/SyncTimeStep),'k','LineWidth',2);
            plot(1:T/dt/SyncTimeStep,Synch(1:T/dt/SyncTimeStep),'k','LineWidth',2);

            ylabel('Synchrony');
            set(gca,'XTick',[0;T/dt/SyncTimeStep/4;T/dt/SyncTimeStep/2;T/dt/SyncTimeStep*3/4;T/dt/SyncTimeStep]);
             set(gca,'XTickLabel',[0;T/4;T/2;T*3/4;T]);
            set(gca,'Position',[0.1 0.08 0.87 0.2]);
            xlabel('Time (ms)');
            hold on 
            
%            subplot 515
%             plot(1:T/dt/SyncTimeStep,Synch(1:T/dt/SyncTimeStep)','k','LineWidth',2);
%             ylabel('Kuramoto');
%              xlabel('Time (ms)');
%             set(gca,'XTick',[0;T/dt/SyncTimeStep/4;T/dt/SyncTimeStep/2;T/dt/SyncTimeStep*3/4;T/dt/SyncTimeStep]);
%             set(gca,'XTickLabel',[0;T/4;T/2;T*3/4;T]);
%             set(gca,'Position',[0.1 0.09 0.87 0.14]);
        end
        
        if(Saveon)
            save DepressionGrid_Stim_Kindle_Burst_-2invest Firing Kuramoto Entropy ElecStimAmp StimPeriodVect instISI CurrentVect
        end
        
        toc
    end
    
    %     save DepressionGrid_Tau Kuramoto Entropy FacTau DepTau
    
    
end
% end



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




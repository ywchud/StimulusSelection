%% Get MU,PEAK,TAF,BLF for all StimC including single and multi

figure('WindowState', 'maximized')
bin = 0.005;
PSTHrange = (-0.05:bin:0.15);
% StimCs={StimC.W4;StimC.W2;StimC.W1;StimC.W3;StimC.W42;StimC.W21;StimC.W23;StimC.W421;StimC.W423;StimC.W213}; %StimC.W4;StimC.W2;
% titleStr={'W4','W2','W1','W3','W42','W21','W23','W421','W423','W213'}; %'W4','W2',
StimCs={StimC.W1;StimC.W2;StimC.W3;StimC.W4}; %StimC.W4;StimC.W2;
titleStr={'W1','W2','W3','W4'}; %'W4','W2',

timepoints=cell(1,length(StimCs));
for i=1:length(StimCs)
    if isempty(StimCs{i}.TouchDeltaT)
        timepoints{i}=StimCs{i}.TouchTimes;
%     else
%         [lateT, temp]=max([StimCs{i}.TouchTimes StimCs{i}.TouchTimes+StimCs{i}.TouchDeltaT],[],2);   %align with the touch that occurs later
%         timepoints{i}=lateT;%(temp==2);
    end
end


% timepoints={StimC2.TouchTimes;StimC3.TouchTimes;StimC23.TouchTimes};
MU=zeros(length(Neuron),length(StimCs));
MU_deBLF=zeros(length(Neuron),length(StimCs));
PEAKS=zeros(length(Neuron),length(StimCs));
TAF=zeros(length(Neuron),length(StimCs));
BLF=zeros(length(Neuron),1);
SpatialSelectivity=nan(length(Neuron),1);

for i=1:length(Neuron)
    subplot(Exp.si,Exp.sj,i)
    c=get(gca,'colororder');
    BLF(i)=Neuron(i).baselineFiring(Exp);
    for j=1:length(timepoints)
        [PSTHavg, SEM,mu,sd,peak,latency,PSTHmat]=Neuron(i).PSTH(timepoints{j},PSTHrange);
        MU(i,j)=mu;
        MU_deBLF(i,j)=mu-BLF(i);
        PEAKS(i,j)=peak;   
        TAF(i,j)=Neuron(i).TrialAvgFiring(Exp,StimCs{j}); 
    end
    
    plot(MU(i,:),'o')
    hold on
    plot(PEAKS(i,:),'r')
    plot(TAF(i,:),'b')
    addline('y',BLF(i),'m','-')
    addline('x',4.5,'g','-')
    addline('x',7.5,'g','-')
    % spatial selectivity
    pksOfStims=PEAKS(i,:)-BLF(i);
    SpatialSelectivity(i)=1-(norm(pksOfStims)/max(pksOfStims)-1)/(length(pksOfStims)^0.5-1); %range [0 1] for all stimC

    xlim([1 10])
    title(sprintf('C: %d, SS: %d',Neuron(i).CID,SpatialSelectivity(i)))
end
NN=nan(length(timepoints),1);
for i=1:length(timepoints)
    NN(i)=length(timepoints{i});
end

suptitle(compose('%s n = %d', string(titleStr(:)), NN(:)))
%check how correlated TAF and touch-induced MU are wityh pearson corr
R=nan(length(Neuron)-1,1);
figure
% for i=1:length(Neuron)
%     subplot(Exp.si,Exp.sj,i)
%     hold on
%     scatter(MU(i,:),TAF(i,:))
% 
%     r=corrcoef(MU(i,:),TAF(i,:));
%     R(i)=r(1,2);
%     lsline
%     xlabel('TouchAF(s_{-1})')
%     ylabel('TrialAF(s_{-1})')
%     title(sprintf('C: %d R=%0.5g',Neuron(i).CID,R(i)))
% end
for i=1:length(Neuron)-1
    hold on
    c=[0.8 0.8 0.8];
    scatter(MU(i,:),TAF(i,:),'MarkerFaceColor',c,'MarkerEdgeColor','none','MarkerEdgeAlpha',0.3)
    r=corrcoef(MU(i,:),TAF(i,:));
    R(i)=r(1,2);
end
hold on
scatter(MU(3,:),TAF(3,:),'MarkerFaceColor','r','MarkerEdgeColor','none')
scatter(MU(18,:),TAF(18,:),'MarkerFaceColor','b','MarkerEdgeColor','none')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

refline(1,0)
xlabel('TouchAF(s^{-1})')
ylabel('TrialAF(s^{-1})')
title('C14:red C79:blue')
%
figure 
histogram(R,'BinEdges',linspace(-1, 1, 20))
title('TchAF\_TrlAF Pearson corrcoef histogram')


%% Find if a neuron is responsive to a certain stimC
Responsivness=(TAF-repmat(BLF,1,length(StimCs)))./repmat(BLF,1,length(StimCs));
isResponsive=zeros(length(Neuron),length(StimCs));
isResponsive(Responsivness>0.5)=1;
isResponsive(Responsivness<-0.5)=-1;
%% Get MWI,Sharpness....for multiwhiskers only
%{'W4','W2','W1','W3','W42','W21','W23','W421','W423','W213'}
%MWI zero=purely additive; postive=superlinear; negative=
%sublinear/saturated
MStimC={StimC.W42;StimC.W21;StimC.W23;StimC.W421;StimC.W423;StimC.W213}; %StimC.W4;StimC.W2;
MWI=zeros(length(Neuron),length(MStimC));
MWI2=zeros(length(Neuron),length(MStimC));
MWI3=zeros(length(Neuron),length(MStimC));
MWI4=zeros(length(Neuron),length(MStimC));

SWall=zeros(length(Neuron),length(MStimC));
MWall=zeros(length(Neuron),length(MStimC));
SWall_bls=zeros(length(Neuron),length(MStimC));
MWall_bls=zeros(length(Neuron),length(MStimC));

Sharpness=zeros(length(Neuron),length(MStimC));
pkOnly=zeros(length(Neuron),length(MStimC));
win_short=zeros(length(Neuron),length(MStimC));
win_long=zeros(length(Neuron),length(MStimC));

MWSelectivity=nan(length(Neuron),1);


for i=1:length(Neuron)

    for j=1:length(MStimC)
        switch j
            case 1 %W42
                SCid=5;SSCid=[1 2];  %StimcCid;singleStimCid
            case 2 %W21
                SCid=6;SSCid=[2 3];
            case 3 %W23
                SCid=7;SSCid=[2 4];
            case 4 %W421
                SCid=8;SSCid=[1 2 3];
            case 5 %W423
                SCid=9;SSCid=[1 2 4];
            case 6 %W2123
                SCid=10;SSCid=[2 3 4];
            otherwise
                
        end
    
%         MWI(i,j)=log(PEAKS(i,SCid)/sum(PEAKS(i,SSCid)));
        MW=MU(i,SCid)-BLF(i);SW=sum(MU(i,SSCid))-BLF(i)*length(SSCid);
%         MW=(MU(i,SCid)-BLF(i))/BLF(i);SW=sum((MU(i,SSCid)-BLF(i))/BLF(i));
        MW_peak=PEAKS(i,SCid)-BLF(i);SW_peak=sum(PEAKS(i,SSCid))-BLF(i)*length(SSCid);
        
%         SWall(i,j)=sum(PEAKS(i,SSCid));
%         MWall(i,j)=PEAKS(i,SCid);
        SWall_bls(i,j)=SW;
        MWall_bls(i,j)=MW;
        
%         if ~isreal(log(MW/SW))
%             disp(log(MW/SW));
%         end
%         MW=MW_peak;SW=SW_peak;

        MWI(i,j)=(MW-SW);
        MWI2(i,j)=(MW-SW)/SW;
        MWI3(i,j)=(MW_peak-SW_peak)/SW_peak;
        MWI4(i,j)=(MW-SW)/(MW+SW);
        
        Sharpness(i,j)=PEAKS(i,SCid)/MU(i,SCid);  %not ideal, use width instead of mu
        pkOnly(i,j)=PEAKS(i,SCid)/BLF(i);
        win_short(i,j)=MU(i,SCid)/BLF(i);
        win_long(i,j)=TAF(i,SCid)/BLF(i);
  
    end
    MWII=MWI4(i,:);
    MWSelectivity(i)=1-(norm(MWII)/max(MWII)-1)/(length(MWII)^0.5-1); %wrong, MWII must be larger than 0
end
figure,histogram(MWSelectivity)
%plot
mymap = [1 0 0
    1 0.4 0.4
    1 0.7 0.7
    1 1 1
    0.7 0.7 1
    0.4 0.4 1
    0 0 1];
% 
% figure('WindowState', 'maximized')
% for i=1:length(Neuron)
%     subplot(Exp.si,Exp.sj,i)
%  	heatmap([MWI3(i,:)']);% Sharpness(i,:)' pkOnly(i,:)' win_short(i,:)' win_long(i,:)'])
%     colormap(mymap),caxis([-1 1])%caxis([-0.8 0.8])
%     title(sprintf('C: %d',Neuron(i).CID))
% end
figure('WindowState', 'maximized')
titleStr={'W42','W21','W23','W421','W423','W213'}; %'W4','W2',

h=heatmap(MWI4);% Sharpness(i,:)' pkOnly(i,:)' win_short(i,:)' win_long(i,:)'])
ax = gca;
ax.XData = titleStr;
colormap(mymap),caxis([-1 1])%caxis([-0.8 0.8])
[~,I]=sorty(h,{'W42','W421'},'descend');
figure('WindowState', 'maximized')
heatmap(Responsivness(str2double(I),5:end));
ax = gca;
ax.YData = I;
colormap(mymap),caxis([-2 2])
figure('WindowState', 'maximized')
for i=1:length(Neuron)
    scatter(MWI2(i,:),Responsivness(i,5:end))
    hold on
end

close all
%% 3d scatter of MU, TAF, BLF looking at single->pair->MW progression
figure
% StimCs={StimC.W4;StimC.W2;StimC.W1;StimC.W3;StimC.W42;StimC.W21;StimC.W23;StimC.W421;StimC.W423;StimC.W213};
% how MU,TAF..are arranged^^^
% MStimC={StimC.W4;StimC.W42;StimC.W421;StimC.W423;};
% 1 5 8 9 
X=PEAKS;
for i=1:length(Neuron)
    subplot(Exp.si,Exp.sj,i)
    x=[BLF(i) X(i,1)];
    y=[BLF(i) TAF(i,1)];
    hold on
    plot(x,y,'k')
    x=[X(i,1) X(i,5)];
    y=[TAF(i,1) TAF(i,5)];
    plot(x,y,'r')
    x=[X(i,5) X(i,8)];
    y=[TAF(i,5) TAF(i,8)];
    plot(x,y,'g')
    x=[X(i,5) X(i,9)];
    y=[TAF(i,5) TAF(i,9)];
    plot(x,y,'b')
    xlabel('PEAK')
ylabel('TAF')
end
%% 2d scatter, sums of SW vs MW mu 
figure
% StimCs={StimC.W4;StimC.W2;StimC.W1;StimC.W3;StimC.W42;StimC.W21;StimC.W23;StimC.W421;StimC.W423;StimC.W213};
% how MU,TAF..are arranged^^^
% MStimC={StimC.W4;StimC.W42;StimC.W421;StimC.W423;};
% 1 5 8 9 
X=MU;
for i=1:length(Neuron)
    subplot(Exp.si,Exp.sj,i)
    x=X(i,1)+X(i,2)-BLF(i)*2;  %W42
    y=X(i,5)-BLF(i);
    scatter(x,y,'k')
    hold on
    x=X(i,2)+X(i,3)-BLF(i)*2;  %W21
    y=X(i,6)-BLF(i);
    scatter(x,y,'r')
    x=X(i,2)+X(i,4)-BLF(i)*2; %W23
    y=X(i,7)-BLF(i);
    scatter(x,y,'b')
    x=X(i,1)+X(i,2)+X(i,3)-BLF(i)*3;  %W421
    y=X(i,8)-BLF(i);
    scatter(x,y,'g')
    x=X(i,1)+X(i,2)+X(i,4)-BLF(i)*3;  %W423
    y=X(i,9)-BLF(i);
    scatter(x,y,'m')
    x=X(i,2)+X(i,3)+X(i,4)-BLF(i)*3; %W213
    y=X(i,10)-BLF(i);
    scatter(x,y,'c')
    refline(1,0)
    xlabel('sum SWs')
    ylabel('MW')
end
suptitle('K-W42 R-W21 b-w23 g-421 m-W423 c=W213')
%% PSTH of chosen neurons

figure('WindowState', 'maximized')
bin = 0.005;
% PSTHrange = (-0.05:bin:0.15);
PSTHrange = (-0.01:bin:0.06);
figWin=[-0.005 0.035];
StimCN=3;
NeuronIdx=[31]; %3,18,31
for k=1:StimCN
    subplot(3,StimCN,k)
switch k
    case 1
        StimCs={StimC.W4;StimC.W2;StimC.W1;StimC.W421}; %StimC.W4;StimC.W2;
        titleStr={'W4','W2','W1','W421'}; 
    case 2
        StimCs={StimC.W4;StimC.W2;StimC.W42}; %StimC.W4;StimC.W2;
        titleStr={'W4','W2','W42'};
%         StimCs={StimC.W4;StimC.W2;StimC.W3;StimC.W423}; %StimC.W4;StimC.W2;
%         titleStr={'W4','W2','W3','W423'};
    case 3
%         StimCs={StimC.W2;StimC.W1;StimC.W3;StimC.W213}; %StimC.W4;StimC.W2;
%         titleStr={'W2','W1','W3','W213'};
        StimCs={StimC.W2;StimC.W1;StimC.W3;StimC.W213}; %StimC.W4;StimC.W2;
        titleStr={'W2','W1','W3','W213'};
    otherwise
end
% NeuronIdx=1:length(Neuron);


timepoints=cell(1,length(StimCs));
for i=1:length(StimCs)
    if isempty(StimCs{i}.TouchDeltaT)
        timepoints{i}=StimCs{i}.TouchTimes;
    else
        [lateT, temp]=max([StimCs{i}.TouchTimes StimCs{i}.TouchTimes+StimCs{i}.TouchDeltaT],[],2);   %align with the touch that occurs later
        timepoints{i}=lateT;%(temp==2);
    end
end


% timepoints={StimC2.TouchTimes;StimC3.TouchTimes;StimC23.TouchTimes};


for i=1:length(NeuronIdx)
%     subplot(Exp.si,Exp.sj,i)
    

    for j=1:length(timepoints)
        [PSTHavg, SEM,mu,sd,peak,latency,PSTHmat]=Neuron(NeuronIdx(i)).PSTH(timepoints{j},PSTHrange);
%         plot(PSTHrange(2:end),PSTHavg);
%         c=get(gca,'colororder');
        if j==length(timepoints)   %last stimC which is the MW
            c=[1 0 0];
        else %SW
            c=[0 0 0.1*j];
        end
        plotWFilledError(PSTHrange(2:end)-bin/2,PSTHavg,c(1,:),PSTHavg-SEM,PSTHavg+SEM,min([[1 1 1];c(1,:)*1.3]),0.5)
        hold on      
        addline('y',Neuron(NeuronIdx(i)).TrialAvgFiring(Exp,StimCs{j}),c(1,:),'--')  
    end
    addline('y',Neuron(NeuronIdx(i)).baselineFiring(Exp),'k','-')
    xlim(figWin)
    NN=nan(length(timepoints),1);
    for m=1:length(timepoints)
        NN(m)=length(timepoints{m});
    end
    if k==1,ylabel('Spikes (s^{-1})'),end
%     title(strcat(sprintf('C: %d ',Neuron(NeuronIdx(i)).CID),string(compose(' %s n = %d', string(titleStr(:)), NN(:)))))
    
    str1=string(compose(' %s n = %d', string(titleStr(:)), NN(:)));
    title(str1)
    %     title(sprintf('C: %d',Neuron(NeuronIdx(i)).CID))

    subplot(3,StimCN,k+StimCN)
%     for j=1:length(timepoints)
%         if j==length(timepoints)   %last stimC which is the MW
%             c=[1 0 0];
%         else %SW
%             c=[0 0 0.1*j];
% 
%         end
%         Neuron(NeuronIdx(i)).raster(timepoints{j},PSTHrange,c(1,:));
%     end
    c=[1 0 0];
    EmptyTouchId=Neuron(NeuronIdx(i)).raster(timepoints{j},PSTHrange,c(1,:));
%     CheckTouch(Exp,Whisker,StimCs{j},EmptyTouchId)

    xlim(figWin)
    if k==1,ylabel('Trials'),end
    subplot(3,StimCN,k+StimCN*2)
    MUs=nan(length(timepoints)-1,length(PSTHrange)-1);
    for j=1:length(timepoints)-1
        [PSTHavg, SEM,mu,sd,peak,latency,PSTHmat]=Neuron(NeuronIdx(i)).PSTH(timepoints{j},PSTHrange);
        MUs(j,:)=PSTHavg;
    end   
    SWsum=sum(MUs,1)-BLF(NeuronIdx(i))*size(MUs,1);
    c=[0.1 0.1 0.1];
    plot(PSTHrange(2:end)-bin/2,SWsum,'Color',c);
    hold on
    [PSTHavg, SEM,mu,sd,peak,latency,PSTHmat]=Neuron(NeuronIdx(i)).PSTH(timepoints{end},PSTHrange);
    c=[1 0 0];
    plot(PSTHrange(2:end)-bin/2,PSTHavg-BLF(NeuronIdx(i)),'Color',c)
    xlim(figWin)
    if k==1,ylabel('Baseline substracted Spikes (s^{-1})'),end
    legend({'SW SUM','MW'})
end
end
suptitle(sprintf('C: %d ',Neuron(NeuronIdx(i)).CID))
suplabel('Time from touch (s)');

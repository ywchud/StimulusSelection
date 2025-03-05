function [res,newStim]=findPSTH(StimCs,Exp,Whisker,Neuron,NeuronIdx,p)
%This function serves as a hub for recruiting subfunctions which determine
%aligntime based on given p,calcuating PSTH based on this aligntime,
%restructuring outputs into a single structure res

%Input StimCs is finalised and not further filtered. Use updateStimC
%beforehand as needed

PSTHrange=p.PSTHrange;
windowtime=p.windowtime;



PSTHs=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
PSTHmat=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
PSTHrelT=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));

firstSpkLats=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
SEMs=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));

    MU=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    PREZERO=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    SD=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    N=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    PEAKS=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    TROUGHS=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));

    LAT=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    ABS_PEAKS=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    sign=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    TAF=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    BLF=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    
    fwhm=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    dhm=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    x1=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    x2=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    
    firstSpkLat=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    firstSpkJitter=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    firstSpkTendency=zeros(length(NeuronIdx),size(StimCs,2),size(StimCs,1));

    countpertrial=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    ratepertrial=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));
    CorrTrial=cell(length(NeuronIdx),size(StimCs,2),size(StimCs,1));

    
for k=1:size(StimCs,1)
%     fprintf('Working on row %d...\n',k)
    % added newStim 08/13/2022 due to ActiveTrial filter in getTimePoints.
    % downstream could be buggy if k>1 or NeuronIdxN >1
    [timepoints,CT,newStim]=getTimePoints(StimCs(k,:),Exp,Whisker,Neuron,NeuronIdx,p);  %Important: consider filtering with updateStimC first. Only Touch PW is filtered here
  
    for i=1:length(NeuronIdx)
        for j=1:size(newStim,2)
            S=newStim{1,j};
            if size(timepoints,1)>1
                tp=timepoints{i,j};
            else
                tp=timepoints{j};
            end

            [PSTHavg,SEM,SumStat]=Neuron(NeuronIdx(i)).PSTH(tp,PSTHrange,windowtime);
            PSTHs{i,j,k}=PSTHavg;
            PSTHmat{i,j,k}=SumStat.PSTHmat;
            PSTHrelT{i,j,k}=SumStat.PSTHrelT;
            SEMs{i,j,k}=SEM;
            taf=Neuron(NeuronIdx(i)).TrialAvgFiring(Exp,S);

                MU(i,j,k)=SumStat.mu;
                PREZERO(i,j,k)=SumStat.preZero;
                SD(i,j,k)=SumStat.sd;
                N(i,j,k)=length(tp);
                PEAKS(i,j,k)=SumStat.peak;
                TROUGHS(i,j,k)=SumStat.trough;
                ABS_PEAKS(i,j,k)=SumStat.abs_peak;
                LAT(i,j,k)=SumStat.latency;
                sign(i,j,k)=SumStat.sign;
                TAF(i,j,k)=taf;
                

                fwhm(i,j,k)=SumStat.fwhm;
                dhm(i,j,k)=SumStat.deltahm;
                x1(i,j,k)=SumStat.x1;
                x2(i,j,k)=SumStat.x2;
                
                firstSpkLat(i,j,k)=SumStat.firstSpkLat;
                firstSpkJitter(i,j,k)=SumStat.firstSpkJitter;
                firstSpkTendency(i,j,k)=SumStat.firstSpkTendency;
                firstSpkLats{i,j,k}=SumStat.firstSpkLats;
                
                if S.Light==1
                    BLF(i,j,k)=Neuron(NeuronIdx(i)).baselineFiringLightOn(Exp);
                else
                    BLF(i,j,k)=Neuron(NeuronIdx(i)).baselineFiring(Exp);
                end

                bin=PSTHrange(2)-PSTHrange(1);
                countpertrial{i,j,k}=sum(SumStat.PSTHmat(:,SumStat.zeroIndex:SumStat.lastIndex),2);
                ratepertrial{i,j,k}=mean(SumStat.PSTHmat(:,SumStat.zeroIndex:SumStat.lastIndex),2)/bin;
                
                
                try
                    CorrTrial{i,j,k}=CT{i,j};
                catch
                    CorrTrial{i,j,k}=CT{1,j};
                end
%                 CorrTrial(:,:,k)=CT;
        end
                
    end
end

res.PSTHs=PSTHs;   %vector 1*nbin
res.PSTHmat=PSTHmat; %mat trN*nbin
res.PSTHrelT=PSTHrelT; %cell trN*1
res.SEMs=SEMs;  %singular
res.MU=MU;   %singular
res.PREZERO=PREZERO;
res.SD=SD;   %singular
res.N=N;   %singular
res.PEAKS=PEAKS;   %singular
res.TROUGHS=TROUGHS; 
res.ABS_PEAKS=ABS_PEAKS; 



res.LAT=LAT;   %singular
res.sign=sign;
res.TAF=TAF;   %singular
res.BLF=BLF;   %singular
res.fwhm=fwhm;   %singular
res.dhm=dhm;   %singular
res.x1=x1;   %singular
res.x2=x2;   %singular

res.firstSpkLats=firstSpkLats;
res.firstSpkLat=firstSpkLat;
res.firstSpkJitter=firstSpkJitter;
res.firstSpkTendency=firstSpkTendency;

res.countpertrial=countpertrial;
res.ratepertrial=ratepertrial;   %vector trN*1
res.NeuronIdx=toColumn(NeuronIdx);   %singular
res.CorrTrial=CorrTrial; %vector trN*1
end

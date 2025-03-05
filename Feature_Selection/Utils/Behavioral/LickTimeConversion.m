LickTimes=cell(Exp.TrN,1);

temp=RawAnalogLicksCC1;array=[1 0 1 0];
II=find(PistonComb(array,Exp.Stim.Piston.Mat));
for i=1:size(temp,1)
    tempp=temp(i,:);
    tempp=tempp(tempp~=0);
    LickTimes{II(i)}=tempp;
end

temp=RawAnalogLicksDD1;array=[0 1 0 1];
II=find(PistonComb(array,Exp.Stim.Piston.Mat));
for i=1:size(temp,1)
    tempp=temp(i,:);
    tempp=tempp(tempp~=0);
    LickTimes{II(i)}=tempp;
end

temp=RawAnalogLicksCD1;array=[1 0 0 1];
II=find(PistonComb(array,Exp.Stim.Piston.Mat));
for i=1:size(temp,1)
    tempp=temp(i,:);
    tempp=tempp(tempp~=0);
    LickTimes{II(i)}=tempp;
end

temp=RawAnalogLicksDC1;array=[0 1 1 0];
II=find(PistonComb(array,Exp.Stim.Piston.Mat));
for i=1:size(temp,1)
    tempp=temp(i,:);
    tempp=tempp(tempp~=0);
    LickTimes{II(i)}=tempp;
end

%determine HITS
LickHitsBuffer=0.3;


HITS=nan(Exp.TrN,1);
for i=1:length(LickTimes)
    LT=LickTimes{i};
    
    if ~isempty(LT)
        firstLick=LT(1);
        if firstLick>Exp.Stim.Piston.Delay+LickHitsBuffer && firstLick<Exp.WaterDelay   
            %first lick occurs during valid hit criteria: i.e. between PistonOut+Arduino buffer(300ms) and water delivery time
            HITS(i)=1;
        elseif sum(LT>=Exp.Stim.Piston.Delay+LickHitsBuffer & LT<Exp.WaterDelay)>=1
            %first lick occurs before hit criteria but there is at least
            %one subsequent hit during valid hit criteria, so still
            %registered as hits
            HITS(i)=1;   % can set as 2 if u need to segregate them
        else
            %no lick during valid hit criteria
            HITS(i)=0;
        end
    else
        HITS(i)=0; 
    end

end

%% Plot Lick and touch raster
close all
figure,TrialTouchRaster(Exp,Whisker,[StimC_CC;StimC_DD],LickTimes,HITS,[],[1 0]),title('Control Hits')
figure,TrialTouchRaster(Exp,Whisker,[StimC_CC_0;StimC_DD_0],LickTimes,HITS,[],[1 0]),title('Light Hits')
figure,TrialTouchRaster(Exp,Whisker,[StimC_CC;StimC_DD;StimC_CC_0;StimC_DD_0],LickTimes,HITS,[],[1 0]),title('All trials Hits')

figure,TrialTouchRaster(Exp,Whisker,[StimC_CC;StimC_DD],LickTimes,~HITS,[],1),title('Homo Misses')
figure,TrialTouchRaster(Exp,Whisker,[StimC_CD;StimC_DC],LickTimes,HITS,[],1),title('Hetero False Positives')
figure,TrialTouchRaster(Exp,Whisker,[StimC_CD;StimC_DC],LickTimes,~HITS,[],1),title('Hetero Correct Rejection')

%% Plot Touch Number before first lick box plots
close all
hits=HITS;
% hits(111:end)=0;  %only consider 1:110
% hits(1:110)=0;  %only consider 111:end

[TNb4L,~,~]=TrialTouchRaster(Exp,Whisker,StimC_CC,LickTimes,hits,[],0);
W1_Ctrl=TNb4L(:,1);W3_Ctrl=TNb4L(:,3); 
[TNb4L,~,~]=TrialTouchRaster(Exp,Whisker,StimC_CC_0,LickTimes,hits,[],0);
W1_Opt=TNb4L(:,1);W3_Opt=TNb4L(:,3);
[TNb4L,~,~]=TrialTouchRaster(Exp,Whisker,StimC_DD,LickTimes,hits,[],0);
W2_Ctrl=TNb4L(:,2);W4_Ctrl=TNb4L(:,4);
[TNb4L,~,~]=TrialTouchRaster(Exp,Whisker,StimC_DD_0,LickTimes,hits,[],0);
W2_Opt=TNb4L(:,2);W4_Opt=TNb4L(:,4);

origin=[ones(length(W3_Ctrl),1)*1;ones(length(W3_Opt),1)*2;ones(length(W4_Ctrl),1)*3;ones(length(W4_Opt),1)*4;ones(length(W1_Ctrl),1)*5;ones(length(W1_Opt),1)*6;ones(length(W2_Ctrl),1)*7;ones(length(W2_Opt),1)*8];
TNb4L=[W3_Ctrl;W3_Opt;W4_Ctrl;W4_Opt;W1_Ctrl;W1_Opt;W2_Ctrl;W2_Opt];


figure
boxplot(TNb4L,origin)
xticklabels({'CL_Ctrl','CL_Opt','DL_Ctrl','DL_Opt','CR_Ctrl','CR_Opt','DR_Ctrl','DR_Opt'})
ylabel('Touch Number before lick')

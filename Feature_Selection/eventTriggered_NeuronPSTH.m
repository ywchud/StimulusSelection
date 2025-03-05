function Neuron=eventTriggered_NeuronPSTH(Exp,Neuron)
% get Post Stimulus Time Histogram (PSTH) of all neuron spikes from user-defined parameters
% spike count matrix stored in Neuron.R and the parameters in Neuron.PARAMS
% read NeuronPSTH.m for detailed param options
tic
for ii=1:10
    fprintf('PARAM %d\n',ii)
clear fp p
%note: any non free whisking StimC(i.e. as long as there is sth in S.WID) uses S.updateTouches to filter p.XXX, meaning only trials with touches will retain.)
%So only trials with touches will be filtered regardless of
%aligntype(p.type) 
%e.g. p.isrun does not affect p.type='PistonStart' 

p.firstfew=[];
p.ActiveTrialsOnly = 1;   %for each neuron, only trials with blf > ActiveTh is used
p.Sequence = [];  %if positive, only >0 deltaT is used; if negative, only <0 deltaT is used, if 0 , only =0 deltaT is used, if empty, ignore
p.windowtime=[0 0.03];
p.ITIth=0.03;
p.isrun=[];
p.Tdur=[];%[0.008 0.015];
p.Toverlap=[];
p.Tgap=[];
p.mainW='PW';
p.omitTrial=[];

switch ii
    case 1
        p.bin = 0.005;
        p.PSTHrange = (-0.1:p.bin:0.06);
        p.isrun=1;
        p.type='Touch';%p.type='Protraction_Tonly';
        temp=1:Exp.TrN;
%         temp(~(logical(~Exp.Lickport.location) & ~Exp.WhiskerTrimTr & (Exp.Lick.Hit | Exp.Lick.CR)))
        p.omitTrial=temp(logical(Exp.Lickport.location) | Exp.WhiskerTrimTr | Exp.Lick.Miss | Exp.Lick.FA | Exp.Lick.earlyLick | Exp.Lick.weirdTrials);  %home hit/CR (omit leave,W trimmed,miss,FA)        
    case 2
        p.bin = 0.005;
        p.PSTHrange = (-0.1:p.bin:0.06);
        p.isrun=1;
        p.type='Touch';%p.type='Protraction_Tonly';
        temp=1:Exp.TrN;
        p.omitTrial=temp(logical(Exp.Lickport.location) | Exp.WhiskerTrimTr | Exp.Lick.Hit | Exp.Lick.CR | Exp.Lick.earlyLick | Exp.Lick.weirdTrials);  %home miss/FA (omit leave,W trimmed,HIT,CR)      
    case 3
        p.bin = 0.005;
        p.PSTHrange = (-0.1:p.bin:0.06);
        p.isrun=1;
        p.type='Touch';%p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;
        p.omitTrial=temp(~Exp.Lickport.location | Exp.WhiskerTrimTr);   %leave (omit home,W trimmed)
    case 4
        p.bin = 0.005;
        p.PSTHrange = (-0.1:p.bin:0.06);
        p.isrun=1;
        p.type='Touch';%p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;  
        p.omitTrial=temp(~ismember(temp',Exp.Lickport.home_trials) | Exp.WhiskerTrimTr);  % home_begin (omit non home_begin trials)
    case 5
        p.bin = 0.005;
        p.PSTHrange = (-0.1:p.bin:0.06);
        p.isrun=1;
        p.type='Touch';%p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;  
        p.omitTrial=temp(~ismember(temp',Exp.Lickport.leave_trials) | Exp.WhiskerTrimTr);  % leave_begin (omit non home_begin trials)  
        
        

    case 6
        p.bin = 0.025;%1/Exp.Fs;
        p.PSTHrange = (-3:p.bin:6);
        p.windowtime=[0 1];
        p.type='PistonStart';    %check getTimePoints.m for a list of options
        temp=1:Exp.TrN;
        p.omitTrial=temp(logical(Exp.Lickport.location) | Exp.WhiskerTrimTr | Exp.Lick.Miss | Exp.Lick.FA | Exp.Lick.earlyLick | Exp.Lick.weirdTrials);  %home (omit leave)
    case 7
        p.bin = 0.025;%1/Exp.Fs;
        p.PSTHrange = (-3:p.bin:6);
        p.windowtime=[0 1];
        p.type='PistonStart'; %p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;
        p.omitTrial=temp(logical(Exp.Lickport.location) | Exp.WhiskerTrimTr | Exp.Lick.Hit | Exp.Lick.CR | Exp.Lick.earlyLick | Exp.Lick.weirdTrials);  %home miss/FA (omit leave,W trimmed,HIT,CR)
    case 8
        p.bin = 0.025;%1/Exp.Fs;
        p.PSTHrange = (-3:p.bin:6);
        p.windowtime=[0 1];
        p.type='PistonStart'; %p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;        
        p.omitTrial=temp(~Exp.Lickport.location | Exp.WhiskerTrimTr);   %leave (omit home)
    case 9
        p.bin = 0.025;%1/Exp.Fs;
        p.PSTHrange = (-3:p.bin:6);
        p.windowtime=[0 1];
        p.type='PistonStart'; %p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;
        p.omitTrial=temp(~ismember(temp',Exp.Lickport.home_trials) | Exp.WhiskerTrimTr);  % home_begin (omit non home_begin trials)
    case 10
        p.bin = 0.025;%1/Exp.Fs;
        p.PSTHrange = (-3:p.bin:6);
        p.windowtime=[0 1];
        p.type='PistonStart'; %p.type='Protraction_Tonly'; 
        temp=1:Exp.TrN;
        p.omitTrial=temp(~ismember(temp',Exp.Lickport.leave_trials) | Exp.WhiskerTrimTr);  % leave_begin (omit non home_begin trials)  

end
StimCi=Exp.StimCi_1D;
titleStr=Exp.titleStr_1D;

p.titleStr=titleStr;
fp(1)=p;  %posttouch
% fp(2)=p;
% fp(2).firstfew=[2];
% fp(3)=p;
% fp(3).firstfew=[3];
PARAMS(ii,1)=p;
NeuronIdx=1:length(Neuron);%[3 18];
% holdplotN=max([length(fp) 1])*max([length(StimCi) 1]);
% c=createColorMap(holdplotN,[0 0 0;1 0 0;0 1 0;0 0 1]);

clear Res res ResBLF resBLF
%%%%%%% 
res=NeuronPSTH(Exp,Whisker,Neuron,NeuronIdx,StimCi,fp);
%%%%%%%NeuronPSTH(Exp,Whisker,Neuron([Neuron.CID]==2),1,StimCi(8, 1),fp);
res_custom(:,:,1,ii)=res(:,1:Exp.MStim.N,:);
end
toc
% disp('Remember to run the next subsection to save into R')
%%  Store RES_custom into individual neurons as R
temp=size(res_custom);
for i=1:length(NeuronIdx)
%         Neuron(NeuronIdx(i)).R=squeeze(res_custom(NeuronIdx(i),:,:,:));
        Neuron(NeuronIdx(i)).R=reshape(res_custom(NeuronIdx(i),:,:,:),temp(2:end)); %remove first dim (neuron)
        Neuron(NeuronIdx(i)).PARAMS=PARAMS;
end

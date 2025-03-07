%% Extract Behavioral features such as whisking kinematics and task performance
%% Clean wipe
% clear, clc, close all
%% Load save file
temp=dir(fullfile('.','Saves\ExampleExperiment','Preprocess_backup_*'));
temp=temp(1);
load(fullfile(temp.folder,temp.name))
%% Extract raw whisker coordinates from .csv files
[Exp,Xmat,Ymat]=ExtractDLC(Exp,1, 30, 4); %(Exp,Fc1, Fc2, order)
%in case any whiskers are not labelled, resize Xmat Ymat
[Xmat,Ymat]=resizeDLC(Exp,Xmat,Ymat);
%% Create whisker objects
Whisker=WhiskerClass.empty;
for i=1:Exp.Stim.Piston.tot
    if Exp.Stim.Piston.active(i)
        Whisker(i) = WhiskerClass(Exp,i,Xmat,Ymat);
    else
        Whisker(i) = WhiskerClass(Exp,0,Xmat,Ymat);
    end
end
%% Frame matching and cleaning with digital camera signals
Exp=checkTrials(Exp,Whisker(Exp.Stim.Piston.active),Sigs);
%if each video has different videofps, recreate whisker.trial objects
if length(Exp.videoFps)~=1
    Whisker=WhiskerClass.empty;
    for i=1:Exp.Stim.Piston.tot
        if Exp.Stim.Piston.active(i)
            Whisker(i) = WhiskerClass(Exp,i,Xmat,Ymat);
        else
            Whisker(i) = WhiskerClass(Exp,0,Xmat,Ymat);
        end
    end
end
%% Find precise individual piston out and in frame(accounts for varying videoFps too, used in findtouchtimeV2)
% might take a couple of minutes
Exp=MechanicalDelayAdjustment_Auto(Exp,[],[],Ymat,Xmat);
checkPistonFrames(Exp) %visually check piston frames
%% Determine ROI where a whisker is considered 'touching' a piston
Exp.Stim.Piston.ROIfraction=4;
Exp.Stim.Piston.ROI=[];
for i=1:Exp.Stim.Piston.tot
    fprintf('Finding ROI for Piston %d...\n',Exp.Stim.Piston.ID(i))
    [Exp,~]=findPistonROI(Exp,i,[0.5 100],1);
end
Exp.Stim.Piston.ROIradius=[17 17 17 17];
%% Find Runspeed per trial with fs=Exp.videoFps
Exp.RunSpeed=findRunSpeed(Exp,Sigs.runspeed);
Exp.Runcmpers=findRunSpeed(Exp,Sigs.runspeedcmpers);
%% Find whisker kinematics:
W=Whisker(Exp.Stim.Piston.active);
W=findWhiskerbase(Exp,W); %find whiskerbase
W=findWhiskerParam(Exp,W);  %find angle, phase, curvature
[W,Exp]=CurvDebase(Exp,W); %find CurvDebase: curvature with baseline removed, may not work for retracted angles with fewer sample
W=findWhiskerEnvelope(Exp,W,[],[],[],1); %Exp,Whisker,th,absMin,maxDur,samesize %find envelope related signals: pks,trs,setpoint,locsP,locsT,Amplitude
Whisker(Exp.Stim.Piston.active)=W;
%% play video to check whisker kinematics
temp=find(PistonComb([0 1 0 0],Exp.Stim.Piston.Mat)); %find the trial to look at
playVid(Exp,[1],temp(randi(length(temp))),[],Whisker(Exp.Stim.Piston.active),"CenAngle",[500 1500],1)  %CenAngle,Curvature,Phase,CurvDebase,DeltaPhase
%% Find Touchtime
W=Whisker(Exp.Stim.Piston.active);
W=findTouchTimeV2(Exp,W,0);   
Whisker(Exp.Stim.Piston.active)=W;
%% Check Touch
WhiskerID=1;tr=1;
CheckWhiskerTouch(Exp,Whisker(WhiskerID),tr,[0.2 1000],0);
%% Further fine-tuning as needed
% %----RefineROI trial by trial as needed and rerun above
% Exp.Stim.Piston.ROIradius=[25 25 12 12];
% W=Whisker(Exp.Stim.Piston.active);
% Exp=refineROI(Exp,W);
% W=findTouchTimeV2(Exp,W,0);   
% Whisker(Exp.Stim.Piston.active)=W;
% %----rough touch determination with cenangle due to poor tracking
% Exp.Stim.Piston.ROIradius=[12 12 12 12];
% W=Whisker(Exp.Stim.Piston.active);
% clear p
% p.poorTrackingContingency=[1 0];
% p.th=0.8;  %threshold(rad) to be considered as touch, determined maunally by checking touch from below
% W=findTouchTimeV2(Exp,W,0,p);   
% Whisker(Exp.Stim.Piston.active)=W;
%% Create StimCondition objects that outline the combination of whiskers and stimuli involved
% e.g. SC task go whisker & no go whisker only
Exp.MStim.N=2;
Exp.MStim.type=1:2;
Exp.MStim.recipe=[1 0;
                  0 1;
                  0 0;
                  0 0];
Exp.MStim.theta=nan(1,Exp.MStim.N);%[1/3 nan -1/3 0 1/3 -1/3 0 1/6 -1/6 nan];
Exp.MStim.SingleN=Exp.Stim.Piston.num;%4
Exp.MStim.MultiN=0;
Exp.MStim.Free=0;
Exp.MStim.Str={'C1','D1'};

ITIth=[];%[]/Exp.videoFps;  %normally 13/Exp.videoFps
FrameRange=[0 max(Exp.FrameN)];
StimC.W1 = StimCondition(Exp,Whisker,"C1",[1],0,[],[],ITIth,[],FrameRange);
StimC.W2 = StimCondition(Exp,Whisker,"D1",[2],0,[],[],ITIth,[],FrameRange);

Exp.StimCi={StimC.W1,StimC.W2};
Exp.titleStr=cellfun(@(x) x.type,Exp.StimCi);
Exp.StimCi_1D=reshape(Exp.StimCi',length(Exp.StimCi(:)),1);
Exp.titleStr_1D=reshape(Exp.titleStr',length(Exp.StimCi(:)),1);
cellfun(@(x) length(x.TouchTimes),Exp.StimCi_1D)

%e.g. C2 centred multiwhisker with optogenetics on
% Exp.MStim.N=11;
% Exp.MStim.type=1:11;
% Exp.MStim.recipe=[1 0 0 0 1 0 0 1 1 0 0;
%                   0 1 0 0 0 1 0 1 0 1 0;
%                   0 0 1 0 0 0 1 0 1 1 0;
%                   0 0 0 1 1 1 1 1 1 1 0];
% Exp.MStim.MultiN=6;
% Exp.MStim.Free=1;
% Exp.MStim.Str={'B1','C1','D1','C2','B1C2','C1C2','D1C2','B1C1C2','B1D1C2','C1D1C2','Free'};
% 
% StimC.W43_O = StimCondition(Exp,Whisker,"C2D_O",[4 3],1,0,[],ITIth,[]);
%% Check StimC averaged touch frame
figure
[si,sj]=subplotDim(length(Exp.StimCi(1,:)));
for i=1:2
    subplot(2,1,i)
    extractframe_avg(Exp,Exp.StimCi{1,i})
end
%     CheckTouch(Exp,Whisker,StimC.W2,[]) %video
%% Classify Trials into hit/cr/miss/fa
Exp.Lick.goW=2;
Exp.Lick.nogoW=1;
Exp=TrialPerformace_discrimination(Exp);
%% Specify Cutoff trial where whiskers are trimmed
Exp.nthTrialCutoff=Exp.TrN;
Exp.WhiskerTrimTr=zeros(Exp.TrN,1);
Exp.WhiskerTrimTr(Exp.nthTrialCutoff:Exp.TrN)=1;
%%
%% Neuron objects
% Get spikes
removeTh=700;%700;
[Exp,Neuron]= getSpikesInfo(Exp,removeTh); %create neuron class, get spiketimes for each neuron
%% Neural stability across trials
[Exp,Neuron]=NeuralActivity(Exp,Neuron);
%% Neuron spike shape
% Determine RS/FS/UN 
clear p
p.dataDir=Exp.Path.data;
p.fileName='amplifier.dat'; %always the same
p.dataType='int16'; %always the same
p.nCh = sum(Exp.KSort.ChN); % Number of total channels that were streamed to disk in .dat file (2*128 ch probe)
p.wfWin = [-40 41];% Number of samples before and after spiketime to include in waveform
p.nWf = 1200;   % Number of waveforms per unit to pull out
p.spikeTimes =readNPY(fullfile(p.dataDir,'\spike_times.npy'));% Vector of cluster spike times (in samples) same length as .spikeClusters
p.spikeClusters = readNPY(fullfile(p.dataDir,'\spike_clusters.npy'));%Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes% loop through each neuron, pseudorandomly take nWf numbers of spike times from each neuron, and get wave forms; 
p.ClusterInfo = tdfread(fullfile(p.dataDir,'\cluster_info.tsv'),'\t');

wf= getWaveForms(p,[Neuron.CID]);
[Neuron,~,~]=getFSRS(Exp,Neuron,Exp.KSort.chanMap-1,wf,p.ClusterInfo); %tdfread(fullfile(p.dataDir,'\cluster_info.tsv'),'\t')
save(fullfile(Exp.Path.save,'wf.mat'),'wf')
%% Neuron PSTHS
% define parameters as needed for each condition in function 
% (e.g. alignType(Touch or PistonStart..), PSTH range and bin size, trial
% used, etc)
Neuron=eventTriggered_NeuronPSTH(Exp,Neuron);
%% Save
Exp.savenote='Behavioral_features_updated';
variableList = {'Exp','Whisker'};
sfilename=['FeatSelect_backup_', datestr(floor(now)) '.mat'];
filename=fullfile(Exp.Path.save,sfilename);
disp('Saving...')
for i = 1:length(variableList)
    disp(i) 
    if exist(variableList{i},'var')         
         if exist(filename,'file')
            save(filename,variableList{i},'-append')
         else
             save(filename,variableList{i})
         end
        
    end
end
whos('-file',filename)
disp('Saved')
%Neurons only
variableList = {'Neuron'};%,'RES','res_custom','PARAMS','HITS','LickTimes','dprime','StimC'};
sfilename=['Neu_backup_', datestr(floor(now)) '.mat'];
filename=fullfile(Exp.Path.save,sfilename);
disp('Saving...')
for i = 1:length(variableList)
    disp(i) 
    if exist(variableList{i},'var')         
         if exist(filename,'file')
            save(filename,variableList{i},'-append')
         else
             save(filename,variableList{i})
         end
        
    end
end
whos('-file',filename)
disp('Saved')

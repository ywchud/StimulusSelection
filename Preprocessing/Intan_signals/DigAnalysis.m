function [Exp,Sigs]=DigAnalysis(Exp,Sigs,isplot)
%% Trial signal processing
t=Exp.t;
disp('Processing Trial Signal...')
%added 04/27/21, sometimes there are flickers in trialstart signal
Sigs.Trial.sig=removeFlicker(Sigs.Trial.sig,1);
%get the trial start and end indices and duration in seconds
TrialStart = Sigs.Trial.onset;
TrialEnd = Sigs.Trial.offset;
% TrialStart = find(diff(DigInput(:,TrialDig))>0);
% TrialEnd = find(diff(DigInput(:,TrialDig))<0);
first_removed=0;
last_removed=0;
    if Sigs.Trial.sig(1)==1
        TrialStart=TrialStart(2:end);
        TrialEnd=TrialEnd(2:end);
        disp('First Trial signal removed due to an abrupt start during camera capture')
        first_removed=1;
    end
    if Sigs.Trial.sig(end)==1
        TrialStart=TrialStart(1:end-1);
        TrialEnd=TrialEnd(1:end-1);
        disp('Last Trial signal removed due to an abrupt end during camera capture')
        last_removed=1;
    end

%     
%trial 190 void
redundantN=0;
%Exp.TrN=no. of mp4s ; length(TrialStart)=no. of intan
if  Exp.TrN==0
    redundantN=0;
    fprintf('No mp4s fpund. Assuming %d trials from Diginput now',length(TrialStart))
elseif Exp.TrN<length(TrialStart)
    redundantN=length(TrialStart)-Exp.TrN;
    fprintf('Number of .mp4(%d) and Trial Signal(%d) mismatch\nRemoving the last %d trial signal\n',Exp.TrN,length(TrialStart),redundantN)
%     Exp.TrN=length(TrialStart);
%     fprintf('Trial Number set to (%d)\n',Exp.TrN)
elseif Exp.TrN>length(TrialStart)
    fprintf('More .mp4(%d) found than signal(%d) can find. Trial signal may be incomplete or corrupted, or unrelated videos found in folder.\n',Exp.TrN,length(TrialStart))
    if last_removed
        fprintf('Last mp4 is likely cut short, removed from analysis but please check the video yourself and confirm\n')
    elseif first_removed
        error('lol how, did your first video started recording midway')
    else
        error('No idea what happened plz check')
    end
end

Exp.TrialStart=TrialStart(1:end-redundantN);
Exp.TrialEnd=TrialEnd(1:end-redundantN);
Exp.TrialStartT=Exp.t(Exp.TrialStart);
Exp.TrialEndT=Exp.t(Exp.TrialEnd);
Exp.TrialDur=Exp.TrialEndT-Exp.TrialStartT;
Exp.TrN=length(Exp.TrialStart);

if isplot
    figure
    hist(Exp.TrialDur,0:0.5:20)
    title('Trial Duration')
end
disp('Obtained TrN, TrialStart, TrialEnd, TrialStartT, TrialEndT and TrialDur')
%% Wheel and camera processing
% 
%running speed
disp('Processing Wheel Signal...')
WheelTics = Sigs.Wheel.offset;
if any(WheelTics>length(t))
    disp('Incomplete signal, probably ran out of memory during recording. Cutting losses...')
    WheelTics=WheelTics(WheelTics<=length(t));
end
Fs=1000;
[WheelSpeed, WheelSpeedbin] = hist(t(WheelTics),0:1/Fs:max(t));
WheelSpeedSmooth = smooth(WheelSpeed,200,'lowess')*1000;

RunSpeedNorm=WheelSpeedSmooth/prctile(WheelSpeedSmooth,90);
if isplot
    figure
    plot(WheelSpeedbin,RunSpeedNorm)
    title('Running speed')
end
Sigs.runspeed=Signal([],[],RunSpeedNorm,Fs);
runtime=WheelSpeedbin(RunSpeedNorm>0.2);
Exp.runStartT=toColumn(runtime(diff([-100 runtime])>1));
Exp.runEndT=toColumn(runtime(diff([runtime runtime(end)+100])>1));
Sigs.runtime=Signal([],[],runtime,Exp.Fs);

%get cm/s instead
if Exp.Intan.WheelPort==15  %old rig
    
    TicksPerRevolution=5000;%400;
    WheelRadius=8.5; %cm
    cmPerSec=2*pi*WheelRadius/TicksPerRevolution*WheelSpeedSmooth;
    Sigs.runspeedcmpers=Signal([],[],cmPerSec,Fs);
    
elseif Exp.Intan.WheelPort==14 %new rig
    TicksPerRevolution=400;
    WheelRadius=9; %cm
    cmPerSec=2*pi*WheelRadius/TicksPerRevolution*WheelSpeedSmooth;
    Sigs.runspeedcmpers=Signal([],[],cmPerSec,Fs);
end
% runtime=WheelSpeedbin(WheelSpeedSmooth/max(WheelSpeedSmooth)>0.2);

% Exp.runStartT=Sigs.runtime.onset;
% Exp.runEndT=Sigs.runtime.offset;

disp('Obtained runtime signal,runStartT and runEndT')
%camera delay from trialstart
disp('Processing Camera Signal...')
sig=Sigs.Camera.sig;
FrameOnsets=Sigs.Camera.onset;


N2Avg=Exp.TrN;
CamDelay=nan(N2Avg,1);
CamFps=nan(N2Avg,1);
CamOffWin=nan(N2Avg,1);
trial=floor(1:Exp.TrN/N2Avg:Exp.TrN);
for i=1:N2Avg  %sample 10 trials for camDelay
    FrameOnseti=FrameOnsets(FrameOnsets>=Exp.TrialStart(trial(i)) & FrameOnsets<=Exp.TrialEnd(trial(i)));
    try
    CamDelay(i)=Exp.t(FrameOnseti(1)-Exp.TrialStart(trial(i)));
    catch
        disp(i)
    end
    CamOffWin(i)=Exp.t(Exp.TrialEnd(trial(i))-FrameOnseti(end));
    CamDelay(i)=find(sig(TrialStart(trial(i)):TrialEnd(trial(i)))==1,1,'first')/Exp.Fs;
    CamFps(i)=Exp.Fs/mode(diff(FrameOnseti));
%      
end
Exp.Cam.Delay=mean(CamDelay);
Exp.videoFps=mean(CamFps);
Exp.Cam.OffWin=mean(CamOffWin);
Exp.Cam.Delays=CamDelay;
Exp.Cam.OffWins=CamOffWin;

if sum(abs(CamDelay-Exp.Cam.Delay)>1/Exp.Fs)~=0
    disp('Camera Delay is not consistent among trials, but minor issue')
end
if sum(abs(CamOffWin-Exp.Cam.OffWin)>1/Exp.Fs)~=0
    disp('Camera Off Window is not consistent among trials')
end
fprintf('Obtained CamDelay %d, videoFps %d\n',Exp.Cam.Delay,Exp.videoFps)
% (length(DigInput(TrialStart(trial):TrialEnd(trial),8))-max(find(DigInput(TrialStart(trial):TrialEnd(trial),8)==1)))/timestampFs
%% Lick info processing
if isfield(Sigs,'Spout')
disp('Processing Spout Signal...')

%%
if Sigs.Spout.Port==11  %old rig where spout signal is linear motor signal
start_move = zeros(length(Exp.TrialStart),1);
location = zeros(length(Exp.TrialStart),1);
ticks=[diff(Sigs.Spout.sig);nan];
temp=[1;Exp.TrialEnd];
for i = 1:length(Exp.TrialStart)
%     tick_times=find(ticks(Exp.TrialStart(i):Exp.TrialEnd(i))~=0); % What times did the drum ticks happen
try
    tick_times=find(ticks(temp(i):temp(i+1))~=0); % extended window to cover between n-1 and n TrialEnd instead since spout move outside of trial for later exp
catch
123    
end
%     tick_times = tick_times + Exp.TrialStart(i);
    no_ticks(i) = length(tick_times);
    if length(tick_times)>5
        start_move(i) = tick_times(1)/Exp.Fs;
    end
end
% cluster ticks into two groups (home & leave)
idx = kmeans([no_ticks(no_ticks>5)]',2);
nt = no_ticks(no_ticks>5)';
% lowest ticks is home
home_thres = min([min(nt(idx==1)) min(nt(idx==2))]);
% highest ticks is leave
leave_thres = max([min(nt(idx==1)) min(nt(idx==2))]);
leave_trials = find(no_ticks>=leave_thres);
home_trials = find(no_ticks>=home_thres & no_ticks<leave_thres);
if length(home_trials)-length(leave_trials)==1 && home_trials(1)==1
    home_trials=home_trials(2:end);
end
% for i = 1:min([length(leave_trials),length(home_trials)])
%     location(leave_trials(i):home_trials(i)-1) = 1;
% end

home=zeros(1,length(home_trials));
leave=ones(1,length(leave_trials));
spoutTriggerTr=[home_trials leave_trials];
[spoutTriggerTr,I]=sort(spoutTriggerTr);
spoutTrigger=[home leave];
spoutTrigger=spoutTrigger(I);
pos=0;
for m=1:Exp.TrN
    if ~ismember(m,spoutTriggerTr)
        location(m)=pos;
        continue
    end
    pos=spoutTrigger(m==spoutTriggerTr);
    location(m)=pos;
end

% In the event that exp ends in a leave block
% if length(leave_trials)>length(home_trials)
%     location(leave_trials(end):end) = 1;
% end
% save variables
Exp.Lickport.location = location;
Exp.Lickport.leave_trials = leave_trials;
Exp.Lickport.home_trials = home_trials;
Exp.Lickport.move_time = start_move;
elseif Sigs.Spout.Port==1  %bilateral rig where spout signal is intan input signal
    for t = 1:length(Exp.TrialStartT)
        lick_port_times = Exp.TrialStart(t)+0.5*Exp.Fs:Exp.TrialEnd(t); % Piston shoots of .65/.85 from start of trial
        Exp.Lickport.location(t,1) = all(Sigs.Spout.sig(lick_port_times)); % DI 9,10,11,12 are piston
        try
            Exp.Lickport.move_time(t,1) = find(diff(Sigs.Spout.sig(Exp.TrialStart(t):Exp.TrialEnd(t)))~=0,1,'first');
        catch
            Exp.Lickport.move_time(t,1) = 0;
        end
    end
%     Exp.Lickport.location_cont = Exp.DigInput(:,1);
    temp=diff([Exp.Lickport.location;nan]);
    Exp.Lickport.leave_trials = find(temp==1)+1;
    Exp.Lickport.home_trials = [1;find(temp==-1)+1];
end
else
    disp('Spout inactive')
end
%%
if isfield(Sigs,'LickPiezo')
disp('Processing LickPiezo Signal...')

s=Sigs.LickPiezo.sig;
v_th=((max(s)-min(s))/6)+min(s); %2
s_binary = Signal( Sigs.LickPiezo.Port,Sigs.LickPiezo.type, logical(s>v_th), Sigs.LickPiezo.fs );
Exp.Lick.time=s_binary.onset/Exp.Fs;
Exp.Lick.duration=s_binary.offset/Exp.Fs-Exp.Lick.time;
Exp.Lick.Frames=TimeToFrame(Exp.Lick.time,Exp.TrialStartT,Exp.TrialEndT,Exp.videoFps,zeros(Exp.TrN,1),Exp.Cam.Delay);
disp('Licktimes are not final, check Exp.Lick.Frames to connect consecutive frames as 1 lick rather than many')

%%
Exp.Lick.Hit=false(Exp.TrN,1);
Exp.Lick.Miss=false(Exp.TrN,1);
Exp.Lick.FA=false(Exp.TrN,1);
Exp.Lick.CR=false(Exp.TrN,1);
for t=1:Exp.TrN
trial_times_ext = Exp.TrialStart(t):Exp.TrialEnd(t)+Exp.Fs;
    if isfield(Sigs,'MissIntan')
        s=Sigs.MissIntan.sig;
        if isempty(s)
            fprintf('Sigs.MissIntan.sig is empty, probably forgot to record. Exp.Lick.Miss set to all zeros for now.\n')
            break;
        end
        if nnz(s(trial_times_ext))~= 0
            Exp.Lick.Miss(t,1) = 1;
        end
        
    elseif isfield(Sigs,'HitIntan')
        s=Sigs.HitIntan.sig;
        if isempty(s)
            fprintf('Sigs.HitIntan.sig is empty, probably forgot to record. Exp.Lick.Hit set to all zeros for now.\n')
            break;
        end
        if sum(s(trial_times_ext))~= 0
            Exp.Lick.Hit(t,1) = 1;
        end
    end

    s=Sigs.FAIntan.sig;
    if nnz(s(trial_times_ext))~= 0
        Exp.Lick.FA(t,1) = 1;
    end
end
end
%%

if isfield(Sigs,'Solenoid')
disp('Processing Solenoid Signal...')

Solenoid_on=Sigs.Solenoid.onset;
Solenoid_off=Sigs.Solenoid.offset;

if length(Solenoid_on)~=length(Solenoid_off)
    error('Unequal solenoid onset offset length')
end

Solenoid_open=zeros(Exp.TrN,1);
Solenoid_open_time=nan(Exp.TrN,1);
for i=1:Exp.TrN
    temp=find(Solenoid_on>=Exp.TrialStart(i) & Solenoid_off<=Exp.TrialEnd(i));
    if length(temp)==1
        Solenoid_open(i)=1;
        Solenoid_open_time(i)=Solenoid_on(temp)-Exp.TrialStart(i);
    elseif length(temp)>1
        fprintf('Trial %d: Solenoid opened for >1 times, picking only the first one\n',i)
        Solenoid_open(i)=1;
        Solenoid_open_time(i)=Solenoid_on(temp(1))-Exp.TrialStart(i);
    end
end
if length(Solenoid_on)~=sum(Solenoid_open)
    fprintf('Solenoid open outside of some recorded trials, should not matter if you only look within trials though\n')
end
% save variables
Exp.Solenoid.open=Solenoid_open;
Exp.Solenoid.open_time=Solenoid_open_time/Exp.Fs;
elseif isfield(Sigs,'LickPiezo') && ~isempty(Sigs.LickPiezo.sig)  %no solenoid signal but there is lick behavior(bilateral rig), manually create Exp.Solenoid
    

    %-----------estimate solenoid open time from peak lick rate------------
%     Solenoid_open=zeros(Exp.TrN,1);
%     Solenoid_open_time=nan(Exp.TrN,1);
%     RelTimes=cell(Exp.TrN,1);
%     for i=1:Exp.TrN
%         Solenoid_open(i)=Exp.Lick.Hit(i);%sum(Exp.Lick.time>=Exp.TrialStartT(i)+0.5  & Exp.Lick.time<Exp.TrialEndT(i))>2;
%         if Solenoid_open(i)
%             RelTimes{i}=Exp.Lick.time(Exp.Lick.time>=Exp.TrialStartT(i)+0.5  & Exp.Lick.time<Exp.TrialEndT(i))-Exp.TrialStartT(i);
%         end
%     end
%     try
%         RTs=cell2mat(RelTimes);
%     catch
%         RTs=cell2mat(RelTimes');
%     end
%     RTs=RTs(:);
%     [N,edges] = histcounts(RTs,50);
%     [~,I]=max(N);
%     opentime=edges(I);
%     for i=1:Exp.TrN    
%         if Solenoid_open(i)
%             Solenoid_open_time(i)=opentime;
%         end
%     end
%     Exp.Solenoid.open=Solenoid_open;
%     Exp.Solenoid.open_time=Solenoid_open_time;
    %-------------------------------
    %---------estimate solenoid open time using onset of hit signal
    if 1
        disp('Estimating solenoid open time using onset of hit signal')
%         Sigs.HitIntan.sig
        Solenoid_on=Sigs.HitIntan.onset;
        Solenoid_off=Sigs.HitIntan.offset;
        
        if length(Solenoid_on)~=length(Solenoid_off)
            error('Unequal solenoid onset offset length')
        end
        
        Solenoid_open=zeros(Exp.TrN,1);
        Solenoid_open_time=nan(Exp.TrN,1);
        for i=1:Exp.TrN
            temp=find(Solenoid_on>=Exp.TrialStart(i) & Solenoid_off<=Exp.TrialEnd(i));
            if length(temp)==1
                Solenoid_open(i)=1;
                Solenoid_open_time(i)=Solenoid_on(temp)-Exp.TrialStart(i);
            elseif length(temp)>1
                fprintf('Trial %d: Solenoid opened for >1 times, picking only the first one\n',i)
                Solenoid_open(i)=1;
                Solenoid_open_time(i)=Solenoid_on(temp(1))-Exp.TrialStart(i);
            end
        end
        if length(Solenoid_on)~=sum(Solenoid_open)
            fprintf('Solenoid open outside of some recorded trials, should not matter if you only look within trials though\n')
        end
        % save variables
        Exp.Solenoid.open=Solenoid_open;
        Exp.Solenoid.open_time=Solenoid_open_time/Exp.Fs;
    end
    
end
disp('Exp.Lick.Hits/Miss/FA/CR are not final, run --Get hit and CR-- in mastercode ~line 1107')

%% plot Sigs
if isplot
    plotIntanSig(Exp,Sigs)
end
% ylim([0 1.8])
%%
% fc=500;
% [b,a] = butter(6,fc/(timestampFs/2));
% 
% for i = 1:num_channels % downsample to 1k
%     FP(i,:) = filtfilt(b,a,FieldPotential(i,:)) ;
% end
% t_FP=resample(t,1,30);
% 
% figure
% for i = 1:num_channels 
%     subplot(4,8,i)
%     plot(FP(i,1000:1000+200));
% end
end
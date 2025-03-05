function Exp=getOpticsInfo(Exp,Sigs)
disp('Getting Optics info...')
    OpticsStart=nan(Exp.TrN,1);
    OpticsEnd=nan(Exp.TrN,1);
    OpticsDelay=nan(Exp.TrN,1);
    OpticsStop=nan(Exp.TrN,1);
    freq=nan(Exp.TrN,1);
    amp=nan(Exp.TrN,1);
    if Exp.Stim.Optic.num==1
        v=Sigs.Optic.sig;
        OptNorm=(v-min(v))/max(v-min(v));
    end
    
    if exist('OptNorm','var') && isempty(OptNorm)
        error('You dummy did you forget to record the analog signal again')
    end
        
    UP=cell(Exp.TrN,1);
    DOWN=cell(Exp.TrN,1);
    for i=1:Exp.TrN
        if Exp.Stim.Optic.num==1
            temp=OptNorm(Exp.TrialStart(i):Exp.TrialEnd(i));  
        else
            temp=nan;
        end
        amp(i)=round(mean(temp(temp>0.2)),2);
        temp2=find(temp>0.2,1,'first');
        temp3=find(temp>0.2,1,'last');
        transition=[]; %to find Hz
        up=1;down=1;temp0=temp;
        while ~isempty(up) && ~isempty(down) 
            up=find(temp0>0.2,1,'first')+down-1;
            temp0=temp(up:end);
            down=find(temp0<0.2,1,'first')+up-1;
            temp0=temp(down:end);
            transition=[transition;up;down];
        end
        try
        UP{i}=round((transition(1:2:end)/Exp.Fs-Exp.Cam.Delay)*Exp.videoFps);  %onset
        DOWN{i}=round((transition(2:2:end)/Exp.Fs-Exp.Cam.Delay)*Exp.videoFps);  %offset
        freq(i)=round(Exp.Fs/mean(diff(transition)*2));
       
        
        catch
            disp('fu')
        end
        
        if ~isempty(temp2)
            OpticsDelay(i)=temp2;
            OpticsStop(i)=temp3;
        end
         if ~isnan(OpticsDelay(i))      %there is optic
            OpticsStart(i)=Exp.TrialStart(i)+OpticsDelay(i);  
            OpticsEnd(i)=Exp.TrialStart(i)+OpticsStop(i);
         end
         fprintf('Completed Trial: %d\n',i)
    end
    OpticsStartTime=OpticsStart/Exp.Fs;
    OpticsEndTime=OpticsEnd/Exp.Fs;
    % OpticsDelay=mean(OpticsStartTime-TrialStartTime,'omitnan');
    
    Exp.Stim.Optic.Start=OpticsStart;
    Exp.Stim.Optic.End=OpticsEnd;
    Exp.Stim.Optic.StartT=OpticsStartTime;
    Exp.Stim.Optic.EndT=OpticsEndTime;
    Exp.Stim.Optic.Delay=mean(OpticsDelay,'omitnan')/Exp.Fs;
    Exp.Stim.Optic.OffWin=mean(Exp.TrialEnd-Exp.TrialStart-OpticsStop,'omitnan')/Exp.Fs;
    
    Exp.Stim.Optic.On=~isnan(OpticsStart);
    Exp.Stim.Optic.f=freq;
    amp(isnan(amp))=0;
    Exp.Stim.Optic.amp=amp;
    Exp.Stim.Optic.OnsetF=UP;
    Exp.Stim.Optic.OffsetF=DOWN;
    disp('Got optics start/end startT/endT Delay and On')
    freqq=unique(freq(~isnan(freq)));
    str="";
    for i=1:length(freqq)
        str=strcat(str,num2str(freqq(i)),' ');
    end
    str=strcat('Frequencies found at ',str, ' Hz');
    disp(str)
    [a,b]=findgroups(Exp.Stim.Optic.amp);
    str="";
    for i=1:max(a)
        str=strcat(str,num2str(b(i)),' ');
    end
    str=strcat('Amplitudes found at ',str, ' Hz');
    disp(str)
end
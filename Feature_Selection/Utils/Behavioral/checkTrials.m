function Exp=checkTrials(Exp,Whisker,Sigs)
     Badtrials=[];
     FrameOnsets=Sigs.Camera.onset;
     FrameDiff=zeros(Exp.TrN,1);
     ObservedFrameN=zeros(Exp.TrN,1);
     ExpectedFrameN=zeros(Exp.TrN,1);
     firstSample=zeros(Exp.TrN,1);
     lastSample=zeros(Exp.TrN,1);
     if ~isfield(Exp,'OldFps')
         Exp.OldFps=Exp.videoFps;
     end
     
     for i=1:Exp.TrN
         Samples=(FrameOnsets>Exp.TrialStart(i) & FrameOnsets<Exp.TrialEnd(i));
         firstSample(i)=FrameOnsets(find(Samples,1,'first'));
         lastSample(i)=FrameOnsets(find(Samples,1,'last'));
         ExpectedFrameN(i)=sum(Samples);
%          (Exp.TrialEndT(i)-Exp.TrialStartT(i)-Exp.Cam.Delay)*Exp.videoFps;
%            %inaccurate, gives 1-3 frame error likely due to inconsistent
%            trialEnds from last frame or irregular camera signal(unlikely
%            but who knows)
         ObservedFrameN(i)=Whisker(1).trial(i).FrameN;
         Exp.FrameN(i,1)=ObservedFrameN(i);
         FrameDiff(i,1)=ObservedFrameN(i)-ExpectedFrameN(i);
         
         if abs(FrameDiff(i,1))>=1
             Badtrials=[Badtrials;i];
             fprintf('Bad trial with uneven frameN found: %d; Missing: %d\n',i,FrameDiff(i,1))
             Exp.FrameT{i,1}(:,1)=nan;
         else
             Exp.FrameT{i,1}(:,1)=Exp.t(FrameOnsets(Samples));
         end
     end
     Exp.FrameDiff=FrameDiff;

     
     %if all trials are bad, there is likely an inconsistency between input
     %Exp.videoFps versus actual recording rate due to hardware issues
     if length(Badtrials)/Exp.TrN>0.05  %many bad >5%
         fprintf('%0.2f trials are bad\n',length(Badtrials)/Exp.TrN)
         if mean(FrameDiff<=0)>=0.8  %mostly negative
             disp('Likely an oversized frame size so camera cannot keepup with input Exp.videoFps')
             disp('Reassigning Exp.videoFps... assuming constant')
             newFps=zeros(Exp.TrN,1);
             OldFps=Exp.OldFps;
             for i=1:Exp.TrN
                 vidtime=sum(FrameOnsets>Exp.TrialStart(i) & FrameOnsets<Exp.TrialEnd(i))/OldFps;
                 newFps(i)=Whisker(1).trial(i).FrameN/vidtime;
             end
%              NewFps=round(mean(newFps));
             NewFps=newFps;

%              sampleNperCycle=Exp.Fs/NewFps;
             fprintf(['New Exp.videoFps set from %d to:' repmat('%0.2f\n',1,Exp.TrN)],OldFps,NewFps);
             Exp.videoFps=NewFps;
             Exp.FrameT=[];
             for i=1:Exp.TrN
                 Exp.FrameN(i,1)=Whisker(1).trial(i).FrameN;
%                  temp=firstSample(i):sampleNperCycle:firstSample(i)+sampleNperCycle*Exp.FrameN(i,1);
%                  Exp.FrameT{i,1}(:,1)=Exp.t(round(temp));
                 usable_size=(lastSample(i)-firstSample(i))+1;   
                 temp=firstSample(i):usable_size/ObservedFrameN(i):lastSample(i);
                 Exp.FrameT{i,1}(:,1)=Exp.t(round(temp));
             end
             Badtrials=[];
         else
             error('Unknown error')
         end
     end
     
     if isfield(Exp,'BadTrials')
         Exp.BadTrials=unique([Exp.BadTrials;Badtrials]);
     else
         Exp.BadTrials=unique(Badtrials);
     end
    disp('Screening completed. Bad trials are stored in Exp.BadTrials. FrameT and FrameN obtained for each trial.')
    

%         Exp.FrameT{i,1}(:,1)=linspace((Exp.TrialStartT(i)+Exp.Cam.Delay),Exp.TrialEndT(i),Exp.FrameN(i,1));

end
function Whisker=findTouchTimeV2(Exp,Whisker,filterSAmp,p)
%(PLOI,TLOI,Tcoor,istouch,maxTouchNo,WhiskerAngle,Setpoint,Amplitude,WhiskerLabels,videoFps,PistonDelay,WhiskerError,TrialStart,t,CamDelay,ROIradius)
% ROIradius=10;

%filterSAmp=0.01 ~[0 1] quantile of small amplitude filtered out
  if exist('p','var') && isfield(p,'ESP')
      ESP=p.ESP;
  else
      ESP=Exp.Stim.Piston;
  end
  if exist('p','var') && isfield(p,'poorTrackingContingency')
      poorTrackingContingency=p.poorTrackingContingency;
      th=p.th;
      percentOverth=0.05;
  else
      poorTrackingContingency=zeros(1,length(Whisker));
      th=nan;
      percentOverth=nan;
  end
  
  
for i=1:length(Whisker)
    W=Whisker(i);
    Wid=W.W_id;
    if exist('p','var') && isfield(p,'onPiston')
        onPiston=p.onPiston; %whisker wid touching different piston onpiston
    else
        onPiston=Wid;  %whisker wid touching same piston wid
    end
    
    WTLOI=W.TLOI;
    if ESP.TLOI(Wid)~=WTLOI
        disp('Detected a different assigned TLOI, replacing Whisker.TLOI with ESP.TLOI')
        Whisker(i).TLOI=ESP.TLOI(Wid);
    end
    
    try
        garbo=ESP.ROI(onPiston).x2;
        ROIpt=2;
    catch
        ROIpt=1;
    end
    
for j = 1:Exp.TrN
    WT=W.trial(j);
    
        LT=WT.LocsT;  %protraction: location of troughs
        LP=WT.LocsP;  %retraction: location of peaks
        if length(LP)~=length(LT)
            error('LocsT and LocsP have different sizes, please rerun findWhiskerEnvelope with TrsPksSamesize=1');
        end
        Pnum=length(LP);

    
    LargeWhisk=nan(Pnum,1); %(F1)
    AfterPiston=nan(Pnum,1);%(F2)
    EnterROI=nan(Pnum,1);%(F3)
    NotManyOutsideROI=nan(Pnum,1);%(F4)
    EnoughInROI =nan(Pnum,1);%(F5)
    NoSlip=nan(Pnum,1);%(F6)    
        
        
    %must have large amplitude (F1)
    if filterSAmp
        Amplitude=WT.Amplitude(LT);
        LargeWhisk=reshape(Amplitude>quantile(Amplitude,filterSAmp),[Pnum,1]);
    else
        LargeWhisk=true(Pnum,1);
    end
%     PistonBufferT=Exp.PistonBuffer(j,i)/Exp.videoFps; % Mechanical delay of piston shoot-out time following intan signal onset (empirically verified) 0.1s start +0.1s fully out
    tolerance=5; %tolerate n frames to account for rare cases when piston arrives after whisker protraction but still touches
    PistonBufferF=Exp.PistonBuffer(j,onPiston)-tolerance;
    %must be after piston shoots out (F2)
    if ~isnan(ESP.Mat(j,onPiston*2))  %only if a piston is present
%         AfterPiston=reshape(LT>Exp.videoFps*(ESP.Delay-Exp.Cam.Delay+PistonBufferT)  & LT<Exp.FrameN(j)-Exp.videoFps*(ESP.OffWin-Exp.Cam.OffWin),[Pnum,1]);
        if length(Exp.videoFps)==1
            fps=Exp.videoFps;
        else
            fps=Exp.videoFps(j);
        end
        if isfield(Exp,'PistonBufferReturn')
            AfterPiston=reshape(LT>=PistonBufferF  & LT<Exp.PistonBufferReturn(j,onPiston),[Pnum,1]);
        else
            error('Obsolete code, use MechanicalDelayAdjustment_Auto to get PistonReturn Frame and rerun')
%             AfterPiston=reshape(LT>=PistonBufferF  & LT<Exp.FrameN(j)-fps*(ESP.OffWin-Exp.Cam.OffWins(j)),[Pnum,1]);
        end
    else
        AfterPiston=false(Pnum,1);
    end

        P_Frame = LT;
        R_Frame = LP;
        
        % Define T_frame as moment of touch and T_End as moment of touch release
        T_Frame=nan(Pnum,1);
        T_End=nan(Pnum,1);
        %% Find ROI
        if ESP.ROI(onPiston).isdrift
            id=find(ESP.ROI(onPiston).index_trial==j);
            if isempty(id)
                ROIx1=0;
                ROIy1=0;
            else
                ROIx1=ESP.ROI(onPiston).x1_intp(id,end);  %use last column, if run refineROI, column 2 will be filled
                ROIy1=ESP.ROI(onPiston).y1_intp(id,end);
            end
        else
            ROIx1=ESP.ROI(onPiston).x1;
            ROIy1=ESP.ROI(onPiston).y1;
        end
        if ROIpt>1
            if ESP.ROI(onPiston).isdrift
                id=find(ESP.ROI(onPiston).index_trial==j);
                if isempty(id)
                    ROIx2=0;
                    ROIy2=0;
                else
                    ROIx2=ESP.ROI(onPiston).x2_intp(id,end);
                    ROIy2=ESP.ROI(onPiston).y2_intp(id,end);
                end
            else
                ROIx2=ESP.ROI(onPiston).x2;
                ROIy2=ESP.ROI(onPiston).y2;
            end
        end
%%
            GapDistance=nan(Pnum,1);
            GapFrame=nan(Pnum,1);
            % filter based on touch label within region of interest
            for k=1:Pnum       
%                 buffer=min([20 length(WhiskerLabels{j}(:,(TLOI)*2))-R_Frame(k)]); %20 extra frames or whatever is left needed for label to stay in ROI
%                 y2=WhiskerLabels{j}(P_Frame(k):R_Frame(k)+buffer,(TLOI)*2);   %touch LOI
%                 x2=WhiskerLabels{j}(P_Frame(k):R_Frame(k)+buffer,(TLOI)*2-1);
               
%changed to fixed buffer frames after protraction since retraction point is
%unreliable when touching
%find touch from window protraction to 20 frames after
                TLOI=rem(WTLOI,Exp.DLC.LblPerW);
                if length(ESP.ROIradius)>1
                    ROIradius=ESP.ROIradius(onPiston);
                else
                    ROIradius=ESP.ROIradius;
                end
                try
                    Pwin=P_Frame(k+1)-P_Frame(k);  %protraction window from current protraction to the next protraction
                catch
                    Pwin=WT.FrameN-P_Frame(k); %protraction window from current protraction to last frame of vid
                end
 %-----------------------               
% added 05/14/21 can now use interpolated TLOI
if rem(TLOI,1)==0  %TLOI is a positive integer/double without .XXX
                y2=WT.Y(TLOI).sig(P_Frame(k):P_Frame(k)+Pwin-1);   %touch LOI
                x2=WT.X(TLOI).sig(P_Frame(k):P_Frame(k)+Pwin-1);
else            %TLOI is an interpolated point
                y2_i=WT.Y(floor(TLOI)).sig(P_Frame(k):P_Frame(k)+Pwin-1);
                y2_f=WT.Y(ceil(TLOI)).sig(P_Frame(k):P_Frame(k)+Pwin-1);
                x2_i=WT.X(floor(TLOI)).sig(P_Frame(k):P_Frame(k)+Pwin-1);
                x2_f=WT.X(ceil(TLOI)).sig(P_Frame(k):P_Frame(k)+Pwin-1);

                y2=((y2_f-y2_i)*rem(TLOI,1))+y2_i;
                x2=((x2_f-x2_i)*rem(TLOI,1))+x2_i;
end
%---------use a different label for slip detection since TLOI sometimes get
%'stuck on piston' by DLC while other labels go thru
n=0.6; %move to the next n label on whisker
if rem(TLOI+n,1)==0  %TLOI is a positive integer/double without .XXX
                y2_shift=WT.Y(TLOI+n).sig(P_Frame(k):P_Frame(k)+Pwin-1);   %touch LOI
                x2_shift=WT.X(TLOI+n).sig(P_Frame(k):P_Frame(k)+Pwin-1);
else            %TLOI is an interpolated point
                y2_i=WT.Y(floor(TLOI+n)).sig(P_Frame(k):P_Frame(k)+Pwin-1);
                y2_f=WT.Y(ceil(TLOI+n)).sig(P_Frame(k):P_Frame(k)+Pwin-1);
                x2_i=WT.X(floor(TLOI+n)).sig(P_Frame(k):P_Frame(k)+Pwin-1);
                x2_f=WT.X(ceil(TLOI+n)).sig(P_Frame(k):P_Frame(k)+Pwin-1);

                y2_shift=((y2_f-y2_i)*rem(TLOI+n,1))+y2_i;
                x2_shift=((x2_f-x2_i)*rem(TLOI+n,1))+x2_i;
end
 %-----------              
                
                try                  
                    gapD=nan(Pwin,1);
                    for m=1:Pwin
                        gapD(m)=pDistance(x2(m),y2(m),ROIx1,ROIy1,ROIx2,ROIy2);
                    end
                    inROI=gapD<ROIradius;
                    TouchID=find(inROI,1,'first');
                    ReleaseID=find(inROI,1,'last');
                    
                catch
                    gapD= DistanceBtwCoor(x2,y2,ROIx1,ROIy1);
                    inROI=gapD<ROIradius;       %ROI size
                    TouchID=find(inROI,1,'first');
                    ReleaseID=find(inROI,1,'last');
                end
                try
                    [GapDistance(k),GapFrame(k)]=min(gapD);
                catch
                end
%                 if GapDistance(k)>12 && j==174
%                     123
%                 end

%                 figure(21),temp=P_Frame(k);
%                 imshow(extractframe(Exp,1,j,temp,1));
%                 hold on,
%                 for m=1:Pwin
%                     scatter(x2(m),y2(m))
%                 end
%                 line([ROIx1 ROIx2],[ROIy1 ROIy2])

%                 figure(21),plot(inROI)
                slope=(ROIy2-ROIy1)/(ROIx2-ROIx1);
                c=ROIy1+ROIradius-slope*ROIx1;
 %-----------------------
                
%more filters:
%has enter ROI (F3)  Pnum length
EnterROI(k)=~isempty(TouchID);
%tolerate 2 frames outside ROI from touch to release due to bad tracking (F4)
NotManyOutsideROI(k)=sum(inROI(TouchID:ReleaseID)==0)<=4;
%have>=3 frames inROI from touch to release (F5)
EnoughInROI(k)=sum(inROI(TouchID:ReleaseID)==1)>=3;
%no slipping off (F6) **if there is obvious slip F4 would have found it,
%this is just for finding slip with a following protraction that occurs below ROI
% NoSlip(k)=~belowROI(ReleaseID+1,x2,y2,slope,c);
NoSlip(k)=~belowROImostly(x2_shift,y2_shift,slope,c);
% work in progress (F7) : check curvature change at touch

%                 if ~isempty(find(inROI,1,'first'))  && sum(inROI(find(inROI,1,'first'):end)==0)<=1  && sum(inROI(find(inROI,1,'first'):find(inROI,1,'last'))==1)>=4
                if LargeWhisk(k) && AfterPiston(k) && EnterROI(k) && NotManyOutsideROI(k) && EnoughInROI(k) && NoSlip(k)
                    T_Frame(k)=P_Frame(k)+TouchID-1;
                    T_End(k)=P_Frame(k)+ReleaseID-1;
%                     if T_End(k)-T_Frame(k)>12
%                         figure
%                         frame=extractframe(Exp,1,j,T_Frame(k),0);
%                         imshow(frame)
%                         playVid(Exp,[1],j,[],[],0,[T_Frame(k)-5 T_End(k)+5],0)  %playVid(Exp,CamN,VidN,scat,Whisker,playAngle,duration,slow)
% 
%                     end
                elseif poorTrackingContingency(i)
                    if AfterPiston(k) && mean(WT.CenAngle.sig(P_Frame(k):P_Frame(k)+Pwin-1)>th)>percentOverth       %if more than (percentOverth%) frames are >th rad
                        temp=WT.CenAngle.sig(P_Frame(k):P_Frame(k)+Pwin-1);
                        TouchID=find(temp>th,1,'first');
                        ReleaseID=find(temp>th,1,'last');
                        T_Frame(k)=P_Frame(k)+TouchID-1;
                        T_End(k)=P_Frame(k)+ReleaseID-1;
                    end
                end
            end
        GoodTouch=~isnan(T_Frame); 
        WT.P_Frame=P_Frame;
        WT.R_Frame=R_Frame;
        WT.TouchFrame=T_Frame(GoodTouch);
        WT.ReleaseFrame=T_End(GoodTouch);
        %not saving these since theyre the same as LocsP(GoodTouch)
%         WT.ProtractionFrame= P_Frame(GoodTouch);         %only protractions with touch are kept
%         WT.RetractionFrame= R_Frame(GoodTouch);
        WT.GoodTouch=GoodTouch;
        WT.GapDistance=GapDistance;
        WT.GapFrame=GapFrame;
        WT.TouchGap=diff([0;find(GoodTouch==1)])-1;
        WT.GoodType=[LargeWhisk AfterPiston EnterROI NotManyOutsideROI EnoughInROI NoSlip];
        Whisker(i).trial(j)=WT;
        fprintf('Touchtime found for Whisker %d Trial %d Ntouch %d/%d\n',i,j,sum(GoodTouch),length(GoodTouch))
end
AllTouchesDur=[];
for j=1:Exp.TrN
    WT=Whisker(i).trial(j);
    AllTouchesDur=[AllTouchesDur;WT.ReleaseFrame-WT.TouchFrame];
end
figure
histogram(AllTouchesDur,0:1:20)
title(sprintf('Whisker %d n=%d',Wid,length(AllTouchesDur)))
end


%% Nested functions
    function B=belowROI(id,x,y,slope,c)
        try
            if id<=length(y)

                if y(id) > slope*x(id)+c
                    B=1;
                else
                    B=0;
                end          
            else
                B=0;
            end
        catch
            B=0;
        end
    end
    function B=belowROImostly(x,y,slope,c)
        Ball=zeros(1,length(y));
        for mm=1:length(y)
            if y(mm) > slope*x(mm)+c
                Ball(mm)=1;
            else
                Ball(mm)=0;
            end
        end
        B=sum(Ball)/length(Ball)>0.3;
    end
end
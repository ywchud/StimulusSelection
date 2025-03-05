function VidParams=playVid(Exp,CamN,VidN,scat,Whisker,playAngle,duration,slow,p)
makeGif=0;
if makeGif==1
    disp('Continue to make a Gif')
end
if exist('p','var') && isfield(p,'wlabel')
    wlabel=p.wlabel;
else
    wlabel=1;
end
%arg
%duration: must be a vector of size 2 in the following format => [startFrame endFrame] or [(0<startPorportion<1) FrameN]
%slow: interger for how much the video is slowed, 0 is no slow
%playAngle: still working on it
%scat: matrix of Nx2 that includes N points of (x,y) to be scattered
NumOfCam=length(CamN);

for i=1:NumOfCam
    datapath=Exp.Path.vid{CamN(i)};
    Labelpath=fullfile(datapath,'Labelled');
    videoname = sprintf(Exp.Path.vidName{CamN(i)}(1:end-4),VidN);
    if isfield(Exp.Path,'csvName') && wlabel
        dataname = sprintf(Exp.Path.csvName{CamN(i)},VidN);
        try
            V{i} = VideoReader(fullfile(Labelpath,[videoname '_labeled.mp4']));
        catch
            write_labed_video(datapath,dataname,[videoname '.mp4'],Exp.DLC);
            V{i} = VideoReader(fullfile(Labelpath,[videoname '_labeled.mp4']));
        end
    elseif isfield(Exp.Path,'csvName')
        V{i} = VideoReader(fullfile(datapath,[videoname '.mp4']));
    else
        disp('No csv found, playing unlabelled original video')
        V{i} = VideoReader(fullfile(datapath,[videoname '.mp4']));
    end

    
end

%get frame number and frame size, assume all cameras have the same
%dimension, so only look at first camera
frameN=ceil(V{1}.Duration*V{1}.FrameRate);
try
    trash=read(V{1},frameN);
catch
    frameN=frameN-1;
end
VidParams.frameN=frameN; 
VidParams.Width=V{1}.Width;
VidParams.Height=V{1}.Height;

if ~isempty(scat)
    if size(scat,2)~=2
        error('Expecting a scatter matrix with N rows and 2 columns (x,y)')
    end
    scat=scat(~isnan(scat(:,1)) & ~isnan(scat(:,2)),:);
end

if ~isempty(duration)   
    if length(duration)~=2
        disp('Expecting a duration vector of size 2, [startFrame endFrame] or [startFrame FrameN] or [(0<startPorportion<1) FrameN]')
        disp('Playing the entire video instead...')
        startF=1;
        endF=frameN;
    elseif duration(1)<1 && duration(1)>0 %[(0<startPorportion<1) FrameN]
        startF=ceil(frameN*duration(1));
        endF=startF+duration(2);
    elseif duration(1)<duration(2)  %[startFrame endFrame]
        startF=max([1 duration(1)]);  %in case user gave 0, change to 1
        endF=duration(2);
    elseif duration(1)>=duration(2)  %[startFrame FrameN]
        startF=max([1 duration(1)]);  %in case user gave 0, change to 1
        endF=startF+duration(2);
    end
    
else  
    startF=1;
    endF=frameN;
end
Wn=0;AngleWin=[-30:30];
if isstring(playAngle) || ischar(playAngle)
    Wn=length(Whisker);
    fprintf('%d whiskers found\n',Wn)
    for k=1:Wn
        breakearly=0;
        WT=Whisker(k).trial(VidN);
        switch playAngle
            case 'Curvature'
                plotSig{1,k}=WT.Curvature;
            case 'CenAngle'
                plotSig{1,k}=WT.CenAngle.sig;
            case 'Phase'
                plotSig{1,k}=WT.Phase;
            case 'Setpoint'
%                 plotSig{1,k}=WT.Setpoint;
                plotSig{1,k}=minmaxnorm(WT.Setpoint);
            case 'CurvDebase'
                plotSig{1,k}=WT.Curvature;
                plotSig{2,k}=WT.CurvDebase;
            case 'CheckTouchAngle'
                plotSig{1,k}=WT.CenAngle.sig;
                CheckPro=WT.LocsT(WT.GoodTouch);
                shiftByOne=[0;WT.GoodTouch(1:end-1)];
                CheckEnd=WT.LocsT(logical(shiftByOne));  %the next protraction frame for use as end of each cycle
                if length(CheckEnd)<length(CheckPro)
                    CheckEnd=[CheckEnd;WT.FrameN];
                end
                CheckRet=WT.LocsP(WT.GoodTouch);
                CheckTouch=WT.TouchFrame;
                CheckRelease=WT.ReleaseFrame;
                HighlightP{1,k}={[CheckPro CheckEnd],'k'};
                HighlightP{2,k}={[CheckTouch CheckRelease],'r'};
            case 'CheckTouchPhase'
                plotSig{1,k}=WT.Phase;
                CheckPro=WT.LocsT(WT.GoodTouch);
                shiftByOne=[0;WT.GoodTouch(1:end-1)];
                CheckEnd=WT.LocsT(logical(shiftByOne));  %the next protraction frame for use as end of each cycle
                if length(CheckEnd)<length(CheckPro)
                    CheckEnd=[CheckEnd;WT.FrameN];
                end
                CheckRet=WT.LocsP(WT.GoodTouch);
                CheckTouch=WT.TouchFrame;
                CheckRelease=WT.ReleaseFrame;
                HighlightP{1,k}={[CheckPro CheckEnd],'k'};
                HighlightP{2,k}={[CheckTouch CheckRelease],'r'};
            case 'DeltaPhase'
                if Wn~=2
                    error('Expecting 2 whiskers only')
                end
                plotSig{1,1}=WT.Phase;                
                WT2=Whisker(k+1).trial(VidN);
                plotSig{1,2}=WT2.Phase;
                plotSig{1,3}=WT.findPhaseDiff(WT2,0);
                breakearly=1;
            case 'RetractionBasedSyncI'
                if Wn~=2
                    error('Expecting 2 whiskers only')
                end
                [SYNCI,PHASE,AMP]=findRetractionBasedSyncIndex(Exp,Whisker,VidN);
                plotSig{1,1}=SYNCI{1};
                plotSig{2,1}=zeros(WT.FrameN,1);
                [SYNCI,PHASE,AMP]=findRetProBasedSyncIndex(Exp,Whisker,VidN);
                plotSig{1,2}=SYNCI{1};
%                 plotSig{1,2}=PHASE{1};
                plotSig{2,2}=zeros(WT.FrameN,1);
%                 plotSig{1,3}=AMP{1};
%                 plotSig{2,3}=zeros(WT.FrameN,1);
                breakearly=1;
            case 'Sync'  %[0 1]
                if Wn~=2
                    error('Expecting 2 whiskers only')
                end
                plotSig{1,1}=Exp.Sync{VidN};
                plotSig{2,1}=zeros(WT.FrameN,1)+0.8;
                yrange=[0 1];
            case 'Directionality' %[-1 1]
                if Wn~=2
                    error('Expecting 2 whiskers only')
                end
                plotSig{1,1}=Exp.Directionality{VidN};
                plotSig{2,1}=zeros(WT.FrameN,1);
                yrange=[-1 1];
            case 'Spatiality' %[-1 1]
                if Wn~=2
                    error('Expecting 2 whiskers only')
                end
                plotSig{1,1}=Exp.Spatiality{VidN};
                plotSig{2,1}=zeros(WT.FrameN,1);
                yrange=[-1 1];
            otherwise
                error('Signal type %s not found',playAngle)
        end
        if exist('Exp.RunSpeed','var')
            sig=find(Exp.RunSpeed{VidN} >0.25);
            %use sig.onset
            ShadeP{1,k}={[RunS RunEnd],'k'};
        end
        if breakearly
            break
        end
    end
    
end

FrameNo=0;
% figure('WindowState', 'maximized')
 for i=startF:endF  
     for j=1:NumOfCam
         try
            frame{j} = read(V{j},i);
         catch
             disp('If you see this msg please check the code. There may a frame number mismatch')
             break;
         end
        hold off
        %determine plotting dimension
        
        
        if ~exist('plotSig','var')  %only vids no signal plot
            switch NumOfCam
                case 1
                    si=1;sj=1;
                    vidPlotSpace=1;
                case 2
                    si=1;sj=2;
                    vidPlotSpace=j;
                otherwise
                    error('Please define plotting dimension for >2 camera')
            end            
        else  %both vids and signal plot
            switch NumOfCam
                case 1
                    si=4;sj=2;
                    vidPlotSpace=[1 3 5 7];
                case 2
                    si=4;sj=2;
                    vidPlotSpace=[1 3]+4*(j-1);
                otherwise
                    error('Please define plotting dimension for >2 camera')
            end      
        end

            
%         if isstring(playAngle) 
%             si=ceil(size(plotSig,2)/NumOfCam)+1;
%             sj=NumOfCam;
%         elseif  NumOfCam>1
%             si=1;
%             sj=NumOfCam;
%         else
%             si=1;
%             sj=1;
%         end
        subplot(si,sj,vidPlotSpace)
        imagesc(frame{j}),axis('image'),title(sprintf('Cam %d Trial %d Frame= %d/%d',CamN(j),VidN,i-startF+1,endF-startF+1))
        hold on
        for k=1:Wn
            ROI_overlay(Exp,Whisker(k).W_id,VidN)
        end
    end
    if ~isempty(scat)
        scatterThis(scat)
    end
    
    
    
    
    if isstring(playAngle) || ischar(playAngle)
        SigSpaceNeeded=size(plotSig,2);
        for k=1:SigSpaceNeeded
            %%%%
            switch SigSpaceNeeded
                case 1
                    SigPlotSpace=[2 4 6 8];
                case 2
                    SigPlotSpace=[2 4]+4*(k-1);
                case {3, 4}
                    SigPlotSpace=[2]+2*(k-1);
                otherwise
                     error('Please define plotting dimension for >4 signals')
            end
            
            subplot(si,sj,SigPlotSpace)
            win=i+AngleWin;
            inbounds=(win>=1 & win<=frameN);
            win=win(inbounds);
            
            hold off
            x=AngleWin(inbounds);
            for m=1:size(plotSig,1)   %plot unique signals
                plot(x,plotSig{m,k}(win))
                hold on
            end
            if exist('HighlightP','var')  %highlight segments of said unique signals
                for m=1:size(HighlightP,1)
                    Highlight(HighlightP{m,k},plotSig{1,k},win,x);
                end
            end
            line([0 0],ylim,'Color','red')
            if exist('yrange','var') %ylim
                ylim(yrange)
            end
            
            if k<=Wn
                title(sprintf('Whisker %d Cam %d',Whisker(k).W_id,Whisker(k).Cam))
            end
        end
    end   
    FrameNo=FrameNo+1;
    if makeGif
        if FrameNo==1
            filename= [Exp.Path.save  '\\TouchCheck.gif'];
            gif(filename,'frame',gcf)
        else
            gif('frame',gcf)
        end
    end
    pause(1/30*slow)

 end
 

end
function scatterThis(scat)

    for ii=1:size(scat,1)
        scatter(scat(ii,1),scat(ii,2));
        hold on
    end
end 

function Highlight(HighlightP,plotSig,win,X)
    
P=HighlightP{1};
color=HighlightP{2};


    Hstart=find(P(:,1)>=win(1) &  P(:,1)<=win(end));

    for i=1:length(Hstart)
        ID=Hstart(i);
        HLseg=P(ID,1):P(ID,2);
        x=HLseg-win(1)+1;
        
        HLseg=HLseg(x<=length(win)); %in case segment exceed window size
        x=x(x<=length(win));
        y=nan(length(win),1);
        y(x)=plotSig(HLseg);
        plot(X,y,'Color',color)
    end

end

function Shade(ShadeP,plotSig,win,X)
 %%not implemented yet
    P=ShadeP{1};
    color=ShadeP{2};


    Hstart=find(P(:,1)>=win(1) &  P(:,1)<=win(end));

    for i=1:length(Hstart)
        ID=Hstart(i);
        HLseg=P(ID,1):P(ID,2);
        x=HLseg-win(1)+1;
        
        HLseg=HLseg(x<=length(win)); %in case segment exceed window size
        x=x(x<=length(win));
        y=nan(length(win),1);
        y(x)=plotSig(HLseg);
        plot(X,y,'Color',color)
    end

end

function ROI_overlay(Exp,w,tr)
    hold on
    trI=find(Exp.Stim.Piston.ROI(w).index_trial==tr);
    if ~isempty(trI)
        try
            line([Exp.Stim.Piston.ROI(w).x1_intp(trI,2) Exp.Stim.Piston.ROI(w).x2_intp(trI,2)],[Exp.Stim.Piston.ROI(w).y1_intp(trI,2) Exp.Stim.Piston.ROI(w).y2_intp(trI,2)])
        catch
            line([Exp.Stim.Piston.ROI(w).x1_intp(trI,1) Exp.Stim.Piston.ROI(w).x2_intp(trI,1)],[Exp.Stim.Piston.ROI(w).y1_intp(trI,1) Exp.Stim.Piston.ROI(w).y2_intp(trI,1)])
            
        end
        
    end

end
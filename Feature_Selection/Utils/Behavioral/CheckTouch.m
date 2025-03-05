function CheckTouch(Exp,Whisker,StimC,I,p)
%(Exp,Whisker,StimC,I)
%I:vector of length 2: between 0 & 1 determines proportion of touch to check, e.g. [0 1] checks all
% or vector of n touch indexes to check
if exist('p','var') && isfield(p,'checkType') 
    checkType=p.checkType;
else
    checkType="CheckTouchAngle";  %CheckTouchPhase
end
if exist('p','var') && isfield(p,'randomized') 
    randomized=p.randomized;
else
    randomized=0;  %CheckTouchPhase
end
if exist('p','var') && isfield(p,'video') 
    video=p.video;
else
    video=1;  %CheckTouchPhase
end
if exist('p','var') && isfield(p,'maxN') 
    maxN=p.maxN;
else
    maxN=[];  %CheckTouchPhase
end                           
 

try
    TF=cell2mat(StimC.TouchFrames);
catch
    TF=cell2mat(TouchFrameFromTime(Exp,StimC));
end
try
    CT=StimC.CorrTrials;
catch
    CT=StimC.CorrTrial{1};
end

if ~isempty(I)
    if length(I)==2 && I(1)>=0 && I(1)<=1 && I(2)>=0 && I(2)<=1   %vector of length 2: between 0 & 1 determines proportion of touch to check, e.g. [0 1] checks all
        II=round(I(1)*length(TF))+1:round(I(2)*length(TF));
        TF=TF(II);
        CT=CT(II);
    else   %vector of n touch indexes to check
        TF=TF(I);
        CT=CT(I);
    end
end
if randomized
    I=randperm(length(TF));
    TF=TF(I);
    CT=CT(I);
end
if ~isempty(maxN) && length(TF)>maxN
    TF=TF(1:maxN);
    CT=CT(1:maxN);
end

try
    CamN=unique(Exp.Stim.Piston.Cam(logical(StimC.WhiskerBool)));
catch
    WB=false(1,Exp.Stim.Piston.num);WB(StimC.WID)=1;
    CamN=unique(Exp.Stim.Piston.Cam(logical(WB)));
end
if video
    vidWin=[-10 30];
    for i=1:length(TF)
        VidN=CT(i);
        Frame=TF(i);
    %     figure

        playVid(Exp,CamN,VidN,[],Whisker,checkType,vidWin+Frame,3);
        pause(0.5);
    %     close
    end
else  %image
    Fs=cell(length(TF),1);
    [si,sj]=subplotDim(length(TF));
    figure
    for i=1:length(TF)
        VidN=CT(i);
        Frame=TF(i);
        Fs{i}=extractframe(Exp,CamN,VidN,Frame,0);
        subplot(si,sj,i)
        image(Fs{i})
    end  
end
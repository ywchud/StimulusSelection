function ProtractionFrame=TimeToFrame(PTtype,TrialStartTime,TrialEndTime,videoFps,WhiskerError,CamDelay)

disp('Use TouchFrameFromTime instead of TimeToFrame(obsolete) for consistency')
    ProtractionFrame=[];
%     CamDelay=0.05;  %old rig
%     CamDelay=0.15;  %bilateral
    for j = 1:length(TrialStartTime)
    temp=[];
    for i = 1:length(PTtype)
        if PTtype(i)>= TrialStartTime(j)+CamDelay && PTtype(i)<= TrialEndTime(j)+2/videoFps  %adding 2 other frametime for TrialEndTime to remove edge effect, e.g. the time coverted from the last frame of video may be larger than TrialEndTime due to resolution error in frames and CamDelay 
            Frame=ceil((PTtype(i)-TrialStartTime(j)-CamDelay)*videoFps-WhiskerError(j));
            if Frame<1
                fprintf('Time(%d),Frame(%d),Trial(%d)\n',PTtype(i),Frame,j)
            else
                temp=[temp; Frame];
            end
        end
    end
    ProtractionFrame{j}=temp;
    end
end
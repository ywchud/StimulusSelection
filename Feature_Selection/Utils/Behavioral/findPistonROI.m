function [Exp,ESP]=findPistonROI(Exp,WOI,duration,slow)
%%  Init
%determine fraction: default is 4
%the number of times (fraction+1) a ROI is sampled across trials for associated piston of WOI (used for interpolation to address drifting piston & to find mean ROI)   
    ESP=Exp.Stim.Piston;
    
    try
        if ~isempty(ESP.ROIfraction)
            fraction=ESP.ROIfraction;
        else
            fraction=4;
        end
    catch
        ESP.ROIfraction=4;
        fraction=ESP.ROIfraction;
    end


    array=zeros(1,ESP.tot);array(WOI)=1;
    n=find(PistonComb(array,ESP.Mat),1,'first');
%     N=find(PistonComb(array,ESP.Mat));   %trials with WOI alone
    NwithWOI= find(~isnan(ESP.Mat(:,WOI*2)));  %trials with WOI w/wo other whiskers
    if isempty(NwithWOI)
        return
    end
    if isempty(n)
        n=NwithWOI(1);
    end
    ESP.ROI(WOI).TrN=length(NwithWOI);
    
    breakpoints=zeros(length(NwithWOI),1);
    breakpointsType=strings(length(NwithWOI),1);
    ROIs_x1=nan(length(NwithWOI),1);
    ROIs_y1=nan(length(NwithWOI),1);
    ROIs_x2=nan(length(NwithWOI),1);
    ROIs_y2=nan(length(NwithWOI),1);

    
    scat=[ROIs_x1, ROIs_y1;ROIs_x2, ROIs_y2];
    
    %% show which piston to find ROI for
    CamN=ESP.Cam(WOI);
    datapath=Exp.Path.vid{CamN};
    dataname = sprintf(Exp.Path.csvName{CamN},n);
    videoname = sprintf(Exp.Path.vidName{CamN}(1:end-4),n);  
    try
        V = VideoReader(fullfile(datapath,[videoname '.mp4']));
    catch
        error('rip');
    end
    frameN=ceil(V.Duration*V.FrameRate);
    frame=read(V,round(frameN/2));
    figure
    imagesc(frame),axis('image'),title('Finding ROI for the following piston, press anywhere to continue')
    [~ ,~]=ginput(1);
    close
    %% loop through fractions of trials to find ROI in sequence
    GrpN=1;drift=0;
for frac=1:fraction+1

    if frac==1
        id=1;
        VidN=NwithWOI(id); 
        PreviousID=id;
    else
        PreviousID=id;
        id=ceil(length(NwithWOI)/fraction*(frac-1));
        VidN=NwithWOI(id);
    end
    breakpoints(id)=1;
    breakpointsType(id)="Frac";
    CamN=ESP.Cam(WOI);
    
%%  Totally could have made a function for picking ROI but I didn't. Pray you don't have to debug this part
    ok=0;d=duration;
    fprintf('Find ROI for fraction: %d/%d\n',frac,fraction+1)
    while(~ok)
        ID_last=find(NwithWOI==VidN);
        ID_first=PreviousID;
        win=ID_first:ID_last;  %window defined as the indexes of trials with WOI being looked at confined in this fraction
        VidParams=playVid(Exp,CamN,VidN,scat,[],0,d,slow);  %playVid(Exp,CamN,VidN,scat,Whisker,playAngle,duration,slow)
        title(sprintf('Pick ROI vertex 1/2 for Piston: %d, Frac: %d/%d',ESP.ID(WOI),frac,fraction+1))
    %     [touchCoor(WOI).x1(frac,1),touchCoor(WOI).y1(frac,1)]=ginput(1);
        [ROIs_x1(ID_last),ROIs_y1(ID_last)]=ginput(1);
        title(sprintf('Pick ROI vertex 2/2 for Piston: %d, Frac: %d/%d',ESP.ID(WOI),frac,fraction+1))
    %     [touchCoor(WOI).x2(frac,1),touchCoor(WOI).y2(frac,1)]=ginput(1);
        [ROIs_x2(ID_last),ROIs_y2(ID_last)]=ginput(1);
        title(sprintf('Click image if OK, Left if piston drifted,Right to shift window and restart'))
        [check1,~]=ginput(1);
          close
        if check1<0 && frac>1   %if drifted, further segment this fraction into groups with grouptype(=same or linear), then find ROI for every last trial in each subgroup
            
            CheckDriftVidN=NwithWOI(win);
            
            scat=[ROIs_x1, ROIs_y1;ROIs_x2, ROIs_y2];
            scatID=[win win+length(NwithWOI)];
            scat_seg=scat(scatID,:);
            
            [groups,groupType]=findDriftTrial(Exp,CamN,CheckDriftVidN,scat_seg);
            groups=toColumn(groups);
            needRefineROI=logical([diff(groups); 0]);
            needRefineROIid=find(needRefineROI);
            RedoTrials=NwithWOI(win(needRefineROI));
            fprintf('Find ROI for additional trials: [')
            fprintf('%g ', RedoTrials);
            fprintf(']\n');
            
            for k = 1:length(RedoTrials)
                ok2=0;d2=duration;
                while(~ok2)
                    VidParams=playVid(Exp,CamN,RedoTrials(k),scat_seg,[],0,d2,slow);  %playVid(Exp,CamN,VidN,scat,Whisker,playAngle,duration,slow)
                    title(sprintf('Pick ROI vertex 1/2 for Piston: %d, Frac: %d/%d, SubFrac: %d/%d',ESP.ID(WOI),frac,fraction+1,k,length(RedoTrials)))
                %     [touchCoor(WOI).x1(frac,1),touchCoor(WOI).y1(frac,1)]=ginput(1);
                    [ROIs_x1(win(needRefineROIid(k))),ROIs_y1(win(needRefineROIid(k)))]=ginput(1);
                    title(sprintf('Pick ROI vertex 1/2 for Piston: %d, Frac: %d/%d, SubFrac: %d/%d',ESP.ID(WOI),frac,fraction+1,k,length(RedoTrials)))
                %     [touchCoor(WOI).x2(frac,1),touchCoor(WOI).y2(frac,1)]=ginput(1);
                    [ROIs_x2(win(needRefineROIid(k))),ROIs_y2(win(needRefineROIid(k)))]=ginput(1);
                    title(sprintf('Click image if OK, click elsewhere to restart'))
                    [check2,~]=ginput(1);
                    if check2<=VidParams.Width && check2>=0
                        scat=[ROIs_x1, ROIs_y1;ROIs_x2, ROIs_y2];
                        scat_seg=scat(scatID,:);
                        ok2=1;
                    else
                        ok2=0;
                        d2=shiftDur(d2);
                    end
                    close
                end
            end
            temp=groups+GrpN-1;
            breakpoints(win(2:end))=temp(2:end);
            breakpointsType(win(2:end))=groupType(2:end);
            GrpN=GrpN+max(groups);
            drift=1;
            ok=1;
        elseif check1<=VidParams.Width && check1>=0  %good ROI,no significant drift, move on
            scat=[ROIs_x1, ROIs_y1;ROIs_x2, ROIs_y2];
            breakpoints(win(2:end))=GrpN;
            breakpointsType(win(2:end))="same";
            GrpN=GrpN+1;  
            ok=1;
        else   %not ok, shift window to look at other segment of the video to re-find ROI
            ok=0;
            d=shiftDur(d);  
        end
        
        
    end  
end
%find trial-dependent ROI based on given breakpoints

fprintf('Interpolating for Piston %d\n',ESP.ID(WOI))
diffB=toColumn(diff(breakpoints));
endID=find([diffB;1]~=0);
startID=find([1;diffB]~=0);


for i=1:length(startID)
    if strcmp(breakpointsType(endID(i)),'same')
        
        if isempty(ROIs_x1(startID(i))) || isnan(ROIs_x1(startID(i)))
            wn=startID(i):endID(i);
        else
            wn=startID(i)+1:endID(i);
        end
        
        ROIs_x1(wn)=ROIs_x1(endID(i));
        ROIs_y1(wn)=ROIs_y1(endID(i));
        ROIs_x2(wn)=ROIs_x2(endID(i));
        ROIs_y2(wn)=ROIs_y2(endID(i));
    elseif strcmp(breakpointsType(endID(i)),'linear')
        if isempty(ROIs_x1(startID(i))) || isnan(ROIs_x1(startID(i)))
            wn=startID(i)-1:endID(i);
        else
            wn=startID(i):endID(i);
        end
        
        ROIs_x1(wn)=interp1([1 length(wn)],[ROIs_x1(wn(1)) ROIs_x1(wn(end))],1:length(wn));
        ROIs_x2(wn)=interp1([1 length(wn)],[ROIs_x2(wn(1)) ROIs_x2(wn(end))],1:length(wn));
        ROIs_y1(wn)=interp1([1 length(wn)],[ROIs_y1(wn(1)) ROIs_y1(wn(end))],1:length(wn));
        ROIs_y2(wn)=interp1([1 length(wn)],[ROIs_y2(wn(1)) ROIs_y2(wn(end))],1:length(wn));
        
    end
end

ESP.ROI(WOI).isdrift=logical(drift);
ESP.ROI(WOI).index_group=breakpoints;
ESP.ROI(WOI).index_type=breakpointsType;
ESP.ROI(WOI).index_trial=NwithWOI;

ESP.ROI(WOI).x1_intp = ROIs_x1;
ESP.ROI(WOI).x2_intp = ROIs_x2;
ESP.ROI(WOI).y1_intp = ROIs_y1;
ESP.ROI(WOI).y2_intp = ROIs_y2;

ESP.ROI(WOI).x1=mean(ROIs_x1);
ESP.ROI(WOI).x2=mean(ROIs_x2);
ESP.ROI(WOI).y1=mean(ROIs_y1);
ESP.ROI(WOI).y2=mean(ROIs_y2);
fprintf('Mean and Interpolated ROIs are saved for Piston %d\n',ESP.ID(WOI))

Exp.Stim.Piston=ESP;
end



function d_new=shiftDur(d)
    if d(1)>0 && d(1)<1
        d_new(1)=d(1)+0.05;
        d_new(2)=d(2);
        if d_new(1)>=1
            error('Out of Frames, just retry another video')
        end
        
    else
        d_new(1)=d(2);
        d_new(2)=2*d(2)-d(1);
    end

end


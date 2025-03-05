function Exp=MechanicalDelayAdjustment_Auto(Exp,bpy,bpx,Ymat,Xmat)
PistonN=Exp.Stim.Piston.tot;
sampleN=30;
if ~isempty(Exp.PistonBuffer)
%     error('PistonBuffer is not empty. To continue, comment this line')
end
Exp.PistonBuffer=[];

if isempty(bpy) || isempty(bpx)
    BPy=cell(1,PistonN);   
    BPx=cell(1,PistonN);

for i=1:PistonN
    if ~Exp.Stim.Piston.active(i)
        continue
    end
    fprintf('Finding cues for piston %d...\n',i)
    %random sampling
    Piston_i_tr=find(~isnan(Exp.Stim.Piston.Mat(:,i*2)));
    randid=randperm(length(Piston_i_tr),min([sampleN length(Piston_i_tr)]));
    selectedtrials=Piston_i_tr(randid);
    %get dim
    tr=selectedtrials(1);
    temp=extractframe(Exp,Exp.Stim.Piston.Cam(i),tr,round(Exp.FrameN(tr)/2),0);
    temp=rgb2gray(temp);
    [H,W]= size(temp);
    F=zeros(H,W,length(selectedtrials));
    %get averaged piston i frame
    for m=1:length(selectedtrials)
        tr=selectedtrials(m);
        f=extractframe(Exp,Exp.Stim.Piston.Cam(i),tr,round(Exp.FrameN(tr)/2),0);
        try
            F(:,:,m)=rgb2gray(f);
        catch
            disp('Problem with extracted frame, skipping...')
        end
        fprintf('Sampling %d/%d\n',m,min([sampleN length(Piston_i_tr)]));
    end
    FF=uint8(mean(F,3));
    %
    x=[]; y=[];
    for m=1:length(selectedtrials)
        tr=selectedtrials(m);
        matID=[1:Exp.DLC.LblPerW]+Exp.DLC.LblPerW*(i-1);
        tip=1;
        x=[x;Xmat{tr}(:,matID(tip))];
        y=[y;Ymat{tr}(:,matID(tip))];
    end
            
    %Get ROI
    figure
    imshow(FF),title(sprintf('Find overlapping roi w/whisker for piston %d',i))
    hold on
    scatter(x,y,'.')
    roi = drawpolygon;
    ypos=roi.Position(:,1);
    xpos=roi.Position(:,2);
    close
    
    %find bounded pixels
    boundedpixels=[];
    for h=1:H
        for w=1:W
            in = inpolygon(w,h,ypos,xpos);
            if in
                boundedpixels=[boundedpixels;h w];
            end
        end
    end
    BPy{i}=boundedpixels(:,1);
    BPx{i}=boundedpixels(:,2);
end
assignin('base','BPy',BPy) 
assignin('base','BPx',BPx) 
else
    BPy=bpy;
    BPx=bpx;
end


for i=1:PistonN
    if ~Exp.Stim.Piston.active(i)
        continue
    end
    tic
    fprintf('Finding delay frames for piston %d...\n',i)
    %find Piston ROI entry frame per trial 
    Piston_i_tr=find(~isnan(Exp.Stim.Piston.Mat(:,i*2)));
    L=length(Piston_i_tr);
    PBuffer=nan(L,1);
    PReturnBuffer=nan(L,1);
    BPyi=BPy{i};
    BPxi=BPx{i};
    
    %get piston out frame
    delayF=round((Exp.Stim.Piston.Delay-Exp.Cam.Delay)*Exp.videoFps);checkwin=150;
    for m=1:L
        tr=Piston_i_tr(m);     
        if length(delayF)==1
            delayFtr=delayF(1);
        else
            delayFtr=delayF(tr);
        end       
        I=[];
        while isempty(I)
            test=extractframe(Exp,Exp.Stim.Piston.Cam(i),tr,[delayFtr:delayFtr+checkwin],0);
            ROIcolor=zeros(1,checkwin);
            for n=1:checkwin
                %                 f=zeros(H,W);
                f=rgb2gray(test(:,:,:,n));
                ind = sub2ind(size(f),BPyi ,BPxi);
                ROIcolor(n)=mean(f(ind));
            end
            [temp,I]=findpeaks(abs(diff(ROIcolor)),'MinPeakHeight',5);
            if isempty(I) 
                %             [temp,I]=max(abs(diff(ROIcolor)));
                checkwin=checkwin+100;
                fprintf('No clear peak, extending window to %d...\n',checkwin)
            end
            if checkwin>400
                error('Check')
                %       	[temp,I]=max(abs(diff(ROIcolor)));
            end
        end
%         if max(temp)<10
%             figure
%             plot(ROIcolor)
%             hold on
%             plot(abs(diff(ROIcolor)))
%         end
        if ~isempty(I)
            z=I(end)+1; %+1 to get the after peak position      
        else
            z=1;
        end
        PBuffer(m)=delayFtr+z-1;
%         fprintf('Piston %d trial %d/%d\n',i,m,L)
        %             figure(1)
        %             plot(ROIcolor),hold on
        %             scatter(z,ROIcolor(z)),hold off
    end
    Exp.PistonBuffer(Piston_i_tr,i)=PBuffer;
    
    %get piston back frame
    fprintf('Finding delay return frames for piston %d...\n',i)
    delayF=round(Exp.FrameN-Exp.videoFps.*(Exp.Stim.Piston.OffWin-Exp.Cam.OffWins));   
    for m=1:L
        checkwin=200;
        tr=Piston_i_tr(m);     
        if length(delayF)==1
            delayFtr=delayF(1);
        else
            delayFtr=delayF(tr);
        end
        if length(Exp.videoFps)>1 && Exp.videoFps(tr)<mode(Exp.videoFps)  %since this tr is missing frames, we look at the whole seg. till end of video
            test=extractframe(Exp,Exp.Stim.Piston.Cam(i),tr,[delayFtr:Exp.FrameN(tr)],0);
        else
            test=extractframe(Exp,Exp.Stim.Piston.Cam(i),tr,[delayFtr:min([delayFtr+checkwin Exp.FrameN(tr)])],0);
        end
        ROIcolor=zeros(1,size(test,4));
        for n=1:size(test,4)
            %                 f=zeros(H,W);
            f=rgb2gray(test(:,:,:,n));
            ind = sub2ind(size(f),BPyi ,BPxi);
            ROIcolor(n)=mean(f(ind));
        end
        try
        [temp,I]=findpeaks(abs(diff(ROIcolor)),'MinPeakHeight',5);
        catch
      	[temp,I]=max(abs(diff(ROIcolor)));
        end
%         if max(temp)<10
%             figure
%             plot(ROIcolor)
%             hold on
%             plot(abs(diff(ROIcolor)))
%         end
        if ~isempty(I)
            z=I(end)+1; %+1 to get the after peak position      
        else
            z=1;
        end
        PBuffer(m)=delayFtr+z-1;
%         fprintf('Piston %d trial %d/%d\n',i,m,L)
        %             figure(1)
        %             plot(ROIcolor),hold on
        %             scatter(z,ROIcolor(z)),hold off
    end
    toc
    Exp.PistonBufferReturn(Piston_i_tr,i)=PBuffer;
    
    
    
end
Exp.PistonBuffer(Exp.PistonBuffer==0)=nan;
Exp.PistonBufferReturn(Exp.PistonBufferReturn==0)=nan;

end

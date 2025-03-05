function touch_frames_sampling(Exp,Whisker)
    

sample_path='E:\Darren\Touch_RNmodel\Samples';


%sampling parameters:

Tr_SampleN=25;   %TrialN sampled for this whisker
TouchPerTr_SampleN=4; %TouchN sampled from each of the above trial
SampleWindowSize=5;   %i.e. combined frameN taken from before and after the current frame collectively as one touch sample

%other parameters:
WhiskerN=length(Whisker);
h=64;
w=64;

SampleN=Tr_SampleN*TouchPerTr_SampleN;  %perwhisker


for i=1:WhiskerN
    W=Whisker(i);CamN=W.Cam;
    Touches=[];CorrTr=[];
    
    for j = 1:Exp.TrN  
        WT=W.trial(j);
        VidN=j;
        Touches=[Touches;WT.TouchFrame]; 
        CorrTr=[CorrTr;j*ones(length(WT.TouchFrame),1)];
    end
    
    AvailT=length(Touches);
    
    if AvailT<SampleN
        error('Not enough touches(%d)',AvailT)
    end
    
    unique()
    
    SampleI=sort(randperm(AvailT,SampleN));
    
    
    
    
end
        
        
        

        
        
        
        
        
    end
end



Sampled_T=Touches(SampleI);
Sampled_CorrTr=CorrTr(SampleI);

Sampled_Tr=unique(Sampled_CorrTr);

for i=1:length(Whisker)
    W=Whisker(i);CamN=W.Cam;
    
    for j=1:length(Sampled_Tr)
        tr=Sampled_Tr(j);
        TFs=Sampled_T(Sampled_CorrTr==tr);
        
        WT=W.trial(tr);
        VidN=tr;
        
        
        
    end
    
    
    
    for j = 1:Exp.TrN  
        WT=W.trial(j);
        VidN=j;

        
        GT=WT.GoodTouch;        
        GT_new=GT;
%         PF=WT.P_Frame;
        TF_old_all=WT.TouchFrame;
        if isempty(TF_old_all)
            continue;
        end
        
        WTP=WT.P_Frame(WT.P_Frame>max(Exp.PistonBuffer(j,:)));
        [~,I]=min(WT.CenAngle.sig(WTP));  %Use most protracted frame as reference_frame
        RefF=WTP(I);
        
        
        
%         ROIx1=0;ROIx2=0;ROIy1=0;ROIy2=0;        
        if Exp.Stim.Piston.ROI(i).isdrift
            ROI_trid=find(Exp.Stim.Piston.ROI(i).index_trial==j);
            ROIx1=Exp.Stim.Piston.ROI(i).x1_intp(ROI_trid);
            ROIx2=Exp.Stim.Piston.ROI(i).x2_intp(ROI_trid);
            ROIy1=Exp.Stim.Piston.ROI(i).y1_intp(ROI_trid);
            ROIy2=Exp.Stim.Piston.ROI(i).y2_intp(ROI_trid);
        else
             ROIx1=Exp.Stim.Piston.ROI(i).x1;
             ROIx2=Exp.Stim.Piston.ROI(i).x2;
             ROIy1=Exp.Stim.Piston.ROI(i).y1;
             ROIy2=Exp.Stim.Piston.ROI(i).y2;
        end
%         
%         for k=1:length(WT.Y)
%             WT.Y(k).sig
%         end
        
        
        
        TF_old_all(randi(length(TF_old_all)))
        
        TF_new=TF_old;
        stillgood=ones(1,length(TF_old));
        for k=1:length(TF_old)
            Frange=TF_old(k)-10:TF_old(k)+10;  %500fps*0.02s=10f
            frames=extractframe(Exp,CamN,VidN,Frange,0);
            frames_gray=[];
            for m=1:size(frames,4)
                frames_gray(:,:,m)=rgb2gray(frames(:,:,:,m));
            end
%             figure(111)
%             imshow(frames_gray(:,:,10),'DisplayRange',[0 250])
            frames_ref=double(rgb2gray(extractframe(Exp,CamN,VidN,RefF,0)));
                        
            [FrameAdj,stillgood(k)]=refineframe(frames_gray,frames_ref,ROIx1,ROIx2,ROIy1,ROIy2,Exp.Stim.Piston.ROIradius(i));
            TF_new(k)=TF_old(k)+FrameAdj;
        end
        WT.TouchFrame=TF_new;
        GT_new(GT)=stillgood;
        WT.GoodTouch=GT_new;
        
        %   CARE
        % gapdistance records closest pixel distance between TLOI and piston line,
        % and gapframe is the corresponding frame starting from current protraction
        % not updated since here we dont use TLOI
        WT.TouchGap=diff([0;find(WT.GoodTouch==1)])-1;
        W.trial(j)=WT;
    end
end




end
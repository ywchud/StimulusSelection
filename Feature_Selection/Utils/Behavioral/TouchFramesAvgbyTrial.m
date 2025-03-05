function Frames_tr=TouchFramesAvgbyTrial(StimC,Exp)
    S=StimC;
    Wid=S.WID(1);
    TouchN=length(S.TouchTimes);
    [G,I]=findgroups(S.CorrTrials);

    [si,sj]=subplotDim(length(I));
    Frames_tr=zeros(Exp.DLC.size(1),Exp.DLC.size(2),length(I));

    figure
    for i=1:length(I)
        subplot(si,sj,i)
        tr=I(i);
        TF=S.TouchFrames{tr};
        if length(TF)~=sum(G==i)
            123; %quality check
        end
        F=zeros(Exp.DLC.size(1),Exp.DLC.size(2),length(TF));
        for j=1:length(TF)
            F(:,:,j)=rgb2gray(extractframe(Exp,Exp.Stim.Piston.Cam(Wid),tr,TF(j),0));
        end
     FF=uint8(mean(F,3,'omitnan'));
     imshow(FF)
     title(sprintf('Trial: %d n=%d',tr,length(TF)))
     Frames_tr(:,:,i)=FF;
    end
end
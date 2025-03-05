function I=FeatureSort(Whisker,StimC)

    TouchN=length(StimC.TouchTimes);
    Curv=nan(TouchN,1);
    RetAngle=nan(TouchN,1);
    
    WID=StimC.WID(1);
    W=Whisker(WID);
    TouchDur=toColumn(StimC.ReleaseTimes-StimC.TouchTimes);
    
    for i=1:TouchN
        tr=StimC.CorrTrials(i);
        RetAngleTi=W.trial(tr).CenAngle.sig;
        CurvTi=W.trial(tr).CurvDebase;
        
        TrialFirstTouch=find(StimC.CorrTrials==tr,1,'first');
        tid=i-TrialFirstTouch+1;
        FrT=StimC.TouchFrames{tr}(tid);
        FrR=StimC.ReleaseFrames{tr}(tid);
        
        RetAngle(i)=max(RetAngleTi(FrT:FrR));
        Curv(i)=max(CurvTi(FrT:FrR));

    end
    mat=[normalize(TouchDur,'range') normalize(RetAngle,'range') normalize(Curv,'range') ];
    
    [~,I]=sort(sum(mat,2));
    figure
%     h=heatmap(mat(I,:));
%     ax = gca;
%     ax.XData = ["TouchDur" "RetAngle" "Curv"];
    imagesc(mat(I,:))
    set(gca,'xtick',[1:3],'xticklabel',{'TouchDur','RetAngle','Curv'})
    colorbar
%     sorty(h,{'RetAngle','Curv','TouchDur'})



end


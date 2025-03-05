function [SAMPLES,run_ratio]=DLCsampling_running(Exp,sampNperCat)
%     sampNperCat=2;
%     runningTrial=(Exp.TrialEndT-Exp.TrialStartT)<7;
    run_ratio=zeros(1,Exp.TrN);
    for i=1:Exp.TrN
        runI=logical(Exp.runStartT>=Exp.TrialStartT(i) & Exp.runStartT<Exp.TrialEndT(i));
        runningT=sum(min([Exp.runEndT(runI) ones(sum(runI),1)*Exp.TrialEndT(i)],[],2)-Exp.runStartT(runI));
        totalT=Exp.TrialEndT(i)-Exp.TrialStartT(i);
        run_ratio(i)=runningT/totalT;
    end
    
    getting_th=1;  th=0.5+0.05;
    while getting_th
        getting_th=0;
    th=th-0.05;
    fprintf('Th: %d\n',th)
    runningTrial=run_ratio>=th;
    I=runningTrial;
    II=find(I);
    mat=Exp.Stim.Piston.Mat;
    mat(~I,:)=[];
    M=~isnan(mat(:,2:2:end));
    G=findgroups(array2table(M));
    samp=zeros(max(G),sampNperCat);   
        for i=1:max(G)
            Gi=find(G==i);
            if length(Gi)<sampNperCat
                fprintf('Fewer samples than %d for group %d\n',sampNperCat,i)
                getting_th=1;
            end
        end
    end
    
    for i=1:max(G)
        Gi=find(G==i);
        if length(Gi)<sampNperCat
            error('Fewer samples than %d for group %d',sampNperCat,i)
        end
        while length(unique(samp(i,:)))~=sampNperCat || sum(samp(i,:))==0
            samp(i,:)=Gi(ceil(rand(sampNperCat,1)*length(Gi)));
        end
    end
    
    SAMPLES=sort(samp(:));
    SAMPLES=II(SAMPLES);
end
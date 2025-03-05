function SAMPLES=DLCsampling(Exp,sampNperCat)
%     sampNperCat=2;
    mat=Exp.Stim.Piston.Mat;
    runningTrial=(Exp.TrialEndT-Exp.TrialStartT)<4;
    I=runningTrial;
    II=find(I);
    mat(~I,:)=[];
    M=~isnan(mat(:,2:2:end));

    
    G=findgroups(array2table(M));
    samp=zeros(max(G),sampNperCat);
    
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
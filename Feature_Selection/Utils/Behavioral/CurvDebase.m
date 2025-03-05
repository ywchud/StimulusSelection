function [Whisker,Exp]=CurvDebase(Exp,Whisker)

for W=1:length(Whisker)
    fprintf('Finding baseline curvature for Whisker %d...\n',W)
 % shift the inherent curvature to zero (varaible to whisker angle)
    baselineBuffer=0.1;
    baselineEnd=floor((Exp.Stim.Piston.Delay-Exp.Cam.Delay-baselineBuffer)*Exp.videoFps );%determine baseline:camera starts to  before piston shoots out (0.1s as buffer)    
    % find BinEdges from histogram, skip this only if BinEdges are known
    if length(Exp.videoFps)==1
        BaseAngle=nan(Exp.TrN*baselineEnd,1);
        for i=1:Exp.TrN  
            if ~ismember(i,Exp.BadTrials)
                st=(Exp.TrN-1)*baselineEnd+1;
                try
                    BaseAngle(st:st+baselineEnd-1)=Whisker(W).trial(i).CenAngle.sig(1:baselineEnd);
                catch
                    fprintf('Not enough baseline. Adding tr %d as bad trial\n',i)
                    Exp.BadTrials(end+1)=i;Exp.BadTrials=unique(Exp.BadTrials);
                end
            end
        end
    else
        BaseAngle=nan(sum(baselineEnd),1);
        for i=1:Exp.TrN  
            if ~ismember(i,Exp.BadTrials)
                st=(Exp.TrN-1)*baselineEnd(i)+1;
                BaseAngle(st:st+baselineEnd(i)-1)=Whisker(W).trial(i).CenAngle.sig(1:baselineEnd(i));
            end
        end
    end
    

    figure
        h = histogram(BaseAngle,25);  %25 nbins
        BinEdges=h.BinEdges;
    close    

    %find all curvs with corresponding BinEdges during free whisk, i.e. before baselineEnd
    % BinEdges=0.45:0.05:0.95;
    CurvBinned=cell(length(BinEdges)); %bin no.
    for i=1:Exp.TrN  %trial no.
        if length(baselineEnd)==1
            baselineE=baselineEnd;
        else
            baselineE=baselineEnd(i);
        end
        if ~ismember(i,Exp.BadTrials)
            for n = 1:baselineE % frames no.    
                Angle=Whisker(W).trial(i).CenAngle.sig;                
                Curvature=Whisker(W).trial(i).Curvature;

                I=find(Angle(n)>BinEdges,1,'last');  %find correct bin
                if isempty(I),I=1;end
                temp=CurvBinned{I};temp=[temp;Curvature(n)];
                CurvBinned{I}=temp;
            end
        end
    end

    %find average baseline curvature for each bin
    CurvBase=zeros(1,length(CurvBinned));
    for i=1:length(CurvBinned)
            temp=CurvBinned{i};
            CurvBase(i)=mean(temp);
            if isnan(CurvBase(i))
                CurvBase(i)=0;
            end

    end

    %remove CurvBase from all frames
    for i=1:Exp.TrN  %trial no.
        if ~ismember(i,Exp.BadTrials)
        N=Whisker(W).trial(i).FrameN;
        C_Debase=nan(N,1);
        for n = 1:N  % frames no.
            Angle=Whisker(W).trial(i).CenAngle.sig; 
            Curvature=Whisker(W).trial(i).Curvature;

            I=find(Angle(n)>BinEdges,1,'last');  %find correct bin
            if isempty(I),I=1;end
            C_Debase(n)=Curvature(n)-CurvBase(I);
        end
        C_Debase=smoothdata(C_Debase,'sgolay',5);
        Whisker(W).trial(i).CurvDebase=C_Debase;
        fprintf('CurvDebase found for Whisker %d Trial %d\n',W,i)
        end
    end
    
    
    
end
end
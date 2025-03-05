function Exp=findResetTime(Exp,th)
    x=Exp.Sync;
    if isempty(th)
        th=0.8;
    end
    isSync=cell(Exp.TrN,1);
    ResetTime=cell(Exp.TrN,1);
    DriftTime=cell(Exp.TrN,1);
    
    
    for i=1:Exp.TrN
        if ~isempty(x{i})
            isSync{i}=logical(x{i}>th);
            xi=Signal( ' ',' ',isSync{i} , Exp.videoFps);
            ResetTime{i}=toColumn(xi.onset);      
            DriftTime{i}=toColumn(xi.offset);
        end
    end
    Exp.ResetTime=ResetTime;
    Exp.DriftTime=DriftTime;
    Exp.isSync=isSync;
end
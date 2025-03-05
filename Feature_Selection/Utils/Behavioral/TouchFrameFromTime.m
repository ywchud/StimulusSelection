function TF=TouchFrameFromTime(Exp,S)
%S can be from both Stim class or Neuron Class as long as they have
%corrtrial and touchtimes

    tt=S.TouchTimes;
    if isfield(S,'CorrTrials')
        ct=S.CorrTrials;
    elseif isfield(S,'CorrTrial')
        ct=S.CorrTrial;
    else
        error('Unknown S type')
    end
    
    if iscell(ct) && length(ct)==1
        ct=cell2mat(ct);
    end
    
    uct=unique(ct);



FT=Exp.FrameT;

TF=cell(Exp.TrN,1);
count=0;
for i=1:length(uct)   
    tr=uct(i);
    TF{tr}=find(ismember(FT{tr},tt));   
    count=count+length(TF{tr});
end
if count~=length(tt)
    error('Inconsistent touchN deduced from time to frames. Wrong Exp.Framet?')
end


end
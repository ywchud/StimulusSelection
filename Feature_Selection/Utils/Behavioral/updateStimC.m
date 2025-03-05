function [S_filtered,goodbools]=updateStimC(StimC,Exp,p)
% p.firstfew
% p.ITIth
% p.isrun
% p.Tdur
% p.omitTrial
% p.includeTrial

%%added 08/07/2022
% p.Toverlap
S_filtered=cell(size(StimC));
for i=1:length(StimC)
    if iscell(StimC(i))
        S=StimC{i};
    else
        S=StimC(i);
    end
L=length(S.TouchTimes);
goodbools=true(1,L);


if isfield(p,'firstfew') && ~isempty(p.firstfew)
    if ~isvector(p.firstfew) 
        error('p.firstfew must be a positive interger or vector')
    end
    goodbool=toColumn(FirstFewupdate(S,p.firstfew))';
    goodbools=[goodbools; goodbool];
end

if isfield(p,'ITIth') && (~isempty(p.ITIth) && ~isempty(S.TouchDeltaT))
    goodbool=toColumn(all(abs(S.TouchDeltaT)<=p.ITIth,2))';
    goodbools=[goodbools; goodbool];
end

if isfield(p,'isrun') && (~isempty(p.isrun) && ~isempty(S.TouchFrames))
    trF=[S.CorrTrials cell2mat(S.TouchFrames)];
    temp=false(1,size(trF,1));
    for m=1:size(trF,1)
        if p.isrun==1
            temp(m)=Exp.Runcmpers{trF(m,1)}(trF(m,2))>=3;
        elseif p.isrun==0
            temp(m)=Exp.Runcmpers{trF(m,1)}(trF(m,2))<=1;
        elseif p.isrun>0
            temp(m)=Exp.Runcmpers{trF(m,1)}(trF(m,2))>=p.isrun;
        end
    end
    goodbool=temp;
    goodbools=[goodbools; goodbool];
end


if isfield(p,'Tdur') && (~isempty(p.Tdur)  && ~isempty(S.TouchFrames))
    dur=S.ReleaseTimes-S.TouchTimes;
    if length(p.Tdur)==2
        goodbool=(dur>p.Tdur(1) & dur<p.Tdur(2));
    elseif length(p.Tdur)==1
        goodbool=(dur>p.Tdur);
    else
        error('uninterpretable p.Tdur')
    end
    goodbool=toColumn(goodbool)';
    goodbools=[goodbools; goodbool];
end
if isfield(p,'omitTrial') && ~isempty(p.omitTrial)
    goodbool=~ismember(S.CorrTrials,p.omitTrial);
    goodbool=toColumn(goodbool)';
    goodbools=[goodbools; goodbool];
end
if isfield(p,'includeTrial') && ~isempty(p.includeTrial)
    goodbool=ismember(S.CorrTrials,p.includeTrial);
    goodbool=toColumn(goodbool)';
    goodbools=[goodbools; goodbool];
end


%     vvvvv  overlap Touch vvv calculated and stored in S regardless of
%     filter
%     |/////|\\\\\\\\|///////|
%     x1    t1       t2      x2
% x1 first touch out of all W
% t1 last touch 
% t2 first release
% x2 last release
dur=S.ReleaseTimes-S.TouchTimes;
mat=[S.TouchDeltaT zeros(L,1)];
x1=min(mat,[],2);
t1=max(mat,[],2);
mat=[S.ReleaseDeltaT zeros(L,1)]+dur;
t2=min(mat,[],2);
x2=max(mat,[],2);  
overlap=(t2-t1)./(x2-x1);
S.Toverlap=overlap;
% if any(overlap<0)
%     123
% end
    
if isfield(p,'Toverlap') && (~isempty(p.Toverlap)  && ~isempty(S.TouchDeltaT))
    if length(p.Toverlap)==1
        goodbool=overlap>=p.Toverlap;
    elseif length(p.Toverlap)==2
        goodbool=overlap>=p.Toverlap(1) & overlap<p.Toverlap(2);
    else
        error('p.Toverlap should be singular of vector of 2')
    end

    %     figure
    %     histogram(overlap(all(goodbools,1)))
    goodbools=[goodbools; goodbool(:)'];
end

if isfield(p,'removeStim') && ~isempty(p.removeStim)
    ITIth=0.03; %touches within ITIth are removed
    
    SR=p.removeStim;
    if any(S.WhiskerBool~=SR.WhiskerBool)
        error('removeStim should have the same Whiskerbool as StimC')
    end
    goodbool=true(1,L);
    if ~isempty(S.TouchTimes) && ~isempty(SR.TouchTimes)
        td=S.TouchTimes-SR.TouchTimes';
        [td0,I]=min(abs(td),[],1);
        goodbool(I(td0<=ITIth))=0;
        fprintf('%s removed: n=%d\n',SR.type,sum(goodbool==0))
    elseif ~isempty(S.TouchTimes) && isempty(SR.TouchTimes)
        fprintf('%s removed: n=%d\n',SR.type,0)
    end
    goodbools=[goodbools; goodbool(:)'];
end



%final update by ccollective goodbools
S_filtered{i}=S.updateTouches(all(goodbools,1));

% test=all(goodbools,1);
% temp=find(overlap<0);
% if sum(test(overlap<0)==1)~=0
%     temp(test(overlap<0)==1)
% end

end
%if length(StimC)==1, output a StimC obj rather than a cell array
if length(S_filtered)==1
    S_filtered=S_filtered{1};
end
end
function goodbool=FirstFewupdate(StimC,corrtouch)
%corrtouch is an integer or range 
N=nan(1,2);
if isvector(corrtouch) && length(corrtouch)==1
    N(1)=corrtouch-1;
    N(2)=corrtouch-1;
elseif isvector(corrtouch) && length(corrtouch)==2
    N(1)=corrtouch(1)-1;
    N(2)=corrtouch(2)-1;
end

    G=StimC.CorrTrials;
    G_shift1=[nan; G(1:end-1)];
    sameasB4=(G==G_shift1);
    firstTouch=find(~sameasB4);
    

    temp=[firstTouch; length(G)+1];

    
    goodbool=false(1,length(G));
%     goodbool(firstTouch)=1;
    
    for i=1:length(temp)-1
        firstindex=min([temp(i)+N(1) temp(i+1)]);
        lastindex=min([temp(i)+N(2) temp(i+1)]);
        
        if firstindex==temp(i+1)  %firstindex is larger or equal to the next first touch, so skip
            continue
        elseif lastindex==temp(i+1) %constrain range to touches friom the same trial
            lastindex=lastindex-1;
        end
        goodbool(firstindex:lastindex)=1;
    end 
    
end
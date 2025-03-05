function Exp=getPistonInfo(Exp,Sigs)

%% Time info for Piston

PistonDig=Exp.Intan.PistonPort(Exp.Stim.Piston.active);
ActivePistonID=find(Exp.Stim.Piston.active);
TrialStartTime=Exp.TrialStartT;
TrialEndTime=Exp.TrialEndT;

PistonOut=[];PistonIn=[];Piston=[];
disp('Finding piston out/in times...')
for i=1:length(PistonDig)
    PistonOut{i}=Sigs.Piston(ActivePistonID(i)).onset/Exp.Fs;
    PistonIn{i}=Sigs.Piston(ActivePistonID(i)).offset/Exp.Fs;
    Piston(:,i*2-1:i*2)=nan(Exp.TrN,2);      %[P1out P1in P2out P2in ....]

end
%
% correction
for i=1:length(PistonDig)
    if PistonIn{i}(1)<PistonOut{i}(1)
        PistonIn{i}=PistonIn{i}(2:end);
        disp('Removed first pistonIn due to unfound pistonOut')
    end
    if length(PistonOut{i})~=length(PistonIn{i})
        PistonOut{i}=PistonOut{i}(1:length(PistonIn{i}));
        disp('Removed last pistonOut due to unfound pistonIn')
    end
    %stray 1s in contaminated/glitchy binary signal can cause misidentification of
    %piston out, here we ensure piston is at least out for 0.1s
    badOnsetSig=PistonIn{i}-PistonOut{i}<0.1;  
    if sum(badOnsetSig)>0
        sprintf('Stray 1s found in Piston %d signal, removing bad pistonOut instances',i)
        PistonIn{i}=PistonIn{i}(~badOnsetSig);
        PistonOut{i}=PistonOut{i}(~badOnsetSig);
    end
    
end
% create 2D matrix from al pistons' out & in
disp('Creating Piston matrix...')
BadTrials=[];
for i=1:Exp.TrN
    for j=1:length(PistonDig)
        k=find(PistonOut{j}> TrialStartTime(i) & PistonOut{j}<TrialEndTime(i));
        if ~isempty(k)
            if length(k)==1  %number of piston out in trial i, normally 1 unless hardware issues
                Piston(i,2*j-1)=PistonOut{j}(k);
                Piston(i,2*j)=PistonIn{j}(k);
            elseif length(k)>1
                fprintf('Trial %d has %d PistonOut, storing only the first piston In/Out, Badtrial assigned\n',i,length(k))
                BadTrials=[BadTrials;i];
                Piston(i,2*j-1)=PistonOut{j}(k(1));
                Piston(i,2*j)=PistonIn{j}(k(1));
            else
                error('Check code')
            end
        end
    end
end
% output
Exp.Stim.Piston.Delay=mean(Piston(:,1)-TrialStartTime,'omitnan');
% Exp.Stim.Piston.Delays=Piston(:,1)-TrialStartTime;
Exp.Stim.Piston.OffWin=mean(TrialEndTime-Piston(:,2),'omitnan');
% Exp.Stim.Piston.OffWins=TrialEndTime-Piston(:,2);

Exp.Stim.Piston.Out=PistonOut;
Exp.Stim.Piston.In=PistonIn;
Exp.Stim.Piston.Mat=Piston;
Exp.Stim.Piston.ID=PistonDig;
for j=1:length(PistonDig)
    Exp.Stim.Piston.Dur{1,j}=Exp.Stim.Piston.In{1,j}-Exp.Stim.Piston.Out{1,j};
end

if any(~Exp.Stim.Piston.active)
    Piston_padded=nan(size(Piston,1),size(Piston,2)+sum(~Exp.Stim.Piston.active)*2);
    PistonOut_padded=cell(1,Exp.Stim.Piston.tot);
    PistonIn_padded=cell(1,Exp.Stim.Piston.tot);
    PistonDig_padded=nan(1,Exp.Stim.Piston.tot);
    PistonDur_padded=cell(1,Exp.Stim.Piston.tot);
    j=1;
    for i=1:Exp.Stim.Piston.tot
        if Exp.Stim.Piston.active(i)
            Piston_padded(:,2*i-1:2*i)=Piston(:,2*j-1:2*j);
            PistonOut_padded(1,i)=PistonOut(1,j);
            PistonIn_padded(1,i)=PistonIn(1,j);
            PistonDig_padded(1,i)=PistonDig(1,j);
            PistonDur_padded(1,i)=Exp.Stim.Piston.Dur(1,j);
        else
            j=j-1;
        end
        j=j+1;
    end
    Exp.Stim.Piston.Mat=Piston_padded;
    Exp.Stim.Piston.Out=PistonOut_padded;
    Exp.Stim.Piston.In=PistonIn_padded;
    Exp.Stim.Piston.ID=PistonDig_padded;
    Exp.Stim.Piston.Dur=PistonDur_padded;
end
if isfield(Exp,'BadTrials')
    Exp.BadTrials=unique([Exp.BadTrials;BadTrials]);
else
    Exp.BadTrials=unique(BadTrials);
end
disp('Obtained Piston Delay,Out,In,Mat,BadTrials')
end
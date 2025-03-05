function Exp=TrialPerformace_discrimination(Exp)
%%  Get hit and CR  (go no-go)
% Exp.Lick.goW=2;
% Exp.Lick.nogoW=1;

if ~isfield(Exp.Lick,'Hit_bkup') 
    Exp.Lick.FA_bkup=Exp.Lick.FA;
    Exp.Lick.Miss_bkup=Exp.Lick.Miss;
    Exp.Lick.Hit_bkup=Exp.Lick.Hit;
%     Exp.Lick.CR_bkup=Exp.Lick.CR;
end

arr=zeros(1,Exp.Stim.Piston.tot);
arr(Exp.Lick.goW)=1;
temp=PistonComb(arr,Exp.Stim.Piston.Mat);
if all((Exp.Lick.Hit_bkup)==0) 
    Exp.Lick.Hit=(temp & ~Exp.Lick.Miss_bkup);
elseif all((Exp.Lick.Miss_bkup)==0)
    Exp.Lick.Miss=(temp & ~Exp.Lick.Hit_bkup);
end
%[Exp,Sigs]=DigAnalysis(Exp,Sigs,1);   %reset if messed up


arr=zeros(1,Exp.Stim.Piston.tot);
arr(Exp.Lick.nogoW)=1;
temp=PistonComb(arr,Exp.Stim.Piston.Mat);
Exp.Lick.CR=(temp & ~Exp.Lick.FA);

Exp.TrN
sum(Exp.Lick.Hit)+sum(Exp.Lick.Miss)+sum(Exp.Lick.FA)+sum(Exp.Lick.CR)
%-------allocate early lick trials
if ~isvector(Exp.PistonBuffer)
temp=Exp.PistonBuffer';
weirdTrialsbool=(all(isnan(temp)));  %not belonging to any category, e.g. no piston trials
tp_PistonOut=nan(Exp.TrN,1);
tp_PistonOut(~weirdTrialsbool)=temp(~isnan(temp));
test=nan(Exp.TrN,1);
for i=1:Exp.TrN
     temp=(Exp.Lick.Frames{i}-tp_PistonOut(i));
     test(i)=sum(temp<0);
end
else %no exact piston buffer yet, maybe cos no video analysis is done yet
    temp=Exp.Stim.Piston.Mat';
    weirdTrialsbool=(all(isnan(temp)));  %not belonging to any category, e.g. no piston trials
    tp_PistonOut=nan(Exp.TrN,1);
    temp2=temp(:,~weirdTrialsbool);
    temp2=temp2(~isnan(temp2));temp2=temp2(1:2:end);
    tp_PistonOut(~weirdTrialsbool)=temp2-Exp.TrialStartT(~weirdTrialsbool)+Exp.PistonBuffer;
    test=nan(Exp.TrN,1);
    for i=1:Exp.TrN
        temp=(Exp.Lick.Frames{i}-tp_PistonOut(i)*Exp.videoFps);
        test(i)=sum(temp<0);
    end
end
%----------------------------
if any(Exp.Lick.Miss_bkup+Exp.Lick.FA_bkup>1)  %newer arduino codes record early licks where both miss and FA=high
    Exp.Lick.earlyLick=(Exp.Lick.Miss_bkup+Exp.Lick.FA_bkup>1);
else
    Exp.Lick.earlyLick=(test>0);
end

Exp.Lick.Hit=(Exp.Lick.Hit & ~Exp.Lick.earlyLick);
Exp.Lick.Miss=(Exp.Lick.Miss & ~Exp.Lick.earlyLick);
Exp.Lick.FA=(Exp.Lick.FA & ~Exp.Lick.earlyLick);
Exp.Lick.CR=(Exp.Lick.CR & ~Exp.Lick.earlyLick);
Exp.Lick.weirdTrials=toColumn(weirdTrialsbool);

%sanity check
Exp.TrN
sum(Exp.Lick.Hit)+sum(Exp.Lick.Miss)+sum(Exp.Lick.FA)+sum(Exp.Lick.CR)+sum(Exp.Lick.earlyLick)+sum(Exp.Lick.weirdTrials)

end
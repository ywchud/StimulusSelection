function Whisker=findWhiskerEnvelope(Exp,Whisker,th,absMin,maxDur,TrsPksSamesize)

    if Exp.Stim.Piston.num~=length(Whisker)
        error('Expecting %d whiskers but found %d',Exp.Stim.Piston.num,length(Whisker))
    end
%     if isempty(th) && isempty(absMin) 
%         th=[];
%         absMin=0.03;
%         fprintf('th and absMin not found, assuming empty th and absMin= %d\n',absMin)
%     end
    if isempty(TrsPksSamesize)
         TrsPksSamesize=0;
     end
    disp('Obtaining whisker envelope')
    for j=1:Exp.TrN
        for i=1:Exp.Stim.Piston.num      
            if ~ismember(j,Exp.BadTrials)
%                 Whisker(i).trial(j)=Whisker(i).trial(j).findEnvelope(th,absMin,maxDur,TrsPksSamesize);
                WT=Whisker(i).trial(j);
                
                absMin=0;
                disp(j)
                [pks,trs,setPt,~,~] = findEnvelope(WT.CenAngle.sig,[],th,absMin,maxDur,TrsPksSamesize);
                amplitude=pks-trs;
                WT.Pks=pks; %correspond to most protracted(touch)
                WT.Trs=trs;
                WT.Setpoint=setPt;
                WT.Amplitude=amplitude;
                
                absMin=0;%0.03;
                [~,~,~,locsP,locsT] = findEnvelope(WT.CenAngle.sig,[],th,absMin,maxDur,TrsPksSamesize);
                WT.LocsP=locsP;
                WT.LocsT=locsT;
                
                Whisker(i).trial(j)=WT;
%                 fprintf('Obtained Envelope for Whisker %d Trial %d\n',i,j)
            end
        end
        try
        fprintf('Trial %d: %d %d %d %d\n',j,length(Whisker(1).trial(j).LocsP),length(Whisker(2).trial(j).LocsP),length(Whisker(3).trial(j).LocsP),length(Whisker(4).trial(j).LocsP))
        catch
        end
        end
    disp('Whisker envelope obtained') 
    
    

disp('Please check envelope, adjust absMin as needed') 
temp=[];
for i=1:Exp.TrN
    if ~ismember(i,Exp.BadTrials)
        temp=[temp;Whisker(1).trial(i).FrameN];
    else
        temp=[temp;nan];
    end
end
[~,MaxI]=max(temp,[],'omitnan');
[~,MinI]=min(temp,[],'omitnan');

for i=1:Exp.Stim.Piston.num
    for j=[MaxI MinI]
        figure

        plot(Whisker(i).trial(j).CenAngle.sig)
        hold on
        % [BS,~]=Whisker(1).trial(i).CenAngle.ButterSmooth(3, 100,4,0);
        % plot(BS)
        % plot(Whisker(1).trial(i).CenAngle.Phase(3, 100,4)/5)
        % plot(Whisker(1).trial(i).CenAngle.Amplitude(3, 100,4))

        pks=Whisker(i).trial(j).Pks;
        trs=Whisker(i).trial(j).Trs;
        locsP=Whisker(i).trial(j).LocsP;
        locsT=Whisker(i).trial(j).LocsT;


        plot(pks)
        plot(trs)
        scatter(locsP,pks(locsP))
        scatter(locsT,trs(locsT))
        title(sprintf('Whisker %d Trial %d',i,j))
% k=1400;
% xlim([k-150 k+150])
    end

end

end
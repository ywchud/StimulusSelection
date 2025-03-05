function plotIntanSig(Exp,Sigs,Trial2plot)
if exist('Trial2plot','var')
    factor=2;
else
    factor=100;
end
    disp('Plotting signals')
    figure
    plot(downsample(Sigs.Trial.sig*1.1,factor),'DisplayName','Trial')  %trial
    hold on
%     plot(downsample(Sigs.Camera.sig*1.05,factor),'DisplayName','Cam')  %trial
    try
     plot(downsample(Sigs.HitIntan.sig*1.4,factor),'DisplayName','Hit')  %trial
    catch
        plot(downsample(Sigs.MissIntan.sig*1.4,factor),'DisplayName','Miss')  %trial
    end
    % plot(downsample(DigInput(:,8)*1,factor)) 

    % plot(downsample(DigInput(:,5)*0.8,factor))
    % plot(downsample(DigInput(:,16)*0.7,factor)) 
    % plot(downsample(DigInput(:,1)*0.6,factor)) 
    % plot(downsample(DigInput(:,2)*0.5,factor)) 
    for i=1:length(Exp.Intan.PistonPort)
        plot(downsample(Sigs.Piston(i).sig*(0.9-0.1*i),factor),'DisplayName',sprintf('Piston %d',Sigs.Piston(i).Port))
    end

    try
        plot(downsample(Sigs.Spout.sig*1.3,factor),'DisplayName','Spout')
    catch
    end
    try
        plot(downsample(Sigs.Solenoid.sig*1.2,factor),'DisplayName','Solenoid')
    catch
    end
    
    if ~isempty(Sigs.Optic.sig)
        v=Sigs.Optic.sig;
        plot(downsample((v-min(v))/max(v-min(v))*0.4,factor),'r','DisplayName','Optic');
    end
    if isfield(Sigs,'LickPiezo')
    if ~isempty(Sigs.LickPiezo.sig)
        v=Sigs.LickPiezo.sig;
        plot(downsample((v-min(v))/max(v-min(v))*0.5,factor),'k','DisplayName','LickPiezo');
%         plot(downsample(v*0.5,factor),'k','DisplayName','LickPiezo');
    end
    end
    legend;
    if exist('Trial2plot','var')  %optional
        k=Trial2plot;
        xlim([Exp.TrialStart(k) Exp.TrialEnd(k)]/factor+[-100 100])
    end
end
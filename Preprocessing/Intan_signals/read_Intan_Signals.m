function [Exp,Sigs]=read_Intan_Signals(Exp)
%% Intan Signals
cd(Exp.codeDir)
% Define signal objects:
Sigs.Camera=Signal( Exp.Intan.CamPort,'binary', [], [] );
clear Piston
for i=1:length(Exp.Intan.PistonPort)
    Sigs.Piston(i)=Signal( Exp.Intan.PistonPort(i),'binary', [], [] );
end

Sigs.Trial=Signal( Exp.Intan.TrialPort,'binary', [], [] );
Sigs.Wheel=Signal( Exp.Intan.WheelPort,'binary', [], [] );

%Reward related
if isfield(Exp.Intan,'SpoutPort') && Exp.Intan.SpoutPort_active
    Sigs.Spout=Signal( Exp.Intan.SpoutPort,'binary', [], [] );
end
if isfield(Exp.Intan,'Solenoid')
    Sigs.Solenoid=Signal( Exp.Intan.Solenoid,'binary', [], [] );
end
if isfield(Exp.Intan,'MissPort')
    Sigs.MissIntan=Signal( Exp.Intan.MissPort,'binary', [], [] );
elseif isfield(Exp.Intan,'HitPort')
    Sigs.HitIntan=Signal( Exp.Intan.HitPort,'binary', [], [] );
end

Sigs.FAIntan=Signal( Exp.Intan.FAPort,'binary', [], [] );
Sigs.LickPiezo=Signal( nan,'binary', [], [] );

% Sigs.Test=Signal( Exp.Intan.test,'binary', [], [] );


% Analyse intan signals
read_Intan_RHD2000_file(Exp.Path.data);
clear reference_channel notes filename board_dig_in_channels
if ~exist('board_adc_channels','var')
    board_adc_channels=nan;
end
end
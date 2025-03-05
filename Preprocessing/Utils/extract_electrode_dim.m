function Exp=extract_electrode_dim(Exp)
%% Get LFP (do it early with empty workspace, need more RAM!)
% get channelmap % Elapsed time is 94.134363 seconds
% Exp.Path.data='E:\Darren\DC_DR97M3_210602_163418';
% Exp.Path.data='E:\Darren\DC_DR98M1_210525_161523';
% Exp.Path.data='E:\Darren\DC_DR97M4_210605_185611';
Exp.KSort.ChDepth=1000; %um
Exp.KSort.tipGap=25;
tic
try
    load(fullfile(Exp.Path.data,'rez2.mat'));  %only after KS ver 2.5 i think
catch
    load(fullfile(Exp.Path.data,'rez.mat'));
end
Exp.KSort.chanMap=rez.ops.chanMap;
Exp.KSort.y=rez.ycoords;
Exp.KSort.x=rez.xcoords;
clear rez
toc


if sum(Exp.KSort.chanMap==0)>1
    error('wat')
elseif sum(Exp.KSort.chanMap==0)==1
    Exp.KSort.chanMap=Exp.KSort.chanMap+1;  %make index starts with 1
    disp('ChanMap Index shifted to 1')
end


G=findgroups(Exp.KSort.x);
Exp.KSort.ShankN=max(G);
for i=1:Exp.KSort.ShankN
    Exp.KSort.ChN(i)=sum(G==i);
end
Exp.KSort.ChGap=abs(mode(diff(Exp.KSort.y)));
if Exp.KSort.ChGap==1  %prob a wrong channel map
    Exp.KSort.ChGap=20;
end
Exp.KSort.top=Exp.KSort.ChDepth-Exp.KSort.ChGap*(Exp.KSort.ChN-1);
Exp.KSort.bot=Exp.KSort.ChDepth-Exp.KSort.tipGap;
end
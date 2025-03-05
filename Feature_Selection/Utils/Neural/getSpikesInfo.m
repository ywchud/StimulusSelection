function [Exp,Neuron]= getSpikesInfo(Exp,removeTh)


Path=Exp.Path.data;
%The cluster ID of every spike
SpikeClusters = readNPY(fullfile(Path,'\spike_clusters.npy'));
%List of every spike in SAMPLES
SpikeSamples = readNPY(fullfile(Path,'\spike_times.npy'));
%     ChanMap = readNPY('E:\Darren\DarrenDualBC_Halo_R_190820_175141\B\channel_map.npy');
%     ChanPos = readNPY('E:\Darren\ScottSC_190628_143553\channel_positions.npy');

% Template = readNPY(fullfile(Path,'\templates.npy'));
% TemplateID = readNPY(fullfile(Path,'\templates_ind.npy'));
% Spike_templates = readNPY(fullfile(Path,'\spike_templates.npy'));
% Similiar_templates = readNPY(fullfile(Path,'\similar_templates.npy'));
% Templates_feat = readNPY(fullfile(Path,'\template_features.npy'));
% Template_feat_ind = readNPY(fullfile(Path,'\template_feature_ind.npy'));
ChannelMap = readNPY(fullfile(Path,'\channel_map.npy'));
ChannelPos = readNPY(fullfile(Path,'\channel_positions.npy'));

%Necessary variables (temporal)
% TrialStartTime = t(TrialStart);
% TrialEndTime = t(TrialEnd); 
% WheelTicTime = t(WheelTics);
% Trials = (1:length(TrialStart))';
% WheelSpeedTrialWin = NaN(2500,length(TrialStart));

% (spatial) %get spike shape/geo info. of each cluster
% try
ClusterInfo = tdfread(fullfile(Path,'\cluster_info.tsv')); 
ListSpkClust = ClusterInfo.id;
Depth = ClusterInfo.depth;  %may be flipped so not used

Ch=ClusterInfo.ch+1;
%remove noise cluster by tag
isNoise=logical(ClusterInfo.group(:,1)=='n');
isGood=logical(ClusterInfo.group(:,1)=='g');
fmt=['%d Noise cluster found and removed:\n' repmat(' %1.0f\n',1,numel(ListSpkClust(isNoise)))];
fprintf(fmt,sum(isNoise),ListSpkClust(isNoise))

if exist('removeTh','var') && removeTh>0
    spkcount=nan(length(ListSpkClust),1);
    for i=1:length(ListSpkClust)
        spkcount(i)=sum(SpikeClusters == ListSpkClust(i));
    end
    lowcount=spkcount<removeTh;
    fmt=['%d Low firing cluster found and removed:\n' repmat(' %1.0f\n',1,numel(ListSpkClust(lowcount)))];
    fprintf(fmt,sum(lowcount),ListSpkClust(lowcount))
    isNoise=isNoise | lowcount;
end


isGood=isGood(~isNoise);
ListSpkClust =ListSpkClust(~isNoise);
Depth=Depth(~isNoise);
Ch=Ch(~isNoise);

X=toColumn(Exp.KSort.x);
Y=toColumn(Exp.KSort.y);
[g,shankX]=findgroups(X);
dy=abs(Y(1)-Y(2));
if dy==1
    fprintf('dy between channel found as 1, probably using an old channel map, changing to 20um\n')
    dy=20;
    Y=Y*dy;  
end


try
    temp=abs(diff(shankX));
    dx=temp(1);
catch
    dx=0; %only 1 shank
end

ShankDepths=cell(1,max(g));
for i=1:length(ShankDepths)
    ShankDepths{i}=Y(g==i);
end

if ShankDepths{1}(1)<ShankDepths{1}(end)  %isflipped
    disp('y coor is ascending, flipping depths given by Ksort')
    for i=1:max(g)
        temp=ShankDepths{i};
        ShankDepths{i}=flip(temp);
        Y(g==i)=flip(temp);    
    end
end


XX=X(Ch);
YY=Y(Ch);
shank=findgroups(XX);
Depth=YY;
%add Exp.KSort.ChDepth to get actual depth
Depth=Depth-(Exp.KSort.ChDepth-Exp.KSort.tipGap);
Y=Y-Exp.KSort.ChDepth;
YY=Depth;

figure
scatter(X,Y,'r*') %probe geometry
hold on
G = findgroups(XX,YY);

for i=1:max(G)
    TXT='   ';
    groupi=find(G==i);
    for j=1:length(groupi)
        txt=sprintf('%s%d',setColor(isGood(groupi(j))),ListSpkClust(groupi(j)));
        TXT=append(TXT,'   ',txt);
    end
    
    
    text(XX(groupi(1)),YY(groupi(1)),TXT)
end
try
    xlim([min(XX)-dx max(XX)+dx])
catch
    xlim([min(XX)-1 max(XX)+1])
end
ylim([min(YY)-dy max(YY)+dy])
title('Cluster location (um)')
filename=fullfile(Exp.Path.save,'NeuronChannelPosition.jpg');
saveas(gcf,filename)
% Channels = ClusterInfo.channel;
% ClustersXY=[];
% ClustersXY(:,2) =  rem(Channels,8);    %y
% ClustersXY(:,1) = (Channels-ClustersXY(:,2))/8;%x 
% ClustersXY=ClustersXY+1;%start from x=1,y=1
% [XX,YY]=ClusterElectrodeMap(ClustersXY(:,1),ClustersXY(:,2),64);
% catch
%     error('cluster_info.tsv not found or incomplete/unmatch struct')
% %     ListSpkClust=unique(SpikeClusters);
% end
[Exp.si,Exp.sj]=subplotDim(length(ListSpkClust));
% 
% figure
% for i=1:length(XX)
%     hold on
%     scatter(XX,YY,'.')
%     c=num2str(ListSpkClust(i));
%     text(XX(i)+1,YY(i)+1,c);
% end
% figure
% for i=1:size(Template,1)
% for j=1:size(Template,3)
%     subplot(si,sj,j)
%     plot(Template(i,:,j))
% end
% pause(2)
% end
%%
%Output is cell matrix of spike ID x Trial, with each cell containing the
%spikes in seconds relative to zero by subtracting the trial start time
NeuronN=length(ListSpkClust);
Neuron=NeuronClass.empty;


for i = 1:NeuronN %Loops through every cluster   
    %Spike time in seconds as single column cell array, each cell contains
    %all spikes for the cluster
    try
        SpikeTime = Exp.t(SpikeSamples(SpikeClusters == ListSpkClust(i)));  
    catch
        temp=SpikeSamples(SpikeClusters == ListSpkClust(i));
        fprintf('Spike Index found larger than max Exp.t, removed. Neuron %d\n',i)
        SpikeTime = Exp.t(temp(temp<=length(Exp.t)));  
        
    end
    Neuron(i)=NeuronClass(Exp,i,ListSpkClust(i),SpikeTime,XX(i),Depth(i),shank(i),isGood(i));    %Depth=0 is cortex surface, negative is deeper
end





disp('done')
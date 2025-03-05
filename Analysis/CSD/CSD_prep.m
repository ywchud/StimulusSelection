function CSD_data=CSD_prep(Exp,LFP,PSTHrange)
LFP_fs=Exp.KSort.LFP_fs;
% PSTHrange = [-0.1 0.2];

%----------choose electrode
%optrode
% ch_map=1+[58 56 54 51 53 47 63 61 48 60 62 52 50 55 57 59 32 34 36 38 40 42 44 46 49 45 43 41 39 37 35 33 30 28 26 24 22 20 18 16 15 19 21 23 25 27 29 31 4 6 8 13 11 17 1 3 14 2 0 10 12 9 7 5];
% temp1=1:32;
% temp2=33:64;
% %128 3-shanks (top 2 bot)
% ch_map=1+[21 22 23 44 45 46 16 47 17 18 19 48 49 50 12 51 13 14 15 52 53 54 8 55 9 10 11 56 57 58 4 59 5 6 7 60 61 62 0 63 1 2 3 ...
%             87 104 86 85 84 103 102 101 91 100 91 98 88 99 98 97 95 96 94 93 92 32 33 34 28 35 29 30 31 36 37 38 24 39 25 26 27 40 41 42 20 43 ...
%             127 126 125 67 124 66 65 64 123 122 121 71 120 70 69 68 119 118 117 75 116 74 73 72 115 114 113 79 112 78 77 76 111 110 109 83 108 82 81 80 107 106 105];
% temp1=1:43;
% temp2=44:85;
% temp3=86:128;
ch_map=Exp.KSort.chanMap;

StimC=Exp.StimCi(1,1:Exp.MStim.SingleN);

% Create touch window
clear window windowSize
for i=1:length(StimC)
    StimCi=StimC{i};
%     switch i
%         case 1 
%             StimCi=StimC.W1;  %B1 whisker
%         case 2
%             StimCi=StimC.W2; %c1 whisker
%         case 3
%             StimCi=StimC.W3; %D1 whisker
%         case 4
%             StimCi=StimC.W4; %C2 whisker
%     end
window{i}=ceil([StimCi.TouchTimes+PSTHrange(1) StimCi.TouchTimes+PSTHrange(2)]*LFP_fs);
windowSize{i}=window{i}(:,2)-window{i}(:,1);
window{i}(:,2)=window{i}(:,2)-(windowSize{i}-mode(windowSize{i}))-1;  %make sure all touches have equal window size
fprintf('Window size is: %d\n',mode(windowSize{i}))
end
%----------choose StimC for each shank
shankN=Exp.KSort.ShankN;
% Output segmented mean LFP
clear CSD_data
CSD_data=cell(1,shankN);
prev=0;
for i=1:shankN   
    activeW=find(Exp.Stim.Piston.active);
    shankPW=find(activeW==Exp.KSort.shankPW(i));
%     fname = ['Shk_',num2str(i)];
%     CSD_data.(fname) = zeros( Exp.KSort.ChN(i),mode(windowSize{shankPW}));
    CSD_data{i} = zeros( Exp.KSort.ChN(i),mode(windowSize{shankPW}));

    %touch window for shank i
    W=window{shankPW};
    %channel used for shank i
    temp=[0 Exp.KSort.ChN];
    Ch_used=prev+1:sum(temp(1:i+1));
    prev=sum(temp(1:i+1));
    pos=ch_map(Ch_used);
    for j=1:size(W,1)
%         CSD_data.(fname)=CSD_data.(fname)+LFP(pos,W(j,1):W(j,2));
         CSD_data{i}=CSD_data{i}+LFP(pos,W(j,1):W(j,2));

    end
    CSD_data{i}=CSD_data{i}/size(W,1);
end

%check broken channel and **interpolate(not implemented yet)
figure 
for j=1:shankN
    subplot(1,shankN,j)
    c=colormap(summer);
    botdepth=Exp.KSort.bot;
    topDepth=Exp.KSort.top(j);
    dCh=topDepth:Exp.KSort.ChGap:botdepth;

for i=1:size(CSD_data{j},1)
    offset=-topDepth-Exp.KSort.ChGap*(i-1);
%     plot(CSD_data{j}(i,:)-mean(CSD_data{j}(i,:))+offset,'color',[0.3 0.3 0.3]);
    plot(CSD_data{j}(i,:)-CSD_data{j}(i,1)+offset,'color',[0.3 0.3 0.3]);
    
    hold on
end
ylim([-botdepth-200 0])


end
filename=fullfile(Exp.Path.save,'touchLFP.jpg');
saveas(gcf,filename)


end
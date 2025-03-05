
%% Example script for curvature feature occupancy control
% to ensure fair comparison of spike data between stimulus
% condition, feature variance(e.g. whisking/running variance) is reduced by selective sampling from the overlapping feature
% distributions
%% curvature data: 
% Ct is the mean of curvature in window 0-300ms, ExpN x StimCN
%% normed curvature debase of each experiment for home and away, histogram
Ct_norm=cell(size(Ct));
% x=nan(length([data.NeuronExp]),1);y=x;
for i=1:15
x=Ct{i,1};
y=Ct{i,3};
z=Ct{i,5};

%norm
mu=mean(x);
sd=std(x,0,1);
x=(x-mu)/sd;
y=(y-mu)/sd;
z=(z-mu)/sd;
Ct_norm{i,1}=x;
Ct_norm{i,3}=y;
Ct_norm{i,5}=z;

% sanity check, is data the right size
Nid=find([data.NeuronExp]==i,1,'first');
if length(x)~=sum(data.Is{Nid,1})
    disp('1')
elseif length(y)~=sum(data.Is{Nid,3})
    disp('3')
elseif length(z)~=sum(data.Is{Nid,5})
    disp('5')
end
end
%%
binSize=1;
edges=-6:binSize:6;
% edges=-0.1:binSize:0.1;
Ecount=zeros(length(edges)-1,15,3);

figure
for i=1:15
    subplot(4,4,i)
    x=Ct_norm{i,1};
    y=Ct_norm{i,3};
    z=Ct_norm{i,5};
%     x=Ct{i,1};
%     y=Ct{i,3};
%     z=Ct{i,5};
% [~,edges]=histcounts([x(:);y(:)],'BinWidth',0.25*10^-3);
[~,edges]=histcounts([x(:);y(:);z(:)],edges);
temp=histogram(x,edges);%,'Normalization','probability'
Ecount(:,i,1)=temp.Values;
hold on
temp=histogram(y,edges);
Ecount(:,i,2)=temp.Values;
temp=histogram(z,edges);
Ecount(:,i,3)=temp.Values;
axis square
% [~,p]=ttest2(x,y);
% title(sprintf('p:%d',p))
end
%% Resample the home and away distribution based on their minimum available trial number per binned curvature 
% such that all subsequent analysis results should be curvature independent

%hit and away go
PSTHrange=data.PSTHrange_master;
winT=PSTHrange<=0.07 & PSTHrange>=0; %transient
winS=PSTHrange<=0.3 & PSTHrange>0.07; %sustained
winA=PSTHrange<=0.3 & PSTHrange>=0; %all
 win=winA;
%  th=0;
binSize=1;
edges=-6:binSize:6;
IterMax=20;
muT_boot_HA=cell(15,2,IterMax);  %exp,home/away, iter
for Iter=1:IterMax
for i=1:15
x1=Ecount(:,i,1);
x2=Ecount(:,i,2);
xmin=min([x1(:) x2(:)],[],2);

%resampling
temp = [data.NeuronExp]==i;% & others;
Nid=find([data.NeuronExp]==i,1,'first');

k=1;
a=data.Is{Nid,k};
L=sum(a);%     L=mode(cellfun(@(x) size(x,1),data.MUs(temp,win,1)));
TrI_resamp=false(L,1);
for j=1:length(edges)-1
    N=xmin(j);
    TrI=find(Ct_norm{i,k}>=edges(j) & Ct_norm{i,k}<=edges(j)+binSize);
    if ~isempty(TrI)
        TrI_resamp(randsample(TrI,N))=1;
    end
end
temp1=cellfun(@(x) mean(x(TrI_resamp)),data_MUs_normed(temp,win,k));%data_MUs_normed
y1=mean(temp1,2);muT_boot_HA{i,1,Iter}=y1;

k=3;
a=data.Is{Nid,k};
L=sum(a);
TrI_resamp=false(L,1);
for j=1:length(edges)-1
    N=xmin(j);
    TrI=find(Ct_norm{i,k}>=edges(j) & Ct_norm{i,k}<=edges(j)+binSize);
    if ~isempty(TrI)
        TrI_resamp(randsample(TrI,N))=1;
    end
end
temp1=cellfun(@(x) mean(x(TrI_resamp)),data_MUs_normed(temp,win,k));%data_MUs_normed
y2=mean(temp1,2);muT_boot_HA{i,2,Iter}=y2;


end
disp(Iter)
end
%% hit and miss
PSTHrange=data.PSTHrange_master;
winT=PSTHrange<=0.07 & PSTHrange>=0; %transient
winS=PSTHrange<=0.3 & PSTHrange>0.07; %sustained
winA=PSTHrange<=0.3 & PSTHrange>=0; %all
 win=winA;
%  th=0;
binSize=1;
edges=-6:binSize:6;
IterMax=20;
muT_boot_HM=cell(15,2,IterMax);  %exp,home/away, iter
for Iter=1:IterMax
for i=1:15
x1=Ecount(:,i,1);
x2=Ecount(:,i,3);
xmin=min([x1(:) x2(:)],[],2);

%resampling
temp = [data.NeuronExp]==i;% & others;
Nid=find([data.NeuronExp]==i,1,'first');

k=1;
a=data.Is{Nid,k};
L=sum(a);%     L=mode(cellfun(@(x) size(x,1),data.MUs(temp,win,1)));
TrI_resamp=false(L,1);
for j=1:length(edges)-1
    N=xmin(j);
    TrI=find(Ct_norm{i,k}>=edges(j) & Ct_norm{i,k}<=edges(j)+binSize);
    if ~isempty(TrI)
        TrI_resamp(randsample(TrI,N))=1;
    end
end
temp1=cellfun(@(x) mean(x(TrI_resamp)),data_MUs_normed(temp,win,k));%data_MUs_normed
y1=mean(temp1,2);muT_boot_HM{i,1,Iter}=y1;

k=5;
% a=data.Is{Nid,k};
% L=sum(a);
L=mode(cellfun(@(x) size(x,1),data.MUs(Nid,win,k)));
TrI_resamp=false(L,1);
for j=1:length(edges)-1
    N=xmin(j);
    TrI=find(Ct_norm{i,k}>=edges(j) & Ct_norm{i,k}<=edges(j)+binSize);
    if ~isempty(TrI) && L~=0
        TrI_resamp(randsample(TrI,N))=1;
    end
end
temp1=cellfun(@(x) mean(x(TrI_resamp)),data_MUs_normed(temp,win,k));%data_MUs_normed
y2=mean(temp1,2);muT_boot_HM{i,2,Iter}=y2;

end
disp(Iter)
end
%%
%% structure into delta(zfr)
%------
Y=cell(IterMax,2);  %iterN,type(HA,HM)
for i=1:IterMax

    x1=cell2mat(muT_boot_HA(:,1,i,1));
    x2=cell2mat(muT_boot_HA(:,2,i,1));
    temp=(x2-x1);
%     temp=(x2-x1)./(x2+x1);
    Y{i,1}=temp(:);

    x1=cell2mat(muT_boot_HM(:,1,i,1));
    x2=cell2mat(muT_boot_HM(:,2,i,1));
    temp=(x2-x1);
%     temp=(x2-x1)./(x2+x1);
    Y{i,2}=temp(:);
    
end
%------
Yall=cell(2,1);  %iterN,type(HA,HM)
%     x1=cellfun(@(x) mean(x),MUmat1(:,win1));%data_MUs_normed
%     x3=cellfun(@(x) mean(x),MUmat2(:,win1));%data_MUs_normed
%     x5=cellfun(@(x) mean(x),data_MUs_normed(:,win1,5));%data_MUs_normed
%     Yall{1}=x3-x1;
%     Yall{2}=x5-x1;
%     
    x1=getWinMean(squeeze(data_MUs_normed(:,:,1)),win(1:end-1),1);
    x3=getWinMean(squeeze(data_MUs_normed(:,:,3)),win(1:end-1),1);
    x5=getWinMean(squeeze(data_MUs_normed(:,:,5)),win(1:end-1),1);
    Yall{1}=x3-x1;
    Yall{2}=x5-x1;
%% ------plot
figure
% binSize=0.1;
% edges=-1:binSize:1;
binSize=1;
edges=-10:binSize:10;
clear yy p
for k=1:4
%     others=MI(:,3)>=0;
others=MI(:,3)<0;

switch k
    case 1
        type=1;
        temp=~data.Nsig & ismember([data.NeuronExp],find(data.area_SC)) & ismember([data.NeuronExp],find(data.area_context)) & others;
        s='away-home dzfr SC';
    case 2
        type=1;
        temp=~data.Nsig & ismember([data.NeuronExp],find(~data.area_SC)) & ismember([data.NeuronExp],find(data.area_context)) & others;
        s='away-home dzfr S1';
    case 3
        type=2;
        temp=~data.Nsig & ismember([data.NeuronExp],find(data.area_SC)) & ~nomiss & others;
        s='miss-home dzfr SC';
    case 4
        type=2;
        temp=~data.Nsig & ismember([data.NeuronExp],find(~data.area_SC)) & ~nomiss & others;      
        s='miss-home dzfr S1';
end
subplot(2,2,k)

yy=zeros(IterMax,length(edges)-1);
for i=1:IterMax
y=histcounts(Y{i,type}(temp),edges,'Normalization','probability'); %y = discretize(mu1(temp,1),edges)
h=plot(edges(1:end-1)+diff(edges)/2,y,'-');
h.Color=[0.6 0.6 0.6 0.3];
yy(i,:)=y;
[~,p(i)]=ttest(Y{i,type}(temp));
hold on
end
axis square,xlabel(s),ylabel('Neurons%*100')
% h=plot(edges(1:end-1)+diff(edges)/2,mean(yy,1),'-');

%--------takes mean of each MI bin's IterN counts
c=[0 0 1];
plotWE(edges(1:end-1)+diff(edges)/2,mean(yy,1,'omitnan'),std(yy,0,1,'omitnan')./size(yy,1)^0.5,c,sprintf('n: %d',size(yy,1)));

%--------takes mean of each neuron's IterN MI
YY=mean(cell2mat(Y(:,type)'),2,'omitnan');
y=histcounts(YY(temp),edges,'Normalization','probability'); %y = discretize(mu1(temp,1),edges)
h=plot(edges(1:end-1)+diff(edges)/2,y,'-');
[~,p1]=ttest(YY(temp));
% h=plot(edges(1:end-1)+diff(edges)/2,mean(yy,1),'-');
h.Color=[0 0 0];
h.Color=[1 0 0];
xline(0)
%-------raw all trials without ctrl
y=histcounts(Yall{type}(temp),edges,'Normalization','probability'); %y = discretize(mu1(temp,1),edges)
h=plot(edges(1:end-1)+diff(edges)/2,y,'k-');
[~,p2]=ttest(Yall{type}(temp));

[~,p3]=ttest(YY(temp),Yall{type}(temp));

title(sprintf('pCtrl:%d pAll:%d pair%d',p1,p2,p3))
legend({sprintf('n:%d',sum(temp)),sprintf('mu_{ctrl}:%d',mean(YY(temp))),sprintf('mu_{all}:%d',mean(Yall{type}(temp)))})
axis square,ylim([0 0.6])
end
%%
%%---
%% how MI home vs away change with curvature (per mouse)
% ----screen good home vs away
win=winA;

clear c
% c(5:8,:)=createColorMap(length(5:8),[0.5 0.5 1;0 0 0]);
figure(20),clf
figure(21),clf
for i=1:5
MI_temp=[];

MI_C=[];
figure(20)
% win=winT;
% win=winS;
others=MI(:,3)>=0;
% temp = ~data.Nsig & ismember([data.NeuronExp],find(data.area_SC)) & ismember([data.NeuronExp],find(data.area_context)) & others;
temp = ~data.Nsig & [data.NeuronExp]==i & others;
L=sum(BinI(i,:));
c=createColorMap(L,[0 0 1;1 0 0]);
for j=1:L
k=find(BinI(i,:));
k=k(j);
binSize=1;
    edges=-6:binSize:6;
TrI=Ct_norm{i,1}>=edges(k) & Ct_norm{i,1}<=edges(k)+binSize;
x1=cellfun(@(x) mean(x(TrI)),data.MUs(temp,win,1));
x1=mean(x1,2);
TrIa=Ct_norm{i,3}>=edges(k) & Ct_norm{i,3}<=edges(k)+binSize;
x3=cellfun(@(x) mean(x(TrIa)),data.MUs(temp,win,3));
x3=mean(x3,2);



MI_temp=(x3-x1)./(x1+x3);
subplot(2,3,i)
binSize=0.2;
edges=-1:binSize:1;
y=histcounts(MI_temp,edges,'Normalization','probability'); %y = discretize(mu1(temp,1),edges)
h=plot(edges(1:end-1)+diff(edges)/2,y,'-');
h.Color=c(j,:);
h.DisplayName=sprintf('n=%d TrHome=%d TrAway=%d',length(MI_temp),sum(TrI),sum(TrIa));
hold on
axis square,xlabel('MI go away-home'),ylabel('Neurons%*100'),title('SC hit-pref')
legend
if j==1
    MI_C(:,1)=MI_temp;
elseif j==L
    MI_C(:,2)=MI_temp;
end
end
xline(0)

figure(21)
subplot(2,3,i)
scatter(MI_C(:,1),MI_C(:,2))
xlabel('Low curv'),ylabel('High curv'),title('MI go away-home')
[~,p]=ttest(MI_C(:,1),MI_C(:,2));
title(sprintf('p=%d',p)),axis square,refline(1,0)

end
%% ----------
figure
for i=1:5

% win=winA;
% win=winT;
% win=winS;
others=MI(:,3)>=0;
% temp = ~data.Nsig & ismember([data.NeuronExp],find(data.area_SC)) & ismember([data.NeuronExp],find(data.area_context)) & others;
temp = ~data.Nsig & [data.NeuronExp]==i & others;
L=sum(BinI(i,:));
c=createColorMap(L,[0 0 1;1 0 0]);
for j=1:L
k=find(BinI(i,:));
k=k(j);
binSize=1;
    edges=-6:binSize:6;
TrI=Ct_norm{i,1}>=edges(k) & Ct_norm{i,1}<=edges(k)+binSize;
x1=cellfun(@(x) mean(x(TrI)),data_MUs_normed(temp,win,1));
x1=mean(x1,2);
TrIa=Ct_norm{i,3}>=edges(k) & Ct_norm{i,3}<=edges(k)+binSize;
x3=cellfun(@(x) mean(x(TrIa)),data_MUs_normed(temp,win,3));
x3=mean(x3,2);
subplot(5,4,j+4*(i-1))
scatter(x1,x3,'.'),hold on
axis square,xlabel('Home z firing'),ylabel('Away z firing'),title(sprintf('Curv strength:%d',j))
xlim([0 20])
ylim([0 20])
refline(1,0)

end
end

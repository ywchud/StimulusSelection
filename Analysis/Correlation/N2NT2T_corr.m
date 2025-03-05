function [Results,lags,corrNeuPairs]=N2NT2T_corr(data,data_MUs_normed,H,MI)
%% 
PSTHrange=data.PSTHrange_master;
winT=PSTHrange<=0.07 & PSTHrange>=0; %transient
winS=PSTHrange<=0.3 & PSTHrange>0.07; %sustained
winA=PSTHrange<=0.3 & PSTHrange>=0; %all
winB=PSTHrange>=-1.5 & PSTHrange<-0.05; %all
winSteps=0.01;
winRange=[-0.2 0.3];
range=PSTHrange>=winRange(1) & PSTHrange<winRange(2);
winSize=sum(range);
t=data.PSTHrange_master(range)+winSteps;
lags=(-(sum(range)-1):sum(range)-1)*winSteps;
exps=unique(data.NeuronExp);

N=length(data.NeuronExp);
expid=[data.NeuronExp];
expN=max(unique([data.NeuronExp]));

H.Hit_facilitated(isnan(H.Hit_facilitated))=0;
G1=~data.Nsig & H.Hit_facilitated & MI(:,3)>=0;%
H.Hit_suppressed(isnan(H.Hit_suppressed))=0;
G2=~data.Nsig & H.Hit_suppressed & MI(:,3)<0;%
H.CR_facilitated(isnan(H.CR_facilitated))=0;
G3=~data.Nsig & H.CR_facilitated & MI(:,3)<0;%
H.CR_suppressed(isnan(H.CR_suppressed))=0;
G4=~data.Nsig & H.CR_suppressed & MI(:,3)>=0;%

Results=cell(expN,2);
for Exp=1:expN
for G=1:2
    switch G
        case 1
            %--------GO-------------------------------
            stimType=1;
            Neu_i=find(G1 & [data.NeuronExp]==Exp);Li=length(Neu_i);
            Neu_j=find(G2 & [data.NeuronExp]==Exp);Lj=length(Neu_j);
        case 2
            %--------NG-------------------------------
            stimType=2;
            Neu_i=find(G3 & [data.NeuronExp]==Exp);Li=length(Neu_i);
            Neu_j=find(G4 & [data.NeuronExp]==Exp);Lj=length(Neu_j);   
    end
Result=nan(Li*Lj,winSize*2-1); 
corrNeuPair=nan(Li*Lj,2); 
for i=1:Li
    for j=1:Lj
        k=(i-1)*Lj+j;
%         i_Exp=expid(i);
%         j_Exp=expid(j);
%         x=squeeze(data.MUs(Neu_i(i),range,stimType));
        x=squeeze(data_MUs_normed(Neu_i(i),range,stimType));
        x1=cell2mat(x);
        [TrN1,~] = size(x1);
        
%         x=squeeze(data.MUs(Neu_j(j),range,stimType));
        x=squeeze(data_MUs_normed(Neu_j(j),range,stimType));
        x2=cell2mat(x);
        
        corrNeuPair(k,:)=[Neu_i(i) Neu_j(j)];
%         [TrN2,~] = size(x2);
        
%         ------trial X trial (slow)
%         correlation_result=zeros(TrN1,winSize*2-1);
%         for m=1:TrN1
%                 n=m;
% %                 correlation_result(m,:)=xcorr(x1(m,:),x2(n,:),[],'normalized');
%                 correlation_result(m,:)=xcorr(x1(m,:),x2(m,:),[],'normalized');
% %                 correlation_result(m,:)=xcorr(x1(m,:),x2(n,:));
%         end
%         Result(k,:)=mean(correlation_result,1,'omitnan');
%         %---plot and check
%             figure(111)
%             subplot(1,2,1),hold off
% %             plot(PSTHrange(range),mean(x1,1,'omitnan'))
% %             hold on
% %             plot(PSTHrange(range),mean(x2,1,'omitnan'))
%             plot(PSTHrange(range),x1(129,:))
%             hold on
%             plot(PSTHrange(range),x2(129,:))
%             subplot(1,2,2),hold off
% %             plot(lags,Result(k,:))
%             plot(lags,correlation_result(129,:))
% %             y=correlation_result;
% %             plotWE(lags,mean(y,1,'omitnan'),std(y,0,1,'omitnan')./size(y,1)^0.5,[0 0 0],sprintf('n: %d',size(y,1)));
%             figure(112),hold off
%             plot(lags,correlation_result')
%             
%             title(sprintf('%d %d',i,j))
%         pause(2)
         %------trialavg X trialavg (faster)
        Result(k,:)=xcorr(mean(x1,1,'omitnan'),mean(x2,1,'omitnan'),[],'normalized');        
        fprintf('Exp%dG%d: %d/%d\n',Exp,G,k,Li*Lj)
    end
end
Results{Exp,G}=Result;
corrNeuPairs{Exp,G}=corrNeuPair;
end
end

end
function I=findLag(A)
% center of mass of only the negative component
A(A>0)=0;
B = cumsum(A);
I=find(B<=sum(A)/2,1,'first');

%min after smooth
% A = smoothdata(A,"sgolay");
% [~,I]=min(A);
end
function [group1,group2]=extractTouchSpkJitter(PSTHmat,p)
PSTHrange=p.PSTHrange;
windowtime=[0 0.04];%p.windowtime;

zeroIndex=find(PSTHrange>=windowtime(1),1,'first');
rightIndex=find(PSTHrange>=windowtime(2),1,'first');
leftIndex=find(PSTHrange>=-windowtime(2),1,'first');

group1=[];
group2=[];
for i=1:size(PSTHmat,1)
    dfrom0=find(PSTHmat(i,:))-zeroIndex;
    df0=dfrom0(dfrom0<rightIndex-zeroIndex  & dfrom0>=leftIndex-zeroIndex);  
    group1=[group1,df0(df0>=0)+0.5];
    group2=[group2,abs(df0(df0<0)+0.5)];
end

% histogram(group1,0:1:rightIndex-zeroIndex,'Normalization','probability')
% hold on
% histogram(group2,0:1:rightIndex-zeroIndex,'Normalization','probability')
% 


end
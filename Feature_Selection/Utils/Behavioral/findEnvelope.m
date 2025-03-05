function [pks,trs,setPt,locsP,locsT] = findEnvelope(data,midpoint,th,absMin,maxDur,TrsPksSamesize)
% Darren
% finds upper envelope (pks) and lower envelope (trs) of signal (data) by
% akimai spline fitting of peaks
% peaks or troughs further than (midpoint) are filtered out
% no filter if left empty
% peak2peak smaller than (0<th<1 of population) are filtered out
%or absMin, whichever is larger
%maxDur determines how continuous small bumps are filtered, if a good
%absMin is set, no limit i.e. maxDur =[] is preferred so all bumps smaller
%than absMin are filtered. Set maxDur=1 only if some small bumps are
%necessary for nice envelope visualisation
%TrsPksSamesize if=1, locsP and locsT have the same array size with locsT
%preceding locsP

y=data;
if length(midpoint)==length(data) && ~isempty(midpoint)
    z=midpoint;
else
    z=NaN(1,length(data))';
end

NumZeros=100;
interPdegree=1;
if sum(isnan(data))==0
    %reduce edge effect
    y=[ones(NumZeros,1)*mean(data(1:round(length(data)/6))); y; ones(NumZeros,1)*mean(data(round(length(data)/6*5):end))];
    z=[ones(NumZeros,1)*mean(data(1:round(length(data)/6))); z; ones(NumZeros,1)*mean(data(round(length(data)/6*5):end))];
    %find peaks and troughs
    [pks,locsP] = findpeaks(y);
    inv_y=max(y)-y;
    [trs,locsT] = findpeaks(inv_y);
    trs=(max(y)-trs);
    %filter small bumps (minD)   
    locs=sort([locsP; locsT]);
    d=[abs(diff(y(locs))); 1];
    D=sort(d(1:end-1));   
%     D=D(D>=5*10^-5);
    
    if isempty(th) || th<=0
        minD=0;          
    else
        try
            minD=D(round(length(D)*th));   
        catch
            minD=0;
        end
    end
    
    minD=max([minD absMin]);  %absolute value absMin(not %) if left empty
    
    temp=(d<=minD);
    [startf, endf]=segmentWminDur(temp',1,maxDur);  
    
%     figure
%     plot(data)
%     hold on 
%     scatter(locs-NumZeros,data(locs-NumZeros),'c')
%     scatter(locs(temp)-NumZeros,data(locs(temp)-NumZeros),'r')
    
    for i=1:length(startf)
        if rem(endf(i)-startf(i),2)==0
            for j=startf(i):(endf(i)+1)
                locs(j)=nan;
            end
        else
            if ismember(locs(startf(i)),locsP)
                [~,I] = max(y(locs(startf(i):endf(i)+1)));
                for j=startf(i):endf(i)+1
                    if ~ismember(j,I-1+startf(i))
                        locs(j)=nan;
                    end
                end
            elseif ismember(locs(startf(i)),locsT)
                [~,I] = min(y(locs(startf(i):endf(i)+1)));
                for j=startf(i):endf(i)+1
                    if ~ismember(j,I-1+startf(i))
                        locs(j)=nan;
                    end
                end
            end                        
        end        
    end

%     figure
%     plot(data)
%     hold on 
%     scatter(rmmissing(locs-NumZeros),data(rmmissing(locs-NumZeros)),'c')

    
    
    
    L=locs(~isnan(locs));
    locsP=locsP(ismember(locsP,L));
    locsT=locsT(ismember(locsT,L));
    %akima spline fitting
    xq = toColumn(1:1/interPdegree:length(y));
    locs=[1; locsP(y(locsP)>z(locsP) | isnan(z(locsP)) ); length(y)];  %peak threshold set as given midpoint, ignore if empty
    pks=y(locs);
%     pks=akimai(locs,pks,xq);    
    pks = interp1(locs,pks,xq,'makima');
    pks = resample(pks,1,length(xq)/length(y));
    pks = pks(NumZeros+1:length(pks)-NumZeros);
    locsP = locs(locs>= NumZeros+1 & locs<=length(y)-NumZeros)-NumZeros;
    %     xq=[0:1/length(locs)/4:1]*length(y);
    %     trs = interp1(locs,trghs,xq,'spline');
    locs=[1; locsT(y(locsT)<z(locsT) | isnan(z(locsT)) ); length(y)];
    trs=y(locs);
    trs = interp1(locs,trs,xq,'makima');
%     trs=akimai(locs,trs,xq);
    trs = resample(trs,1,length(xq)/length(y));
    trs = trs(NumZeros+1:length(trs)-NumZeros);
    locsT = locs(locs>= NumZeros+1 & locs<=length(y)-NumZeros)-NumZeros;
%     

    if TrsPksSamesize
        LT=locsT;  %protraction: location of troughs
        LP=locsP;  %retraction: location of peaks
        if(LP(1)<LT(1))
            LP=LP(2:end);%('oops protraction(troughs) should come first');
            
        end
        Pnum=length(LP);
        try
            LT=LT(1:Pnum); %ensure T-P pairs
        catch
            disp('Pks than trs number diffrence larger than 1, likely bad signal')
            for k=1:length(LT)
                lp(k)=LP(find(LP>LT(k),1,'first'));
            end
            LP=lp;
        end
        locsT=LT;
        locsP=LP;
    end
    
    
else
    pks=0;
    trs=0;
end

% pks=max([pks data],[],2);   %ensures envelope wraps outside of data range
% trs=min([trs data],[],2);
setPt=(pks+trs)/2;

%re-evaluate if there is sudden setpoint change




end
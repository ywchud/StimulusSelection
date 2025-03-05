function [Iidx,deltaT]=findClosest(T1s,T2s,th)
% find the closest value T2 for each T1 and output an array of length(T1s)
% of the respective T2 pairIndexes (Iidx) and value difference (deltaT), 
% otherwise nan if abs(dt)>th

%T1s and T2s are array of same/different lengths, th is a double

  
    deltaT=nan(length(T1s),1);
    Iidx=nan(length(T1s),1);
    hasTh=~isempty(th);
    for i=1:length(T1s)
%         dt=T2s-T1s(i);        %too slow
%         [~,I] = min(abs(dt));
%         dt=dt(I);  
        t1=T1s(i);
        IL=find(T2s>=t1,1,'first');
        IS=find(T2s<=t1,1,'last');
        if ~isempty(IL) && isempty(IS)
            I=IL;
        elseif  isempty(IL) && ~isempty(IS)
            I=IS;
        elseif T2s(IL)-t1>t1-T2s(IS)
            I=IS;
        else
            I=IL;
        end
        dt=T2s(I)-t1;
        if ~isempty(dt)           
            if hasTh && abs(dt)<=th || ~hasTh
                Iidx(i) = I;
                deltaT(i)=dt;
            end
        end
    end
end
function BOOL = isbetween(signalT,StartT,EndT)   %finds if a data time point lies within the given combinations of StartTs and EndTs
    BOOL=zeros(length(signalT),1);
    for j=1:length(signalT)
        for i=1:length(StartT)
            if signalT(j)>=StartT(i) && signalT(j)<=EndT(i)
                BOOL(j)=1;break;
            end
        end
    end
end

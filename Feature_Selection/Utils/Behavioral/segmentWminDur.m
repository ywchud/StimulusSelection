function [startf, endf]=segmentWminDur(x,minDur,maxDur)
% Darren 06/18/19
% given bool array x of 0 and 1, it finds the start and end index(es) of a
% series of 1 with minimum duration, mindur and max duration, maxDur, if
% left empty, assume inf
if x(1)==0
    y=[0 diff(x) -1];
else
    y=[1 diff(x) -1];
end
startf=find(y==1);
endf=find(y==-1)-1;

while length(endf)>length(startf)
    endf=endf(1:end-1);
end

lengthf=endf-startf+1;

startf=startf(lengthf>=max([1 minDur]) & lengthf<=min([length(x) maxDur]));
endf=endf(lengthf>=max([1 minDur]) & lengthf<=min([length(x) maxDur]));
% endf=endf(lengthf>=minDur & lengthf<=maxDur);


end
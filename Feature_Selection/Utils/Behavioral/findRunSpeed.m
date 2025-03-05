function RS=findRunSpeed(Exp,RSsig)
%resample runspeed signal into frame resolution based on trial videos


RS=cell(Exp.TrN,1);
try
RS500 = resample(RSsig.sig,Exp.videoFps,RSsig.fs);
Speedbin=0:1/Exp.videoFps:max(Exp.t);

catch
    RS500 = resample(RSsig.sig,500,RSsig.fs);
    Speedbin=0:1/500:max(Exp.t);
    fprintf('Using 500fps for all vids')
end
for i=1:Exp.TrN
    [Iidx,deltaT]=findClosest(Exp.FrameT{i}(1),Speedbin,[]);
	fprintf('Trial %d: error sec %d\n',i,deltaT)
    if ~isnan(Iidx)
        RS{i}=RS500(Iidx:Iidx+Exp.FrameN(i)-1);
    end
end

end
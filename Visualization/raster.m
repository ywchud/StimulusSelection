function [data,emptyTouch]=raster(SpikeTime,timepoints,PSTHrange,color,py)
px=zeros(1,2);
if isempty(py)
    py=[]; %keep it empty
end
spkT=toColumn(SpikeTime);
px(1)=PSTHrange(1);px(2)=PSTHrange(end);
data=cell(1,length(timepoints));
for i=1:length(timepoints)
    data{i}=spkT(timepoints(i)+px(1)<spkT & timepoints(i)+px(2)>spkT)-timepoints(i);
end
emptyTouch=rasterplot(data,color,py,[PSTHrange(1) PSTHrange(end)]);
end
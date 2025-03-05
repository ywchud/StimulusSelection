function EmptyTrial=rasterplot(data,color,py,px)
    
%Data dimension: 1D cell array of trials,  each cell with 1D timepoints
%array for raster plot

%output:
% EmptyTrial gives the indexes of trials that are empty with nth plotted


if ~iscell(data)
    error('Expecting a cell data format')
end

if ~isempty(color) && size(color,1)==length(data)
    C=color;
elseif ~isempty(color) && size(color,1)==1
    C=repmat(color,length(data),1);
else
    C=repmat([0 0 0],length(data),1);
end

%determine plot scale
if isempty(py)
    py=[1 length(data)];
elseif length(py)~=2
    error('py must be a vector of 2')
elseif py(2)-py(1)+1~=length(data)
    error('py must be a vector of 2 binding the same size as data')
end
if isempty(px)
    maxD=0;
    minD=0;
    for i=1:length(data)
        maxD=max([toColumn(data{i});maxD]);
        minD=min([toColumn(data{i});minD]);
    end
    px=[minD maxD];
end
%scatter
EmptyTrial=[];
for i=1:length(data)
    hold on
    if ~isempty(data{i})
        h=plot(data{i},(py(1)+i-1)*ones(1,length(data{i})),'Color',C(i,:),'Marker','.','LineStyle','none');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h.MarkerSize=2;
    else
        EmptyTrial=[EmptyTrial;i];
%         plot(data{i},i*ones(1,length(data{i})),'Color',[0.3 0.3 0.3],'Marker','none','LineStyle','-');
    end
%     s=scatter(,5,color);
%     s.Marker='.';
end

try
    xlim(px)
catch
    xlim([0 1]);
end
% try
%     ylim(py)
% catch
%     ylim([0 1]);
% end

% fprintf('%d\n%d\n',px,py)

% disp('No data found for the following samples:')
% disp(EmptyTrial)

end
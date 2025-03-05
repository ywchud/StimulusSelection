function h=plotWE(x,y,SEM,C1,legstr)
    
%Data dimension: x:1D timepoints; y:1D data points; SEM:1D error data
%points; C1: rgd 1x3 vector; legstr: string for legend

%output:
% handle of plotted fig

    if ~exist('x','var') || isempty(x)
        x=1:length(y);
    end
    if ~exist('C1','var') || isempty(C1)
        C1=[0 0 1];
    end
    C2=C1*0.8;
    ALPHA=0.3;
    h=plotWFilledError(x,y,C1,y-SEM,y+SEM,C2,ALPHA);
    if exist('legstr','var') && ~isempty(legstr)
        h.DisplayName=legstr;
    end
end

function h=plotWFilledError(x,y,C1,y1,y2,C2,ALPHA)
    
x=toColumn(x)';
y=toColumn(y)';
y1=toColumn(y1)';
y2=toColumn(y2)';

    X=[x,fliplr(x)]; Y=[y1,fliplr(y2)]; 
    
    X(isnan(X))=0;
    Y(isnan(Y))=0;
    
    if ~isempty(Y)
        h=fill(X,Y,C2,'EdgeColor','none'); alpha(ALPHA)
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        hold on
        if ~isempty(y)
            try
                h=plot(x,y,'Color',C1,'Linestyle','-','Marker','none');
            catch
                h=plot(x,y,'Color','k','Linestyle','-','Marker','none');
            end
        end
    else
        h=fill(X,zeros(length(X),1),C2,'EdgeColor','none'); alpha(ALPHA)
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h=plot(x,zeros(length(x),1),'Color',C1);
        disp('Y is empty, plotting y=0')
    end
end
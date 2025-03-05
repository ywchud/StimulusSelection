function [Groups,Edges,corrTime] =DiscretizeBehvr(WhkBhvr,corrTime,binN)
    if iscell(WhkBhvr)
        try
            WhkBhvr=toColumn(cell2mat(WhkBhvr));
        catch
            WhkBhvr=toColumn(cell2mat(WhkBhvr'));
        end
    end
    if ~isempty(corrTime)    %if corrTime is not empty, validate equal dimension with WhkBhvr for finding Neuron.tuningcurve
        if iscell(corrTime)
            try
                corrTime=toColumn(cell2mat(corrTime));
            catch
                corrTime=toColumn(cell2mat(corrTime'));
            end
        end
        if length(WhkBhvr)~=length(corrTime)
            error('Unmatch input dimension')
        end
    end
    [Groups,Edges] = discretize(WhkBhvr,binN);
%     histogram(WhkBhvr,Edges)
end
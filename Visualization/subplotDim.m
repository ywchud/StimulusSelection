function [i,j]=subplotDim(ClustNum)
if isempty(ClustNum) || ClustNum==0
    error('ClustNum cannot be empty or 0')
end

    fprintf('%d subplots\n',ClustNum)

    i=ceil(ClustNum^(0.5));
    j=i;
    while i*j>=ClustNum
        j=j-1;
    end
    j=j+1;

end

% length(ListSpkClust)
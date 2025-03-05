function x=cell2vector(varargin)
x=varargin{1};
    if iscell(x)
        
        cellcontent=x{1};
        if isvector(cellcontent)
            try
                x=toColumn(cell2mat(x));
            catch
                x=toColumn(cell2mat(x'));
            end
        else
            error('Cell content is not a vector')
        end
    elseif isvector(x)
        disp('Input is already a vector')
    else
        error('Input is ot a cell type')
    end

end
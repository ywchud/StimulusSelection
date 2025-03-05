function Y=toColumn(X)
    if isvector(X)
         Y=reshape(X,[length(X),1]);
    elseif size(X,2)==length(X)   %set longer vectors as columns
        Y=reshape(X,[length(X),size(X,1)]);
    else
        Y=X;
    end   
end
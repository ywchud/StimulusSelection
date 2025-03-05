function X_padded=NANpad(X,padL,offset)
   %(X,padL,offset)
   %shift vector X by offset and pad with nan into a column vector X_padded(final length=padL)
   if ~isvector(X)
       error('X should be a vector for now, change code to include matrix if desired')
   end
   
   if isempty(offset)
       offset=0;
   end
   
%    if padL<length(X)+offset
%        error('padL should be larger than or equal to array length +offsest')
%    end
   
X=toColumn(X);

if padL+offset>=length(X)
    backnan=nan(padL-length(X)+offset,1);
    X=[X;backnan];
end
if offset<0
    X=[nan(-offset,1);X];
    X_padded=X(1:padL);
else
    X_padded=X(offset+1:offset+padL);
end

end
function sigf=removeFlicker(sig,n)
%takes a binary 1D signal and removes quick flicker of size n (default n=1)
ExemptFirst=1;

if isempty(sig)
    sigf=[];
    return;
end

if sum(~(sig==1 | sig==0))>0
    error('Input must be binary')
end

if isempty(n)
    n=1;
end

if ExemptFirst
    sigfirst=sig(1);
end

same=0;
for i=(1+n):length(sig)
    if sig(i-1)==sig(i)
        same=same+1;
    else
        if same<n
            for j=1:n
                sig(i-j)=sig(i);
                same=same+1;
            end
        else
            same=0;
        end
    end
end

if ExemptFirst
    sig(1)=sigfirst;
end
    sigf=sig;
end
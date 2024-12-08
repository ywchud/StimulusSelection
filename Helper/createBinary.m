function array=createBinary(totL,points,onesL)
% creates an array of zeros with length: totL,
% with segments of ones at all given points of length onesL
% if onesL is a vector of size==size(points), points(i):onesL(i)=1
% elseif onesL is singular, points(i):points(i)+onesL-1=1
if max(points)>totL
    error('Given points should be within the array length')
end

array=zeros(totL,1);


for i=1:length(points)
    
    if length(onesL)==length(points)
        endpt=min([onesL(i) totL]);  
    elseif length(onesL)==1
        endpt=min([points(i)+onesL-1 totL]);  
    else
        endpt=min([points(i) totL]);  
    end
        
    
     
    array(points(i):endpt)=1;
end

array=logical(array);
end
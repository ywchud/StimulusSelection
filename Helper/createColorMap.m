function cMap=createColorMap(NPoint,colorPivotMat)
        if isempty(colorPivotMat)
            colorPivotMat=[1 0 0;
                            0 0 0];
        end
        [L1,L2]=size(colorPivotMat);
        if L2~=3
            error('3 columns needed for rgb')
        end
        if L1<2
            error('More than 2 rows needed for color pivots')
        end
        
        c=colorPivotMat;
        
        gradient=[];
        for i=1:L1-1
            gradient=[gradient; linspace(c(i,1),c(i+1,1),10)' linspace(c(i,2),c(i+1,2),10)' linspace(c(i,3),c(i+1,3),10)'];
        end

%         gradient=resample(gradient,NPoint,size(gradient,1));
if NPoint<=size(colorPivotMat,1)
    cMap=colorPivotMat(1:NPoint,:);
else
    for i=1:3
        vecSize=size(gradient,1);
        cMap(:,i)=interp1(1:vecSize,gradient(:,i),linspace(1,vecSize,NPoint));
    end
%         gradient=round(gradient,2)*255;
end
end
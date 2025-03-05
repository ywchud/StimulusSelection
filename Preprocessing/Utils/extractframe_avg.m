function extractframe_avg(Exp,StimC,I)

i=StimC.WID;
arr=[0 0 0 0];arr(i)=1;
% Piston_i_tr=find(~isnan(Exp.Stim.Piston.Mat(:,i*2)));
Piston_i_tr=find(PistonComb(arr,Exp.Stim.Piston.Mat));
selectedtrials=Piston_i_tr;
if ~isprop(StimC,'TouchFrames') && ~isfield(StimC,'TouchFrames') %not property from obj nor field from struct
    StimC.TouchFrames=TouchFrameFromTime(Exp,StimC);
end
if exist('I','var')
    StimC=StimC.updateTouches(I);
end
    %get dim
    tr=selectedtrials(1);
    temp=extractframe(Exp,Exp.Stim.Piston.Cam(i(1)),tr,round(Exp.FrameN(tr)/2),0);
    temp=rgb2gray(temp);
    [H,W]= size(temp);
    F=nan(H,W,length(selectedtrials));
    %get averaged piston i frame
    for m=1:length(selectedtrials)
        tr=selectedtrials(m);
        try
        FrameRange=StimC.TouchFrames{tr}(1);
        
        f=extractframe(Exp,Exp.Stim.Piston.Cam(i(1)),tr,FrameRange,0);
        F(:,:,m)=rgb2gray(f);
        catch
        end

    end
    FF=uint8(mean(F,3,'omitnan'));
%     figure
    imshow(FF)
    try
    title(sprintf('%s n=%d',StimC.type,length(StimC.TouchTimes)))
    catch
    end
end
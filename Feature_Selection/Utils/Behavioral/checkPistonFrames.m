function checkPistonFrames(Exp)
%visually validate the frames when piston is considered out in
%Exp.PistonBuffer

temp=Exp.PistonBuffer';
weirdTrialsbool=(all(isnan(temp)));  %not belonging to any category, e.g. no piston trials
tp_PistonOutF=nan(Exp.TrN,1);
tp_PistonOutF(~weirdTrialsbool)=temp(~isnan(temp));

% tp_PistonOutT=cellfun(@(x,y) x(y),Exp.FrameT,num2cell(tp_PistonOutF));
clear tp_PistonOutT
for i=1:length(tp_PistonOutF)
    if ~isnan(Exp.FrameT{i}) & ~isnan(tp_PistonOutF(i))
        tp_PistonOutT(i,1)=Exp.FrameT{i}(tp_PistonOutF(i));
    else
        tp_PistonOutT(i,1)=nan;
    end
end
tp_PistonOutrelT=tp_PistonOutT-Exp.TrialStartT-Exp.Cam.Delay;

%
frames=nan([Exp.DLC.size(1:2),7,20]);
range=[-5:1:5];
for i=1:20 %first th trial
if ~isnan(tp_PistonOutF(i))
Frame=tp_PistonOutF(i)+range;
for k=1:length(Frame)
frames(:,:,k,i) = rgb2gray(extractframe(Exp, 1, i, Frame(k), 0));
end
end
i
end
figure
for k=1:10
    subplot(5,2,k)
    test=mean(frames(:,:,k,:),4,'omitnan');
    imshow(test/255)
    title(sprintf('t(s): %0.3f',range(k)/Exp.videoFps))
end
end
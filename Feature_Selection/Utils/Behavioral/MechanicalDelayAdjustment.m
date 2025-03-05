function Exp=MechanicalDelayAdjustment(Exp)

PistonN=Exp.Stim.Piston.num;
sampleN=3;
delayF=round(Exp.Stim.Piston.Delay*Exp.videoFps);

for i=1:PistonN
    fprintf('Finding delay frame for piston %d...\n',i)
    array=zeros(1,PistonN);array(i)=1;
    Piston_i_tr=find(PistonComb(array,Exp.Stim.Piston.Mat));
    randid=randperm(length(Piston_i_tr),sampleN);
    selectedtrials=Piston_i_tr(randid);
    PiDelayF=zeros(sampleN,1);
    z=0;
    for j=1:sampleN
        tr=selectedtrials(j);
        while z<=0  %only play for first vid
            playVid(Exp,Exp.Stim.Piston.Cam(i),tr,[],[],[],[delayF 150],2) 
            commandwindow
            z = input(sprintf('FrameN when piston is out: (<=0 to replay)'));
            close
        end
        ok=0;
        figure, hold off
        while ~ok
            try
            f=extractframe(Exp,Exp.Stim.Piston.Cam(i),tr,delayF+z);  
            catch
              error('Frame not found')  
            end
            imagesc(f),axis('image'),title(sprintf('Adjust z=%d',z)) 
            [x,~] = ginput(1);
            if x<0
                z=z-1;
            elseif x>size(f,2)
                z=z+1;
            else
                ok=1;
            end
        end
        close
        PiDelayF(j)=z;
    end
    Exp.PistonBuffer(i)=mean(PiDelayF);
end


end
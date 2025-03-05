function Exp=findPistonParam(Exp,duration)
    Exp.Stim.Piston.Cam=[];
    Exp.Stim.Piston.PWOI=[];
    Exp.Stim.Piston.TLOI=[];

for i=1:Exp.Stim.Piston.tot
    
    if ~Exp.Stim.Piston.active(i)
        fprintf('Piston %d, aka Whisker %d is inactive, skipping...\n',Exp.Stim.Piston.ID(i),i)
        continue
    end
    
    fprintf('Getting parameters for Piston %d, aka Whisker %d\n',Exp.Stim.Piston.ID(i),i)
    
    PistonArray=zeros(1,Exp.Stim.Piston.tot);
    PistonArray(i)=1;
    x=0;y=0;z=0;
    while x==0 || y==0 || z==0
        k=find(PistonComb(PistonArray,Exp.Stim.Piston.Mat));
        if isempty(k)  %if there is no single whisker trial
            disp('No single whisker trial found, using multi touch trials (may not work)...')
            k=find(~isnan(Exp.Stim.Piston.Mat(:,i*2)));
            if isempty(k)
                error('Trial with said Piston not found. Check if piston exists.')
            end
        end
        k=k(randperm(length(k),1));   %get random trial with piston array
        for j=1:Exp.Cam.num
            playVid(Exp,j,k,[],[],0,duration,0)  %playVid(Exp,CamN,VidN,scat,Whisker,playAngle,duration,slow)
            title(sprintf('Camera %d Trial %d Piston %d',j,k,Exp.Stim.Piston.ID(i)))
        end
        x = input(sprintf('Which Camera is Piston %d in? Enter 0 to look at another trial',Exp.Stim.Piston.ID(i)));
        if x~=0
            y = input(sprintf('Which Whisker is Piston %d touching? Enter 0 to look at another trial',Exp.Stim.Piston.ID(i)));
        end
        if x~=0 && y~=0
            z = input(sprintf('Which DLC Label will you check for touch ROI (Red(1)->Blue(n))? Enter 0 to look at another trial'));
        end
        for j=1:Exp.Cam.num
            close
        end
    end
    
    Exp.Stim.Piston.Cam(i)=x;
    Exp.Stim.Piston.PWOI(i)=y;
    Exp.Stim.Piston.TLOI(i)=z;
    fprintf('Piston %d .Cam .PWOI .TLOI saved\n',Exp.Stim.Piston.ID(i))
    %find ROI for TLOI
    
        
    
end


disp('Piston .Cam .PWOI .TLOI  obtained')


end
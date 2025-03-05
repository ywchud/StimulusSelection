function Whisker=findWhiskerParam(Exp,Whisker)
    if Exp.Stim.Piston.num~=length(Whisker)
        error('Expecting %d whiskers but found %d',Exp.Stim.Piston.num,length(Whisker))
    end   
   
    
    disp('Obtaining whisker angle')
    for i=1:Exp.Stim.Piston.num
        x0=Whisker(i).whiskerbaseX;
        y0=Whisker(i).whiskerbaseY;
        LabelN=Whisker(i).LabelN;
        for j=1:Exp.TrN
            if ismember(j,Exp.BadTrials)
                fprintf('Skipping bad Trial %d\n',j)
            else
                Whisker(i).trial(j)=Whisker(i).trial(j).findAngle(x0,y0,LabelN) ;
                fprintf('Obtained Angle for Whisker %d Trial %d\n',i,j)
            end
        end
    end
    disp('Whisker angle obtained') 
    
    disp('Obtaining whisker phase')
    for i=1:Exp.Stim.Piston.num
        Fc1=8;Fc2=40;order=4;
        for j=1:Exp.TrN
            if ismember(j,Exp.BadTrials)
                fprintf('Skipping bad Trial %d\n',j)
            else
                Whisker(i).trial(j)=Whisker(i).trial(j).findPhase(Fc1, Fc2, order) ;
                fprintf('Obtained Phase for Whisker %d Trial %d\n',i,j)
            end
        end
    end
    disp('Whisker phase obtained') 
    
    disp('Obtaining whisker curvature')
    for i=1:Exp.Stim.Piston.num
        WBx=Whisker(i).whiskerbaseX;
        WBy=Whisker(i).whiskerbaseY;
        LabelN=Whisker(i).LabelN;
        for j=1:Exp.TrN
            if ismember(j,Exp.BadTrials)
                fprintf('Skipping bad Trial %d\n',j)
            else
                Whisker(i).trial(j)=Whisker(i).trial(j).findCurvature(LabelN,[],WBx,WBy) ;%(obj,LabelN,Points,WBx,WBy)
                fprintf('Obtained Curv for Whisker %d Trial %d\n',i,j)
            end
        end
    end
    disp('Whisker curvature obtained') 
    
    
    
    
    disp('Angle/Phase/Curv obtained')
end
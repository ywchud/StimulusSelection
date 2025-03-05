classdef Trial
    properties( Access = public )
         % All your arguments here with their default values
         X;%signal
         Y;%signal
         FrameN;
         Angle;  %signal
         CenAngle; %signal
         
         %The followings are not set as a signal object to save memory, if
         %you need to use signal object functions, set it as object.sig in another script
         
         Phase;  %not signal %from CenAngle.Phase  
         Curvature; %not signal
         CurvDebase; %not signal
         
         
         Pks; 
         Trs;
         Setpoint;
         LocsP;
         LocsT;%correspond to most protracted(touch)
         Amplitude;  %Amplitude can be found either in .Amplitude (using butterworth) or .Envelope (by substracting upper envelope with lower envelope) in signal functions 
         
         P_Frame;
         R_Frame;         
         TouchFrame;
         ReleaseFrame;
%          ProtractionFrame;         %same as LocsP(GoodTouch)
%          RetractionFrame;
         GoodTouch;
         GapDistance;
         GapFrame;
         TouchGap;
         GoodType;
    end
    methods
        function obj = Trial( XX,YY,labelN,fs)
          % (constructor)
          obj.X=Signal.empty;
          obj.Y=Signal.empty;
          obj.Angle=Signal.empty;
          obj.CenAngle=Signal.empty;
          
          for i=1:labelN
            obj.X(i)    = Signal( [],[], XX(:,i), fs ) ;
            obj.Y(i)    = Signal( [],[], YY(:,i), fs ) ;
            obj.Angle(i)= Signal( [],[], [], fs ) ;
          end
          obj.CenAngle=Signal( [],[], [], fs ) ;
          obj.FrameN=length(obj.X(1).sig);
        end
    
        function obj=findAngle(obj,x0,y0,LabelN) 
            if isempty(x0) || isempty(y0)
                x0=0;y0=0;
                disp('x0 and y0 not found, assuming reference point at (0,0)')
            end
            AngleMat=nan(obj.FrameN,LabelN);
            for k=1:LabelN
%                 x=obj.X(k).smooth;  
%                 y=obj.Y(k).smooth;
                x=obj.X(k).sig;  
                y=obj.Y(k).sig;     
                AngleSig=atan((y-y0)./(x-x0));
                if sum(x<x0)>obj.FrameN/2  %if more than half of the labels are to the left of whiskerbase
                    AngleSig=-AngleSig;     %flip signs so all peaks are protraction
                end
                obj.Angle(k)=Signal([],[],AngleSig,[]);           
                AngleMat(:,k)=obj.Angle(k).sig;
            end         
            obj.CenAngle.sig=mean(AngleMat,2);
        end
        
        function obj=findPhase(obj,Fc1, Fc2, order) %for band pass butter filter
            if isempty(Fc1) || isempty(Fc2) 
                Fc1=8;
                Fc2=40;
                fprintf('Fc1 and Fc2 not found, assuming bandpass from %d to %d\n',Fc1,Fc2)
            end
            if isempty(order)
                order=4;
                fprintf('Butter order not found, assuming order as %d\n',order)
            end
            obj.Phase=obj.CenAngle.Phase(Fc1, Fc2, order);
        end
        
        function DeltaPhase=findPhaseDiff(obj,obj2,isabs)
            if isempty(obj.Phase) || isempty(obj2.Phase)
                error('Phase is empty, use findPhase function')
            end
            
            if isabs
%                 DeltaPhase=exp(1i*obj.Phase)-exp(1i*obj2.Phase)
                DeltaPhase=abs(obj.Phase-obj2.Phase);
            else 
                DeltaPhase=obj.Phase-obj2.Phase;
            end
        end
        
        
        function obj=findCurvature(obj,LabelN,Points,WBx,WBy) 
            Curv=nan(obj.FrameN,1);   
            x=nan(obj.FrameN,LabelN+1); %last column is whiskerbase
            y=nan(obj.FrameN,LabelN+1);
            for k=1:LabelN
                x(:,k)=obj.X(k).sig;
                y(:,k)=obj.Y(k).sig;
            end
            x(:,LabelN+1)=WBx*ones(obj.FrameN,1);
            y(:,LabelN+1)=WBy*ones(obj.FrameN,1);                 
            if isempty(Points)      
                Points =1:LabelN;
            end
            for n = 1:obj.FrameN  % frames no.
                P = combnk(Points,3);
                curv=NaN(size(P,1),1);
                for k=1:size(P,1)  %combination no.
                    x1=x(n,P(k,1));
                    y1=y(n,P(k,1));
                    x2=x(n,P(k,2));
                    y2=y(n,P(k,2));
                    x3=x(n,P(k,3));
                    y3=y(n,P(k,3));

                    A=1/2*abs(det([x1 y1 1;x2 y2 1;x3 y3 1]));
                    curv(k)=4*A/norm([x1-x2 y1-y2])/norm([x1-x3 y1-y3])/norm([x3-x2 y3-y2]); %magnitude always positive
                    %check direction (i.e. mid point above(-) or below(+) line
                    %of edge points
                    m=(y3-y1)/(x3-x1);
                    y_prime=m*(x2-x1)+y1;
                    if y2<y_prime  %change sign if above line
                        curv(k)=-curv(k);
                    end
                end
                Curv(n) = mean(curv);  %average of all possible combs based on given labels per whisker
            end
            Curv = smoothdata(Curv,'sgolay',12);
            obj.Curvature=Curv;
        end
        
        function obj=findEnvelope(obj,th,absMin,maxDur,TrsPksSamesize)
            if isempty(th) && isempty(absMin) 
                th=[];
                absMin=0.03;
                fprintf('th and absMin not found, assuming empty th and absMin= %d\n',absMin)
            end

            [pks,trs,setPt,locsP,locsT,amplitude]=obj.CenAngle.Envelope(th,absMin,maxDur,TrsPksSamesize);  %th,absMin=0.03(empirical for rads with whiskerbase),maxDur
            
            obj.Pks=pks; %correspond to most protracted(touch)
            obj.Trs=trs;
            obj.Setpoint=setPt;
            obj.LocsP=locsP;
            obj.LocsT=locsT;
            obj.Amplitude=amplitude;
        end
        
        function sth(obj)
        end
        
    end
end
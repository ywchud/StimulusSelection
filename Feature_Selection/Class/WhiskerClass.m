classdef WhiskerClass
     properties( Access = public )
         % All your arguments here with their default values
         W_id  %whisker id
         ID  %piston id
         Cam
         PWOI   %as in the nth whisker labelled in DLC just for this camera, different from ID
         TLOI
         LabelN

         trial  %class Trial, contains angles and stuff
         whiskerbaseX
         whiskerbaseY

         
         
         
     end
     methods
         
          function obj = WhiskerClass(Exp,W_id,Xmat,Ymat)

                        % (constructor)
            if W_id==0
                return;
            end
          obj.W_id=W_id;
          if ~isempty(Exp.Stim.Piston.ID)
                obj.ID    = Exp.Stim.Piston.ID(W_id) ;
          end
          if ~isempty(Exp.Stim.Piston.Cam)
                obj.Cam      = Exp.Stim.Piston.Cam(W_id) ; 
          end
          if ~isempty(Exp.Stim.Piston.PWOI)
                obj.PWOI = Exp.Stim.Piston.PWOI(W_id);
          end
          if ~isempty(Exp.Stim.Piston.TLOI)
                obj.TLOI = Exp.Stim.Piston.TLOI(W_id) ;
          end
          if ~isempty(Exp.DLC.LblPerW)
                obj.LabelN = Exp.DLC.LblPerW ;
          end
          if ~isempty(Xmat) && ~isempty(Ymat)  
              obj.trial=Trial.empty;
%               matID=[1:obj.LabelN]+obj.LabelN*(obj.PWOI-1);
              matID=[1:obj.LabelN]+obj.LabelN*(W_id-1);
              for i=1:Exp.TrN
                  if length(Exp.videoFps)==1
                    fs=Exp.videoFps;
                  else
                    fs=Exp.videoFps(i);
                  end
                obj.trial(i)=Trial(Xmat{obj.Cam,i}(:,matID),Ymat{obj.Cam,i}(:,matID),obj.LabelN,fs);
              end 
          end
          
          
          end
         
         function obj=Whiskerbase(obj,Exp,sampleN)   %find whiskerbase X & Y
            if isempty(sampleN)
                sampleN=3;
            end
            x=nan(sampleN,1);
            y=nan(sampleN,1);
            for k=1:sampleN
                
                VidN=rand;
                VidN=ceil(VidN*Exp.TrN);
                
                playVid(Exp,obj.Cam,VidN,[],[],0,[0.5 100],1) %playVid(Exp,CamN,VidN,scat,Whisker,playAngle,duration,slow)
                title(sprintf('Pick (%d/%d) Whiskerbase for Whisker %d in Camera %d',k,sampleN,obj.W_id,obj.Cam))
                [x(k),y(k)]=ginput(1);
                disp(x(k))
                disp(y(k))
                hold on
                scatter(x(:),y(:))
                title('Click anywhere to close')
                [~,~]=ginput(1);
                close
            end
            obj.whiskerbaseX=mean(x);
            obj.whiskerbaseY=mean(y);
         end

          
         function TouchTime(obj)
         end

         
     end
end
classdef Signal
     properties( Access = public )
         % All your arguments here with their default values
         Port
         type='unknown'
         sig
         fs
         
         sig_filtered
    
     end

     methods
         function obj = Signal( Port,type, sig, fs )
          % (constructor)
          if ~isempty(Port)
            obj.Port    = Port ;
          end
          if ~isempty(sig)
            obj.sig      = sig ; 
            obj=obj.dimcheck;
          end
          
          
          if ~isempty(type)
            obj.type    = type ;  %binary,continuous
          else
              obj = obj.typecheck;
          end
          if ~isempty(fs)
            obj.fs = fs ;     
          end
         end
         
         function obj = typecheck(obj)      %assign type
             if sum(obj.sig ~= 0 & obj.sig ~= 1)==0         
                 obj.type='binary';
             else
                 obj.type='cont';
             end
 
         end
         
         function obj = dimcheck(obj)       %reshape 1D signal to column
             if isvector(obj.sig)
                 obj.sig=reshape(obj.sig,[length(obj.sig),1]);
             end    
         end
         
         function onset=onset(obj)
             if ~strcmp(obj.type,'binary') 
                obj = typecheck(obj);
             end          
             if strcmp(obj.type,'binary') 
                 onset = find(diff([0; obj.sig; 0])>0);
             else
                 error('signal must be binary')
             end
         end
         
         function offset=offset(obj)
             if ~strcmp(obj.type,'binary') 
                obj = typecheck(obj);
             end
             if strcmp(obj.type,'binary') 
                 offset = find(diff([0; obj.sig; 0])<0)-1;
             else
                 error('signal must be binary')
             end
         end
         
         function [Sig_smooth,obj]=MASmooth(obj,window,replace)
             if isempty(window)
                 window=12;
             end
             Sig_smooth=smoothdata(obj.sig,'sgolay',window);
             if replace
                obj.sig=Sig_smooth;
             end
         end
         
         function [Sig_smooth,obj]=ButterSmooth(obj,Fc1, Fc2, order,replace)
             if isempty(Fc1) || isempty(Fc2) || isempty(order)
                 error('Fc1,Fc2,order cannot be empty')
             end
             Sig_smooth=genButterFilter(obj.sig,Fc1, Fc2, order,'butter_acausal',obj.fs);
             if replace
                obj.sig=Sig_smooth;
             end
         end
         
         function [pks,trs,setPt,locsP,locsT,amplitude]=Envelope(obj,th,absMin,maxDur,TrsPksSamesize)
             if isempty(absMin)
                absMin=0.03;
                maxDur=[];
                th=[];                
             end
             if isempty(TrsPksSamesize)
                 TrsPksSamesize=0;
             end
             [pks,trs,setPt,locsP,locsT] = findEnvelope(obj.sig,[],th,absMin,maxDur,TrsPksSamesize);
             amplitude=pks-trs;
         end
         

         function phase=Phase(obj,Fc1, Fc2, order)       
             [Sig_smooth,~]=ButterSmooth(obj,Fc1, Fc2, order,0);
             phase=angle(hilbert(Sig_smooth));
         end
         
         function amplitude=Amplitude(obj,Fc1, Fc2, order)       
             [Sig_smooth,~]=ButterSmooth(obj,Fc1, Fc2, order,0);
             amplitude=abs(hilbert(Sig_smooth));
         end
         
         
         function AutoCorr(obj)
         end
         
     end
end
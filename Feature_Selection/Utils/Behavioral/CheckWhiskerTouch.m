% ROIradius=10;
function CheckWhiskerTouch(Exp,Whisker,Trial,duration,slow)

if length(Whisker)~=1
    error('Currently only support one whisker at a time')
end

if isempty(Trial)
    Trial=1:Exp.TrN;
end
if isempty(slow)
    slow=0;
end
if isempty(duration)
    duration=[];
end
for i=Trial
        playVid(Exp,Whisker.Cam,i,[],Whisker,"CheckTouchAngle",duration,slow)                   
end
end

%                
%                         
%                         if size(Tcoor,2)==2
%                             viscircles([Tcoor(WOIS,1) Tcoor(WOIS,2)], ROIradius,'Color','b');
%                             viscircles([Tcoor(PP,1) Tcoor(PP,2)], ROIradius,'Color','b');
%                             viscircles([Tcoor(PPP,1) Tcoor(PPP,2)], ROIradius,'Color','b');
%                         elseif size(Tcoor,2)==4
%                             equidistalPlot(Tcoor(WOIS,1),Tcoor(WOIS,2),Tcoor(WOIS,3),Tcoor(WOIS,4),ROIradius,'b')
%                             equidistalPlot(Tcoor(PP,1),Tcoor(PP,2),Tcoor(PP,3),Tcoor(PP,4),ROIradius,'b')
%                             equidistalPlot(Tcoor(PPP,1),Tcoor(PPP,2),Tcoor(PPP,3),Tcoor(PPP,4),ROIradius,'b')
%                         end
% 

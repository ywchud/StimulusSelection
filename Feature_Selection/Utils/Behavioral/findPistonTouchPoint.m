 %determine touch region (x,y)


% n=n+sum(n>=missing)
% n=91;
% %rmb to change path in write_labed_video
datapath = 'E:\Darren\ScottBC_Halo_200724_150610\ScottBC_Halo_200724_150610\ScottBC_Halo_200724_vids';
d_name='ScottBC_Halo_200724_%dDLC_resnet50_BC_Halo_200724Jul31shuffle1_300000.csv';
v_name='ScottBC_Halo_200724_%d';
fraction=5; %number of times a ROI is sampled acros trials for each piston (used for interpolation to address drifting piston)

clear touchCoor
clear Whiskerbase
%%
for WOIS=1:nofWhiskers
    array=zeros(1,nofWhiskers);array(WOIS)=1;
%     n=find(PistonComb(array,Piston),1,'first');
    N=find(PistonComb(array,Piston));
    touchCoor(WOIS)=TouchCoor;
    touchCoor(WOIS).trialNum=length(N);
for frac=1:fraction+1
id=ceil(length(N)/fraction*(frac-1));
if id==0, n=N(1); 
else
    n=N(id);
end
dataname = sprintf(d_name,n);
videoname = sprintf(v_name,n);  

figure(2)

try
    V = VideoReader(fullfile(datapath,[videoname '_labeled.mp4']));
catch
    fprintf('Writing labelled video: %d\n',n)
    title(sprintf('Writing labelled video: %d',n))
    write_labed_video(datapath,dataname,[videoname '.mp4']);%generate labeled video first
    V = VideoReader(fullfile(datapath,[videoname '_labeled.mp4']));
end
%%
% find(PistonComb([1 0 0 0],Piston))
videotime=0*60+10;  %time(s) in mp4 video
slow=1;
	
FrameI=videotime*mp4Fps;

hold off;

% coor=[nan nan nan nan];
ok=0;
while(~ok)
    for i=60*2:60*4
    frame = read(V,FrameI+i);
    imagesc(frame),axis('image')
        if exist('Whiskerbase'),hold on,scatter(Whiskerbase(1),Whiskerbase(2),'o'),hold off,  
        else
            title('Choose Whiskerbase')
            [Whiskerbase(1),Whiskerbase(2)]=ginput(1);
        end
%         if exist('touch_1'),hold on,scatter([touch_1(1) touch_1(3)],[touch_1(2) touch_1(4)],'o'),hold off,  end
%         if exist('touch_2'),hold on,scatter([touch_2(1) touch_2(3)],[touch_2(2) touch_2(4)],'o'),hold off,  end
%         if exist('touch_3'),hold on,scatter([touch_3(1) touch_3(3)],[touch_3(2) touch_3(4)],'o'),hold off,  end
%         if exist('touch_4'),hold on,scatter([touch_4(1) touch_4(3)],[touch_4(2) touch_4(4)],'o'),hold off,  end
        for w=1:WOIS
            if ~isempty(touchCoor(w).x1),hold on,scatter([touchCoor(w).x1 touchCoor(w).x2],[touchCoor(w).y1 touchCoor(w).y2],'o'),hold off,  end 
        end
%         if ~isempty(touchCoor(2).x1),hold on,scatter([touchCoor(2).x1 touchCoor(2).x2],[touchCoor(2).y1 touchCoor(2).y2],'o'),hold off,  end
%         if ~isempty(touchCoor(3).x1),hold on,scatter([touchCoor(3).x1 touchCoor(3).x2],[touchCoor(3).y1 touchCoor(3).y2],'o'),hold off,  end
%         if ~isempty(touchCoor(4).x1),hold on,scatter([touchCoor(4).x1 touchCoor(4).x2],[touchCoor(4).y1 touchCoor(4).y2],'o'),hold off,  end


                title(sprintf('%d',i+600))
%             hold on,scatter(Tcoor(4,1),Tcoor(4,2),'o'),hold off
%              viscircles([Tcoor(4,1) Tcoor(4,2)], 15,'Color','b');

    pause(1/mp4Fps*2*slow)
    end
    
    title(sprintf('Piston: %d, Frac: %d',WOIS,frac))
    [touchCoor(WOIS).x1(frac),touchCoor(WOIS).y1(frac)]=ginput(1);
    [touchCoor(WOIS).x2(frac),touchCoor(WOIS).y2(frac)]=ginput(1);
    title(sprintf('Click left upper corner if OK'))
    [check1,check2]=ginput(1);
    ok=norm([check1 check2])<150;
%     [coor(2,1),coor(2,2)]=ginput(1);    %30x30 pixels seem ok
end
% touch_1=coor;       %touch coor
% touch_2=coor;
% touch_3=coor;
%  touch_4=coor;
% Whiskerbase=coor(1:2);


end
title(sprintf('Interpolating for Piston %d',WOIS))
n=[1 ceil(length(N)/fraction*(1:fraction))];
touchCoor(WOIS).x1 = interp1(n,touchCoor(WOIS).x1,1:touchCoor(WOIS).trialNum);
touchCoor(WOIS).x2 = interp1(n,touchCoor(WOIS).x2,1:touchCoor(WOIS).trialNum);
touchCoor(WOIS).y1 = interp1(n,touchCoor(WOIS).y1,1:touchCoor(WOIS).trialNum);
touchCoor(WOIS).y2 = interp1(n,touchCoor(WOIS).y2,1:touchCoor(WOIS).trialNum);
end
%plz rmb to save Tcoor when you are done
% Tcoor=[touch_1;touch_2;touch_3;touch_4];



filename = [Path '\COOR.mat'];
% save(filename,'Tcoor')
save(filename,'touchCoor','Whiskerbase')
title(sprintf('saved'))


% pist_1=coor;        %piston coor (so if the labeled dot is below the piston coor, we know it is a bad touch)
% pist_2=coor;
% pist_3=coor;
% pist_4=coor;
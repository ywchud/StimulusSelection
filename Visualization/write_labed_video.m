function write_labed_video(datapath,dataname,videoname,p)
if ~exist('p','var')
     % parameter does not exist, so default it to something
      offsetx = 0;offsety = 0;
elseif p.iscrop
    offsetx=p.crop_x1;
    offsety=p.crop_y1;
else
    offsetx = 0;offsety = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Darren 5/10/2019

% arg n is the trial number to be labelled
% Now remove points out of frame
% auto generation of color gradient and number of labeled points from csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load file
% n = 25
% datapath = 'Z:\data\Hyein\ExpVids\HyeinScottBC_EM107_200528_CAM1';
% dataname = sprintf('HyeinScottBC_EM107_200528_1_%dDeepCut_resnet50_Bilateral200528Jun1shuffle1_560000.csv',n);
% videoname = sprintf('HyeinScottBC_EM107_200528_1_%d.mp4',n);  
% 
% datapath = 'E:\ExpVids\DarrenBCPistonJAWS2_200123'; 
% dataname = sprintf('DarrenBCPistonJAWS2_%dDLC_resnet50_DarrenBCPistonJAWS2_200123Jan28shuffle1_400000.csv',n);
% videoname = sprintf('DarrenBCPistonJAWS2_%d.mp4',n);      %['ScottNubz_200312\' ]
newfilename = strcat(videoname(1:end-4),'_labeled');
Labelpath = fullfile(datapath,'Labelled');

W = csvread(fullfile(datapath,dataname),3,1);
v = VideoReader(fullfile(datapath,videoname));
try
    w = VideoWriter(fullfile(Labelpath,newfilename),'MPEG-4');
catch
    temp=cd;
    cd(datapath)
    mkdir Labelled
    disp('Created Labelled folder')
    cd(temp)
    w = VideoWriter(fullfile(Labelpath,newfilename),'MPEG-4');
end

%% rearrange the coordinate matrix
W = reshape(W,[size(W,1),3,size(W,2)/3]);

W(:,2,:)=W(:,2,:)+offsety;
W(:,1,:)=W(:,1,:)+offsetx;
%%
%figure (1)
pixelHW=1;
disp('Writing labeled video...')
    open(w);
    gradient=colorcode(size(W,3));
for i = 1 : size(W,1)
    video = read(v,i);
    try
    for j=1:size(W,3)
        video(NegToZero(W(i,2,j)-pixelHW:W(i,2,j)+pixelHW),NegToZero(W(i,1,j)-pixelHW:W(i,1,j)+pixelHW),1) = gradient(j,1); 
        video(NegToZero(W(i,2,j)-pixelHW:W(i,2,j)+pixelHW),NegToZero(W(i,1,j)-pixelHW:W(i,1,j)+pixelHW),2) = gradient(j,2);
        video(NegToZero(W(i,2,j)-pixelHW:W(i,2,j)+pixelHW),NegToZero(W(i,1,j)-pixelHW:W(i,1,j)+pixelHW),3) = gradient(j,3);        
    end
    catch
        fprintf('Frame Error: %d\nLabel: %d\n',i,j)
    end
%    imshow(video);
% 
%     F = getframe(gcf);
    video=video(1:v.Height,1:v.Width,:);
    writeVideo(w,video);
%     pause(1/40);

end
    close(w);
    disp('Writing done')
%% Nested functions
    function Y=NegToZero(X)
        Y=ceil((abs(X)+X)/2);
        Y(Y<1)=1;
    end
    function Y=InsideFrame(X,limit)
        if X>=limit
            Y=limit-1;
        else
            Y=X;
        end
    end
    function gradient=colorcode(NPoint)
        c1=[1 0 0];
        c2=[1 1 0];
        c3=[0 1 0];
        c4=[0 1 1];
        c5=[0 0 1];
        c6=[1 0 1];
        gradient=[];
        gradient=[gradient; linspace(c1(1),c2(1),10)' linspace(c1(2),c2(2),10)' linspace(c1(3),c2(3),10)'];
        gradient=[gradient; linspace(c2(1),c3(1),10)' linspace(c2(2),c3(2),10)' linspace(c2(3),c3(3),10)'];
        gradient=[gradient; linspace(c3(1),c4(1),10)' linspace(c3(2),c4(2),10)' linspace(c3(3),c4(3),10)'];
        gradient=[gradient; linspace(c4(1),c5(1),10)' linspace(c4(2),c5(2),10)' linspace(c4(3),c5(3),10)'];
        gradient=[gradient; linspace(c5(1),c6(1),10)' linspace(c5(2),c6(2),10)' linspace(c5(3),c6(3),10)'];
        gradient=resample(gradient,NPoint,size(gradient,1));
        gradient=round(gradient,2)*255;

    end
end

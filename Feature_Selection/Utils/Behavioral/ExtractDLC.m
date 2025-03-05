function [Exp,Xmat,Ymat]=ExtractDLC(Exp,Fc1, Fc2, order)

WhiskerLabelsX=[];
WhiskerLabelsY=[];
Xmat=[];
Ymat=[];
%set offset 
if Exp.DLC.iscrop
    offsetx=Exp.DLC.crop_x1;
    offsety=Exp.DLC.crop_y1;
else
    offsetx = 0;offsety = 0;
end

for CamN=1:Exp.Cam.num
datapath=Exp.Path.vid{CamN};

WhiskerRawData=[];

fprintf('Extracting .csv data from camera %d...\n',CamN)
for VidN = 1:Exp.TrN
    dataname = sprintf(Exp.Path.csvName{CamN},VidN);
    WhiskerRawData{1,VidN} = csvread(fullfile(datapath,dataname),3,1);
end
% missingfile
%% Whiskerbase
% find Whiskerbase with findPistonTouchPoint.m if not labeled
% may be used for finding curvature
%%
%Creates separate cell arrays for reference label and whisker labels
%Removes unncessary information from raw data (e.g. label likelihood)



nofLabels = size(WhiskerRawData{1},2)/3;
ExpectedLbllN=Exp.DLC.LblPerW*(Exp.Stim.Piston.num/Exp.Cam.num);
Exp.DLC.issnout=false;
if nofLabels~=ExpectedLbllN
    if nofLabels-ExpectedLbllN==1
        disp('1 extra label point found, assuming last label as snout position')
        Exp.DLC.issnout=true;
        Exp.DLC.extraLabelsN=1;
    elseif nofLabels-ExpectedLbllN>1
        fprintf('Labels number mismatched, expecting %d but found %d\n',ExpectedLbllN,nofLabels)
        Exp.DLC.issnout=true;
        Exp.DLC.extraLabelsN=nofLabels-ExpectedLbllN;
    end
end
temp = [];tempX=[];tempY=[];
for i = 0:nofLabels - 1   
%     temp = [temp, 3*i+1, 3*i+2];
    tempX = [tempX, 3*i+1];
    tempY = [tempY, 3*i+2];
end
% WhiskerLabels=[];
disp('Restructuring into X and Y...')

for i = 1:size(WhiskerRawData,2) 
%     WhiskerLabels{CamN,i} = WhiskerRawData{i}(:, temp);
    WhiskerLabelsX{CamN,i} = WhiskerRawData{i}(:, tempX)+offsetx;
    WhiskerLabelsY{CamN,i} = WhiskerRawData{i}(:, tempY)+offsety;
    for j = 0:nofLabels - 1
        for k=2:size(WhiskerRawData{i},1)
            if  WhiskerRawData{i}(k, 3*j+3)<0.8
%                 WhiskerLabels{CamN,i}(k,2*j+1:2*j+2)=WhiskerLabels{CamN,i}(k-1,2*j+1:2*j+2);
                WhiskerLabelsX{CamN,i}(k,j+1)=WhiskerLabelsX{CamN,i}(k-1,j+1);
                WhiskerLabelsY{CamN,i}(k,j+1)=WhiskerLabelsY{CamN,i}(k-1,j+1);
            end
        end
    end
end
disp('Smoothing data...')
% butter  kinda stupid dont do this, it removes 0 DC
% for i = 1:size(WhiskerLabelsX,2) 
%     for j=1:size(WhiskerLabelsX{CamN,i},2)         
%         Xmat{CamN,i}(:,j)=genButterFilter(WhiskerLabelsX{CamN,i}(:,j),Fc1, Fc2, order,'butter_acausal',Exp.videoFps);
%         Ymat{CamN,i}(:,j)=genButterFilter(WhiskerLabelsY{CamN,i}(:,j),Fc1, Fc2, order,'butter_acausal',Exp.videoFps);
%     end
% end
%moving average
for i = 1:size(WhiskerLabelsX,2) 
    for j=1:size(WhiskerLabelsX{CamN,i},2) 
%         WhiskerLabelsSmooth{CamN,i}(:,j) = smoothdata(WhiskerLabels{CamN,i}(:,j),'sgolay',12);
        Xmat{CamN,i}(:,j) = smoothdata(WhiskerLabelsX{CamN,i}(:,j),'sgolay',12);
        Ymat{CamN,i}(:,j) = smoothdata(WhiskerLabelsY{CamN,i}(:,j),'sgolay',12);              
    end
end
fprintf('Camera %d done\n',CamN)

end


end
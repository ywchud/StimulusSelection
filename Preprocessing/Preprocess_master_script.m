% This script saves config parameters into Exp and
% extracts related digital and analog signals
%% Clean wipe
clear, clc, close all
%% code dir
codeDir='E:\StimulusSelection';
cd(codeDir);
addpath(genpath(fullfile(codeDir,'')))  
%% General info config
Exp.name = 'ExampleMouse8';
Exp.date ='250101';
Exp.ID=1;
Exp.Description = 'Stimulus Selection, LSC ephys';
Exp.Animal.cage='ExampleCage';
Exp.Animal.ID=[1];
%Camera
Exp.Cam.num=1;

%Path
Exp.Path.data = uigetdir(fullfile(codeDir,'Data'),'Select data folder');
for i=1:Exp.Cam.num
    Exp.Path.vid{i} = uigetdir(Exp.Path.data,sprintf('Select camera %d folder',i));
end
Exp.Path.save = fullfile(Exp.codeDir,'Saves',Exp.name);
mkdir(Exp.Path.save)
folder(Exp)
Exp.codeDir=codeDir;
% deets(Exp)
%% Check video files
Exp=get_video_files(Exp,codeDir);
%----------- Write/Play DLC labelled video
n=79;
write_labed_video(Exp.Path.vid{1},sprintf(Exp.Path.csvName{1},n),sprintf(Exp.Path.vidName{1},n))
playVid(Exp,[1],n,[],[],0,[0 150],2)
%%  LFP 
Exp=extract_electrode_dim(Exp);
Exp.KSort.LFP_fc=200;  %low pass freq
Exp.KSort.LFP_fs=1000; %sampling freq

% % amplifier.dat is not provided in this example due to its large size
% LFP=getLFPbyChannel(Exp.Path.data,128,20000,Exp.KSort.LFP_fc,Exp.KSort.LFP_fs);
% Exp.KSort.LFP_fs=length(LFP)/length(Exp.t)*Exp.Fs;

%% Digital Input Config
% % unilateral rig
Exp.Intan.PistonPort=[5 16 1 2];
% Exp.Intan.BarPort=6;
Exp.Intan.CamPort=8;
Exp.Intan.TrialPort=7;
Exp.Intan.WheelPort=15; 
Exp.Intan.SpoutPort=11;
Exp.Intan.Solenoid=13;
Exp.Intan.MissPort=9;
Exp.Intan.FAPort=10;

% % bilateral
% Exp.Intan.PistonPort=9:12;
% Exp.Intan.CamPort=8;
% Exp.Intan.TrialPort=13;
% Exp.Intan.WheelPort=14; 
% Exp.Intan.SpoutPort=1;
% % Exp.Intan.Solenoid;  %non-existent, manauly make an artificial HIGH 2.2s from trialstart
% Exp.Intan.HitPort=16;
% Exp.Intan.FAPort=15;

%Whisker DLC
Exp.DLC.WkrN=2;   %total number inc. all cameras
Exp.DLC.LblPerW=4;

%Stimuli (light)
Exp.Stim.Optic.num=0;
Exp.Stim.Optic.type='BlueLBC';

% Other known parameters
Exp.Intan.SpoutPort_active=0;
Exp.mp4Fps = 30;
Exp.PistonBuffer=0.2;% Mechanical delay of piston shoot-out time following intan signal onset (empirically verified) 0.1s start +0.1s fully out

%%  Read and extract digital signals
[Exp,Sigs]=read_Intan_Signals(Exp);
[Exp,Sigs]=IntanChannels(Exp,Sigs,frequency_parameters,board_adc_channels);
clear frequency_parameters board_adc_channels
[Exp,Sigs]=DigAnalysis(Exp,Sigs,1);   %Exp,Sigs,isplot


% (debug) if mp4N~=intanTrialN
% startvid=7;
% Exp=sort_mp4(startvid,Exp)

%% Analog Input Config
Exp.Stim.Optic.num=1; %0 if there is no light trials
Exp=getOpticsInfo(Exp,Sigs);
%-----------
%for experiments with inconsistent freq
% temp=Exp.Stim.Optic.f;
% temp(temp>=15 & temp<=31)=20;
% Exp.Stim.Optic.f=temp;
%-------If lick piezo is analog and there is no optic, assign opticanalog
%to lickpiezo
Exp.Stim.Optic.num=0;
Sigs.LickPiezo=Sigs.Optic;
Exp=getOpticsInfo(Exp,Sigs);
[Exp,Sigs]=DigAnalysis(Exp,Sigs,1);
%-------------
% If no analog recorded, salvage optic info from command window log (if available)
%create an excel with the log
%-------------
% pathfile='E:\Darren\DarrenM3LMC_210319\DarrenM3LMCrecording_210319_164440\light_log.xlsx';
% Exp=salvage_analog_from_script(Exp,pathfile);
%% Piston/deeplabcut signals Config
close all
Exp.Stim.Piston.active=logical([1 1 0 0]);  %if any pistons are not used, define here
Exp.Stim.Piston.tot=length(Exp.Intan.PistonPort);
Exp.Stim.Piston.num=sum(Exp.Stim.Piston.active);
Exp=getPistonInfo(Exp,Sigs);  %find piston out/in/delay

frame=extractframe(Exp,1,100,10,0);
Exp.DLC.size=size(frame);
figure,imshow(frame),axis on

%if you cropped the video for DLC analysis, specify here
Exp.DLC.iscrop=false;    
Exp.DLC.crop_x1=150;
Exp.DLC.crop_y1=150;    

Exp=findPistonParam(Exp,[0.5 100]);  %find piston .Cam .PWOI .TLOI 
%if no single whisker touch trial exists, you may want to enter them
%manually below by checking Piston.mat
%     Exp.Stim.Piston.Cam=[1 1 nan nan];
%     Exp.Stim.Piston.PWOI=[1 2 nan nan];
%     Exp.Stim.Piston.TLOI=[1.5 5.5 nan nan];
%  
%     Exp.Stim.Piston.Cam=[1 1 nan nan];
%     Exp.Stim.Piston.PWOI=[2 1 nan nan];
%     Exp.Stim.Piston.TLOI=[5 1 nan nan];
%% Save
Exp.savenote='Exp/Sigs';
variableList = {'Exp','Sigs'};
sfilename=['Preprocess_backup_', datestr(floor(now)) '.mat'];
filename=fullfile(Exp.Path.save,sfilename);
disp('Saving...')
for i = 1:length(variableList)
    disp(i) 
    if exist(variableList{i},'var')         
         if exist(filename,'file')
            save(filename,variableList{i},'-append')
         else
             save(filename,variableList{i})
         end
        
    end
end
whos('-file',filename)
disp('Saved')
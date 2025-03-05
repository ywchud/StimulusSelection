%%run the read_Intan_RHD2000 first
function [Exp,Sigs]=IntanChannels(Exp,Sigs,frequency_parameters,board_adc_channels)
%%
try
    path=Exp.Path.data;
catch
    path=[];
end
%Timestamp data file: time.dat

if isempty(path)

    [file, path, filterindex] = ...
        uigetfile('*.dat', 'Select a dat Time File', 'MultiSelect', 'off');
    if (file == 0)
        return;
    end

else
    file='\time.dat';
end
filename = [path,file];
try
        disp('extracting time.dat...')
fid = fopen(filename, 'r');
fileinfo = dir(filename);
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
t = fread(fid, num_samples, 'int32');
fclose(fid);
t = t / frequency_parameters.amplifier_sample_rate; % sample rate from header file
Exp.Fs = 1/(t(2)-t(1));
Exp.t=t;
disp('created t')
catch
    disp('time.dat not found, skipping file')
    Exp.Fs =[];
    Exp.t =t;
end
%%
%Board ADC input data file: analogin.dat
if isempty(path)
[file, path, filterindex] = ...
    uigetfile('*.dat', 'Select a board ADC dat file', 'MultiSelect', 'off');
if (file == 0)
    return;
end

else
    file='\analogin.dat';
end
filename = [path,file];
num_channels = length(board_adc_channels); % ADC input info from header file

try
fileinfo = dir(filename);
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
disp('extracting analogin.dat...')
fid = fopen(filename, 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
fclose(fid);
v = v * 0.000050354; % convert to volts
if ~isvector(v)
    disp('More than 1 optic analog channel found. Assuming optics in 1st and lick 2nd. Change codes in IntanChannels.m as needed')
    if size(v,1)>size(v,2)
        Sigs.Optic=Signal(1,[],v(:,1),Exp.Fs);
        Sigs.LickPiezo=Signal(1,[],v(:,2),Exp.Fs);
    else
        Sigs.Optic=Signal(1,[],v(1,:),Exp.Fs);
        Sigs.LickPiezo=Signal(1,[],v(2,:),Exp.Fs);
    end
else
    Sigs.Optic=Signal(1,[],v,Exp.Fs);
end
Exp.Stim.Optic.num=1;
fprintf('created Optic signal \nExp.Stim.Optic.num set to 1\n')
catch
    disp('analogin.dat not found, skipping file')
    Sigs.Optic=Signal(0,[],[],[]);
    Exp.Stim.Optic.num=0;
end
%%
%Board digital input data file: digitalin.dat
if isempty(path)
    [file, path, filterindex] = ...
        uigetfile('*.dat', 'Select a digital-in dat file', 'MultiSelect', 'off');
    if (file == 0)
        return;
    end

else
    file='\digitalin.dat';
end
try
filename = [path,file];
fileinfo = dir(filename);
num_samples = fileinfo.bytes/2; % uint16 = 2 bytes
disp('extracting digitalin.dat...')
fid = fopen(filename, 'r');
digital_word = fread(fid, num_samples, 'uint16');
fclose(fid);

% FN=fieldnames(Sigs);
disp('decoding digital_word...')
Sigs.Camera=ExtractDig(digital_word,Sigs.Camera);
Sigs.Piston=ExtractDig(digital_word,Sigs.Piston);
Sigs.Trial=ExtractDig(digital_word,Sigs.Trial);
Sigs.Wheel=ExtractDig(digital_word,Sigs.Wheel);
% Sigs.Test=ExtractDig(digital_word,Sigs.Test);
try
    Sigs.Solenoid=ExtractDig(digital_word,Sigs.Solenoid);
catch
    disp('No Solenoid related signals')
end
try
    Sigs.Spout=ExtractDig(digital_word,Sigs.Spout);
catch
    disp('No spout related signals')
end
if isfield(Sigs,'MissIntan')
    Sigs.MissIntan=ExtractDig(digital_word,Sigs.MissIntan);
elseif isfield(Sigs,'HitIntan')
    Sigs.HitIntan=ExtractDig(digital_word,Sigs.HitIntan);
end
Sigs.FAIntan=ExtractDig(digital_word,Sigs.FAIntan);

disp('Created digital signals')

catch
    disp('digitalin.dat not found, skipping file')
end

function object=ExtractDig(digital_word,object)
    for i=1:length(object)
        N=object(i).Port;
        object(i).sig = bitand(digital_word,2^(N-1)) > 0;
        object(i)=object(i).typecheck;
    end
end
%%
%amplifier data in microvolts
% %added 03.06.19 Shulan

% num_channels = length(amplifier_channels); % amplifier channel info from header file
% [file, path, filterindex] = ...
%     uigetfile('*.dat', 'Select the amplifier.dat', 'MultiSelect', 'off');
% if (file == 0)
%     return;
% end
% filename = [path,file];
% s = dir(filename);
% num_samples = s.bytes/(num_channels * 2); % int16 = 2 bytes
% clear s;
% fid = fopen(filename, 'r');
% RawAnalog = fread(fid, [num_channels, num_samples], '*int16');
% fclose(fid);
% RawAnalog = RawAnalog * 0.195; % convert to microvolts
% for i = 1:num_channels % downsample to 1k
%     FieldPotential(i,:) = resample(double(RawAnalog(i,:))*0.0195,1,30);
% end
% t_FP=resample(t,1,30);
% clear RawAnalog;

end
function LFP=getLFPbyChannel(path,channel_num,fs,fc,p)
% fs=1000;
% fc=200;
%downsampling params
% p=1000;
num_channels = channel_num; 
file='\amplifier.dat';
filename = [path,file];
s = dir(filename);
num_samples = s.bytes/(num_channels * 2); % int16 = 2 bytes
clear s;

disp('Spliting raw amplifier data by channels...')
tic
%file splitter
FID = fopen(filename, 'r');
FID_List = cell(num_channels, 1);
% RawAnalog = fread(fid, [num_channels, num_samples], '*int16');
ChPath=fullfile(path,'AmpPerCh');
mkdir(ChPath)
for j = 1:num_channels
    FID_List{j} = fopen(fullfile(ChPath,sprintf('CH_%d.dat', j)), 'w');
end
Data = fread(FID, [num_channels num_samples], '*int16');
for k = 1:num_channels
    fwrite(FID_List{k}, Data(k, :), 'int16');
end
fclose all;
disp('Stored in folder AmpPerCh')
toc

%low pass params
q=fs;  %original sampling frequency
[b,a] = butter(6,fc/(q/2)) ;
% figure
% freqz(b,a)
disp('Obtaining LFP...')
LFP=zeros(num_channels,ceil(num_samples*(p/q)));
for k = 1:num_channels
    tic
    fid=fopen(fullfile(ChPath,sprintf('CH_%d.dat', k)), 'r');
    RawAnalog = fread(fid, [1 num_samples], '*int16');
    RawAnalog=RawAnalog * 0.195; % convert to microvolts
    fclose(fid);
    
    LFP(k,:)=resample(filter(b,a,RawAnalog),p,q);
    fprintf('Channel %d/%d done\n',k,num_channels)
    toc
end

savepath=fullfile(path,'LFP.mat');
 save(savepath,'LFP')
end



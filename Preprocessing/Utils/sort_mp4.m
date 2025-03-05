function Exp=sort_mp4(startvid,Exp)
% say u did the experiment around midnight so the video names change from
% 23-XX-XX to 00-XX-XX then u somehow sorted the videos by names and
% converted mp4 using that order, u messed up the actual order and u
% decided to create a new folder with the correct order of videos
% thats what we r doing here

%input:
%startvid:  the corresponding actual video number of the first video in the
%current wrongly ordered video, e.g. if startvid=83, reorder vids by
%83,84..end,1,2...82


% startvid=83;
if length(Exp.Path.vid)>1
    error('Checkpoint. Only for one camera folder now, adjust code as needed')
end

oldpath=Exp.Path.vid{:};
newpath=append(oldpath,'_reordered');
if exist(newpath, 'dir')
    error('Delete reordered vid folder first')
end
mkdir(newpath);

vidname=Exp.Path.vidName{:};
csvname=Exp.Path.csvName{:};

temp=dir(fullfile(oldpath,'*.mp4'));

if Exp.TrN~=length(temp)
    fprintf('Mismatch mp4(%d) and signal(%d), del additional files manually\n',length(temp),Exp.TrN)
    mp4N=length(temp);
else
    mp4N=Exp.TrN;
end

reorderedIDs=[startvid:mp4N 1:startvid-1];


for k = 1 : length(temp)
  thisFileName = sprintf(vidname,reorderedIDs(k));
  % Prepare the input filename.
  inputFullFileName = fullfile(oldpath, thisFileName);
  % Prepare the output filename.
  outputBaseFileName = sprintf(vidname,k);
  outputFullFileName = fullfile(newpath, outputBaseFileName);
  % Do the copying and renaming all at once.
  copyfile(inputFullFileName, outputFullFileName);
  fprintf('%d/%d\n',k,length(temp))
end
for k = 1 : length(temp)
  thisFileName = sprintf(csvname,reorderedIDs(k));
  % Prepare the input filename.
  inputFullFileName = fullfile(oldpath, thisFileName);
  % Prepare the output filename.
  outputBaseFileName = sprintf(csvname,k);
  outputFullFileName = fullfile(newpath, outputBaseFileName);
  % Do the copying and renaming all at once.
  copyfile(inputFullFileName, outputFullFileName);
  fprintf('%d/%d\n',k,length(temp))
end
Exp.Path.vid{:}=newpath;

end
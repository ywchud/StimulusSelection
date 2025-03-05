function Exp=get_video_files(Exp,codeDir)
%% find camera and dlc filename
for i=1:Exp.Cam.num
    cd(Exp.Path.vid{i})

    temp=dir(fullfile('*.mp4'));
    Exp.TrN=length(temp);
    temp=temp.name;
    if isempty(temp)
        error('Video not found')
    end
    while temp(end)~='_'   %remove everything after the last '_'
        temp=temp(1:end-1);
    end
    Exp.Path.vidName{i}=strcat(temp,'%d.mp4');
    
    temp=dir(fullfile('*_filtered.csv'));
    if isempty(temp)
        disp('No _filtered csv found, using all other csvs in the folder')
        temp=dir(fullfile('*.csv'));
    else
        disp('Using _filtered csv');
    end
    if length(temp)~=Exp.TrN
        error('Number of .csv(%d) and .mp4(%d) mismatch',length(temp),Exp.TrN)
    end
    temp=temp.name;
    if isempty(temp)
        error('DLC .csv not found')
    end
    k = strfind(temp,'DLC');
    while temp(k-1)~='_'   %remove everything between '_' and DLC
        temp(k-1)=[];
        k=k-1;
    end
    Exp.Path.csvName{i}=insertAfter(temp,k-1,'%d');   
    clear temp
end

cd(codeDir)
end
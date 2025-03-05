function frame=extractframe(Exp,CamN,VidN,Frame,withLabel)
    if length(CamN)>1 && isvector(CamN)
        error('Choose only 1 cam')
    end
    if length(VidN)>1 && isvector(VidN)
        error('Choose only 1 vid')
    end
    i=1;
    
    datapath=Exp.Path.vid{CamN(i)};
    Labelpath=fullfile(datapath,'Labelled');
    videoname = sprintf(Exp.Path.vidName{CamN(i)}(1:end-4),VidN);  
    if withLabel
        try
         V{i} = VideoReader(fullfile(Labelpath,[videoname '_labeled.mp4']));
        catch
         dataname = sprintf(Exp.Path.csvName{CamN(i)},VidN);
         write_labed_video(datapath,dataname,[videoname '.mp4']);  
         V{i} = VideoReader(fullfile(Labelpath,[videoname '_labeled.mp4']));
        end    
    else
        V{i} = VideoReader(fullfile(datapath,[videoname '.mp4']));
    end
    Width=V{i}.Width;
    Height=V{i}.Height;
    frame=zeros(Height,Width,3,length(Frame));
    
    for k=1:length(Frame)
        frame(:,:,:,k) = read(V{i},Frame(k));
    end
    frame=uint8(frame);
end
function deets(Exp)
    path=dir(fullfile(Exp.Path.data,'deet*.txt'));
    winopen(fullfile(path.folder,path.name));
end
%% CSD prep  (make sure LFP is in the workspace/folder)
codeDir='Your_path';
cd(codeDir);
try
    load(fullfile(Exp.Path.data,'LFP.mat'))
catch
    error('No LFP.mat found in folder')
end
Exp.KSort.LFP_fs=length(LFP)/length(Exp.t)*Exp.Fs;
CSD_data=CSD_prep(Exp,LFP,[-0.02 0.06]);
CSDplot(Exp,Neuron,CSD_data)
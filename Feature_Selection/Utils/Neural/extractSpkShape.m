function Neuron=extractSpkShape(Exp,Neuron)
Path=Exp.Path.data;
%% generate results to .csv from python code   
%DOESNT WORK, NEED more recent MATLAB version

% system('python SpikeFormSortingPY.py Path')
% pyrunfile("E:\DarrenGen2\borrowed_codes\SpikeFormSortingPY.py", Path)
%     cd('E:\DarrenGen2\borrowed_codes');
%     pyversion C:\Users\PLUTALAB\AppData\Local\Continuum\anaconda3\python.exe
%     py.importlib.import_module('SpikeFormSortingPY')
% pythonFolderPath='C:\Users\PLUTALAB\AppData\Local\Continuum\anaconda3';
% cd(pythonFolderPath)
% addpath('C:\Users\PLUTALAB\AppData\Local\Continuum\anaconda3\DLLs');
% systemCommand = ['C:\Users\PLUTALAB\AppData\Local\Continuum\anaconda3\python.exe E:\DarrenGen2\borrowed_codes\SpikeFormSortingPY.py ',Path];
% system(systemCommand)

%% get results from .csv
try
    fs = readmatrix(fullfile(Path,'cluster_fs.csv'));
    rs = readmatrix(fullfile(Path,'cluster_rs.csv'));
    un = readmatrix(fullfile(Path,'cluster_un.csv'));

catch
    systemCommand=['python E:\DarrenGen2\borrowed_codes\SpikeFormSortingPY.py ' Path];
    error('Csv files not found. Run the following line in anaconda first\n%s',systemCommand)
    
end

% Templates = readNPY(fullfile(Path,'\templates.npy'));
SpikeClusters = readNPY(fullfile(Path,'\spike_clusters.npy'));
Spike_templates = readNPY(fullfile(Path,'\spike_templates.npy'));

tempPerClu = findTempForEachClu(SpikeClusters, Spike_templates);
CID=unique([Neuron.CID]);
RepTemp=tempPerClu(CID+1);


for i=1:length(Neuron)
    
    if ismember(RepTemp(i),fs)
        Neuron(i).Shape="fs";
    elseif ismember(RepTemp(i),rs)
        Neuron(i).Shape="rs";
    else
        Neuron(i).Shape="un";
    end

end
 Depths=[Neuron.Depth]';
binrng= min([Neuron.Depth]):Exp.KSort.ChGap:max([Neuron.Depth]);
binWidth=mode(diff(binrng));
counts1 = histcounts(Depths([Neuron.Shape]=="fs"),binrng);                                   
counts2 = histcounts(Depths([Neuron.Shape]=="rs"),binrng);   
counts3 = histcounts(Depths([Neuron.Shape]=="un"),binrng);

barValues=[counts1;counts2;counts3]';
figure
% barh(binrng(1:end-1),barValues,'histc')
barh(binrng(1:end-1)+binWidth/2,barValues,'stacked')

ylabel('Depth below pia mater (um)')
 xlabel('NeuronN')
legend({'fs','rs','un'})

filename=fullfile(Exp.Path.save,'NeuronShape2Depth.jpg');
saveas(gcf,filename)
% Depths=repmat([Neuron.Depth]',[1,size(RES,2)]);
% figure, histogram(Depths([Neuron.Shape]=="fs"),min([Neuron.Depth]):Exp.KSort.ChGap:max([Neuron.Depth]), 'Orientation', 'horizontal'),hold on
%  histogram(Depths([Neuron.Shape]=="rs"),min([Neuron.Depth]):Exp.KSort.ChGap:max([Neuron.Depth]), 'Orientation', 'horizontal')
%  histogram(Depths([Neuron.Shape]=="un"),min([Neuron.Depth]):Exp.KSort.ChGap:max([Neuron.Depth]), 'Orientation', 'horizontal')
%  
%  ylabel('Depth below pia mater (um)')
%  xlabel('N Observation')
%  legend({'fs','rs','un'})


end
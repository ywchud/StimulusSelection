function y=getWinMean(data,win,avgObs)
%data: 2d matrix or cell(containing observations) with dim1 NeuN and dim2 win length
%win: binary vector
%avgObs: bool (0/1) average all observations
% takes mean value of data within win and return vector of length N. if data is
% a cell type with observations, it averages the obs first per cell

if length(size(data))~=2 || length(win)~=size(data,2)
    error('Wrong arg dim')
end  

if ~exist('avgObs','var') 
    avgObs=1;
end
    

temp=data(:,win);
if iscell(temp)
    if avgObs %average observations within each cell
        y1=cellfun(@(x) sum(x(~isnan(x))),temp);
        L1=cellfun(@(x) sum(~isnan(x)),temp);
        y2=0;L2=0; %if want to merge other observations from another cell
        mat=(y1+y2)./(L1+L2);
        y=mean(mat,2); %avg across windows
    else %keep observations within each cell
        obsN=sum(squeeze(cellfun(@(x) length(squeeze(x)),data(:,1))),2);
        y=cell(size(data,1),1);
        for i=1:size(data,1)
            mat=zeros(obsN(i),sum(win));
            for j=1:sum(win)
                mat(:,j)=cell2mat(squeeze(temp(i,j)));
            end
            y{i,1}=mean(mat,2); %avg across windows
        end
        
    end
else
    mat=temp;
    y=mean(mat,2); %avg across windows
end

end
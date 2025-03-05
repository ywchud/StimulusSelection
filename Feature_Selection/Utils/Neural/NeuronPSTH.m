function res=NeuronPSTH(Exp,Whisker,Neuron,NeuronIdx,StimCi,fp)
%Note: This function is the final umbrella of subfunctions requiring many preprocessing of stimulus info,
%i.e. Exp,Whisker,StimCi

%This function gives PSTH results, including 
%PSTH mu, std, peaks, latency, fwhm, etc, retaining the raw PSTHmat for
%hypothesis testing if necessary
%Output: res - a structure matrix of all stimulus triggered results, dimension is
%Neuron*StimCi*fp
% Input: 
% Exp,Whisker,StimCi for finding trial,piston,whisker,touch related times
% Neuron, neuronclass to be analysed
% NeuronIdx, id of Neuron to be analysed, default 1:length(Neuron)
% fp, parameters of PSTH, such as range, aligntype, and filters for said
% aligntype. dim - array of fp(s)
%vvvvv e.g. vvv
% fp.firstfew=[];
% fp.ActiveTrialsOnly = 1;
% fp.Sequence = [];
% fp.ITIth=0.03;
% fp.isrun=[];
% fp.Tdur=[];
% fp.bin = 0.005;
% fp.PSTHrange = (-0.05:p.bin:0.05);
% fp.windowtime=[0 0.04];     %within PSTHrange, for finding mu,std
% fp.isrun=1;
% fp.type='Touch'       %check getTimePoints.m for a list of options

%added 08/07
% fp.Toverlap
% fp.mainW

% fp.Tgap (to be done)
if isvector(StimCi)
    StimCi=StimCi(:)';  %1*n vector
end

for i=1:size(StimCi,2)
    fprintf('%d/%d\n',i,size(StimCi,2))
    S=StimCi{1,i};
    for j=1:length(NeuronIdx)
        
        [S_adj]=S.PWadjusted(Exp,Whisker,Neuron(NeuronIdx(j)));   %adjustPW        
        
%         test1=updateStimC(S,Exp,fp(1));
%         test2=updateStimC(S_adj,Exp,fp(1));
%         figure
%         histogram(test1.TouchDeltaT,'binedge',-0.04:0.004:0.04)
%         hold on
%         histogram(test2.TouchDeltaT,'binedge',-0.04:0.004:0.04)

        
        if strcmp(fp.mainW,'FW')
            [S_adj]=S_adj.FWadjusted(Exp);  %adjust first whisker touch
        end
        
        if ~isempty(fp)
            StimCs=cell(1,length(fp));
            fpp=fp;  %replicated fp for making adjustment to free-whisking param
            for k=1:length(fp)
                if isempty(S_adj.WID)  %free whisking no touch
                    StimCs{k}=S_adj;   %avoid touch updating or there will be no validtrials
                    if strcmp(fpp(k).type,'Retraction_Tonly')   %replace touch related features as needed
                        fpp(k).type='Retraction_sampled';
                    elseif strcmp(fpp(k).type,'Protraction_Tonly')
                        fpp(k).type='Protraction_sampled';
                    end
                else
                    StimCs{k}=updateStimC(S_adj,Exp,fp(k));  %touch screening
                end
            end
        else
            StimCs{1}=S_adj;
        end
        for k=1:length(fpp)
%             outshellID=(i-1)*length(StimCs)+k;
            [temp,Snew]=findPSTH(StimCs(k),Exp,Whisker,Neuron,NeuronIdx(j),fpp(k));   %calculate PSTH for all fp
            temp.TouchDeltaT=Snew{1}.TouchDeltaT;
            temp.ReleaseDeltaT=Snew{1}.ReleaseDeltaT;
            temp.TouchTimes=Snew{1}.TouchTimes;
            temp.ReleaseTimes=Snew{1}.ReleaseTimes;          
            temp.WID=Snew{1}.WID;
            temp.WID_FW=Snew{1}.WID_FW;
            temp.TouchDur=Snew{1}.ReleaseTimes-Snew{1}.TouchTimes;
            temp.Toverlap=Snew{1}.Toverlap;
            temp.stat=[]; %placeholder
            
%             if length(temp.CorrTrial{1})~=length(temp.TouchTimes)
%                 
% %                 temp.CorrTrial_touch=Snew{1, 1}.CorrTrials;
%             end
            
            if strcmp(fpp.type,'Touch') && length(temp.TouchTimes)~=size(temp.PSTHmat{1},1)
                error('Touch number obtained is different from sampled number, check buggy code')
            end
            
            res(j,i,k)=temp;
        end
        
    end
end
end
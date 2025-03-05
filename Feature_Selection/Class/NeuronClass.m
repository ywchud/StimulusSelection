classdef NeuronClass
     properties( Access = public )
         % All your arguments here with their default values
         ExpID;  %to link to the corresponding exp parmas, imagine an array of neurons combined from all experiments with subset of them linking to a particular exp array
         ID=nan;
         CID=nan;
         SpikeTime=nan;  %1D 
         SpikeCount=nan;
         Shape='un';
         ShapeForm=nan;
         XX=nan;
         Shank=nan;
         ShankPW=nan;
         Depth=nan;
         Depth_b4correction=nan;
         depthType=nan;  % 0=superficial 1=deep for now, can expand
         isGood=nan;
         
         pre_post_Touch_Corr=nan;
         TrBLF_TrTouchAvg_Corr=nan;
         
              
        touchMU=nan;  
        pretouchMU=nan;
        peak=nan;
        lat=nan;
        sign=nan;
        BLF=nan;
        fwhm=nan;
        deltahm=nan;
        x1=nan;
        x2=nan;

        firstSpkLat=nan;
        firstSpkJitter=nan;
        firstSpkTendency=nan;
        
        R=nan;   %Store RES_custom
        PARAMS=nan;   %params used to get R
        Rstat=nan;  %anova within ctrl(noLight) stimcs
        Rstat_wilcoxon; %wilcoxon ranksum of ctrl vs (light/postlight)
         
         PW=nan;  %all whisker activity sorted
         PW_H=nan;  %sig trial baseline vs trial piston out
         PW_H_Touch=nan;  %sig PSTH before 0 vs PSTH after zero
         PW_H_SpkShape=nan; %sig spk time before 0 vs spk time after zero
         H_light_rate=nan;
         H_light_firstspk=nan;
         
         pw=nan;  %principle whisker, i.e. PW(1)
         ActiveTrials=nan; %Trial indexes when blf is > ActiveTh
        ActiveTh=nan;  %median based outliers
        
     end
     methods
         function obj = NeuronClass(Exp,NeuronID,ClusterID,SpikeTime,XX,YY,shank,isGood)
            obj.ExpID = Exp.ID;
            obj.ID    = NeuronID;
            obj.CID    = ClusterID;
            obj.SpikeTime=SpikeTime;   
            obj.XX=XX; 
            obj.Depth=YY;   
            obj.Shank=shank;
            obj.isGood=isGood;
            obj.SpikeCount=length(SpikeTime);
         end
         function str=Info(obj)
%              str=sprintf('ID%dC%d %s pw:%d shkpw:%d %dum',obj.ID,obj.CID,obj.Shape,obj.pw,obj.ShankPW,obj.Depth);
             str=sprintf('C%d %s pw:%d shkpw:%d %dum',obj.CID,obj.Shape,obj.pw,obj.ShankPW,obj.Depth);
         end
         
         function SpikeTInT = SpikeTimeInTrials(obj,Exp)  %segregate SpikeTime into cells of trials. Each cell is an array of spiketime in that trial.
             SpikeTInT=cell(Exp.TrN,1);
             for i=1:Exp.TrN
                SpikeTInT{i,1}(:,1)=obj.SpikeTime(obj.SpikeTime>Exp.FrameT{i}(1) && obj.SpikeTime<Exp.FrameT{i}(end));
             end
         end
            
         function [SpikeBInT,binTime] = SpikeBinInTrials(obj,Exp,Fs)  %bin SpikeTime into cells of trials. Each cell is an array of binned spikeNumber with length of Trial Frame No.(default 500Fps)             
             SpikeBInT=cell(Exp.TrN,1);  %number of spikes per bin in each trial
             binTime=cell(Exp.TrN,1);   %time of bin in each trial
             for i=1:Exp.TrN
                 
                 if isempty(Fs)         %(default Fs=videoFps=500Fps)
                     FN=Exp.FrameN(i);
                     FT=Exp.FrameT{i};
                 else
                     FN=(Exp.FrameT{i}(end)-Exp.FrameT{i}(1))*Fs;
                     FT=linspace(Exp.FrameT{i}(1),Exp.FrameT{i}(end),FN);
                 end
                 
                SpikeTinT=obj.SpikeTime(obj.SpikeTime>FT(1) & obj.SpikeTime<FT(end));
                Nspikes=length(SpikeTinT);
                SpikeMat=zeros(FN,1);
                for j=1:Nspikes
                    spikeFrame=find(FT>=SpikeTinT(j),1,'first')-1;
                    SpikeMat(spikeFrame)=SpikeMat(spikeFrame)+1;
                end
                SpikeBInT{i,1}(:,1)=SpikeMat;
                binTime{i,1}(:,1)=FT;
             end
         end
         
         function obj=findActiveTrials(obj,Exp,isplot)
             
             th=0.05;  %0-1
             MAwindow=0.05;  %sec
             if ~exist('isplot','var')
                 isplot=0;
             end
             [N,edges]=histcounts(obj.SpikeTime,1:Exp.t(end)/10000:Exp.t(end));
             t=edges(1:end-1);
             N=smooth(N,Exp.Fs*MAwindow);
             if quantile(N,0.95)~=0
                 maxV=quantile(N,0.95);
                 temp=rescale(N,'InputMax',maxV); %remove outliers
             else
                 maxV=max(N);
                 temp=rescale(N);
             end

             t_first=t(find(temp>=th,1,'first'));
             t_last=t(find(temp>=th,1,'last'));
             
             if isplot
                 plot(t,temp)
                 hold on
                 yline(th);
                 temp2=temp;
                 temp2(temp>=th)=th;
                 area(t,temp2)
                 xline(t_first,'g');
                 xline(t_last,'r');
                 yyaxis right
                 histogram(obj.SpikeTime,1:Exp.t(end)/10000:Exp.t(end));
             end
             
             obj.ActiveTrials=Exp.TrialEndT<t_last & Exp.TrialStartT>=t_first; %Trial indexes when blf is > ActiveTh
             obj.ActiveTh=maxV*th;
         end
         
         
         function [PSTHavg, SEM, SumStat]=PSTH(obj,timepoints,PSTHrange,windowtime)
             %output has a dimension of length(PSTHrange)-1
             %use PSTHrange(2:end)-bin/2 for plotting x axis
             if length(PSTHrange)<3
                 error('PSTHrange shoulf be a vector of size >2')
             end
             if isempty(windowtime)
                windowtime=[PSTHrange(1) PSTHrange(end)];  %for finding mu and sd
             end
%              windowtime=[0 PSTHrange(end)];
             bin=diff(PSTHrange(1:2));
             PSTHavg = zeros(1,length(PSTHrange)-1);
             PSTHmat=zeros(length(timepoints),length(PSTHrange)-1);
             PSTHrelT=cell(length(timepoints),1);
             spikeCounts=zeros(length(timepoints),1);
             FirstSpikeLat=zeros(length(timepoints),1);
             zeroIndex=find(PSTHrange>=windowtime(1),1,'first');
             lastIndex=find(PSTHrange>=windowtime(2),1,'first')-1;
             
             for j = 1:length(timepoints)
                 [PSTHmat(j,:),~] =  histcounts(toColumn(obj.SpikeTime),PSTHrange + timepoints(j));
                 PSTHavg = PSTHavg + PSTHmat(j,:);
                 spikeCounts(j)=sum(PSTHmat(j,zeroIndex:lastIndex));
                 FirstSpikeLat_I=find(obj.SpikeTime>=timepoints(j),1,'first');
                 FirstSpikeLat(j)=nan;
                 if ~isempty(FirstSpikeLat_I)
                    T=obj.SpikeTime(FirstSpikeLat_I)-timepoints(j); 
                    if T<=windowtime(2)  %makes sure first spike isnt too far late
                        FirstSpikeLat(j)=T;
                    end
                 end
             end
             
            
             for j = 1:length(timepoints)
                 relTime=obj.SpikeTime-timepoints(j);
                 PSTHrelT{j,1}=relTime(relTime<PSTHrange(end) & relTime>=PSTHrange(1));
             end
             %debug    
             floatTolerance=5.0000e-06; 
             for j = 1:length(timepoints)
                 if sum(PSTHmat(j,:))>length(PSTHrelT{j,1})   %float error      
                     fT=floatTolerance; %included into bounds due to float error when doing deduction
                     relTime=obj.SpikeTime-timepoints(j);
                     PSTHrelT{j,1}=relTime(relTime<PSTHrange(end)+fT & relTime>=PSTHrange(1)-fT);      
                     if sum(PSTHmat(j,:))~=length(PSTHrelT{j,1})
                         123;
                     end
                 elseif sum(PSTHmat(j,:))<length(PSTHrelT{j,1})
                     123;
                 end
             end
             
             SEM=std(PSTHmat/bin,0,'omitnan')/(size(PSTHmat,1))^0.5;
             PSTHavg = PSTHavg/length(timepoints)/bin;             
             SumStat.firstSpkLats=FirstSpikeLat;
             SumStat.firstSpkLat=mean(FirstSpikeLat,'omitnan');
             SumStat.firstSpkJitter=std(FirstSpikeLat,'omitnan');               
             SumStat.firstSpkTendency=sum(FirstSpikeLat<=windowtime(2))/length(FirstSpikeLat);
%              [peak,I]=max(abs(PSTHavg));
%                 latency=(I-zeroIndex)*bin;
%               sign=(PSTHavg(I)<0)*-2+1;

            baseline=mean(PSTHavg(1:zeroIndex));
            baseline_subtracted=PSTHavg(zeroIndex:lastIndex)-baseline;
            [peak_temp,I]=max(abs(baseline_subtracted));
            pksign=(baseline_subtracted(I)<0)*-2+1;
             SumStat.latency=PSTHrange(I+zeroIndex-1)+bin/2;
             SumStat.abs_peak=pksign*peak_temp+baseline;
             SumStat.sign=pksign;
             
             SumStat.peak=max(PSTHavg(zeroIndex:lastIndex));
             SumStat.trough=min(PSTHavg(zeroIndex:lastIndex));
             
             [fwhm_temp,delta_halfHeight,x1_temp,x2_temp]=FWHM(PSTHavg(zeroIndex:lastIndex),PSTHrange(zeroIndex:lastIndex));
             SumStat.fwhm=fwhm_temp;
             SumStat.deltahm=delta_halfHeight;
             SumStat.x1=x1_temp;
             SumStat.x2=x2_temp;
             
             
             spikeRate=spikeCounts/(windowtime(2)-windowtime(1));
             
             SumStat.mu=mean(spikeRate);
             SumStat.preZero=baseline;   %calculated from first possible index to windowtime(1)
             SumStat.sd=std(spikeRate);
             SumStat.PSTHmat=PSTHmat;
             SumStat.PSTHrelT=PSTHrelT;
             SumStat.zeroIndex=zeroIndex;
             SumStat.lastIndex=lastIndex;      
             
%              if SumStat.firstSpkLat-SumStat.latency>0.003 && windowtime(2)==0.04 && SumStat.firstSpkTendency>=0.8 && any(obj.PW_H==1)
%                  figure
%                  plot(PSTHrange(1:end-1),PSTHavg)
%                  figure
%                  histogram(FirstSpikeLat,0:0.005:0.04)
%                  SumStat.firstSpkLat
%                  SumStat.latency
%                  SumStat.firstSpkTendency   
%              end
             
         end
         
         function obj=stats(obj,Rbool,Rcolor,type,stat2do,p)
             
             %stat2do
%              1: BLF vs TAF  (not implemented)
%              2: ttest2, i.e. unpaired ttest of pre vs post touch per stim
%              individually
%              3: anova of all given stim
             
             %p params:
             % sample_window: vector of 2 of sample range,
                %e.g. [-0.04 0.04], will take -0.04:0 vs 0:0.04 if within
                %PSTHrange
             % alpha: test significance alpha
             % tail: 'both'(default) 'right' 'left'
             
             [Rslice,Pslice,Rbool,Rcolor,fp_id]=obj.checkinput(Rbool,Rcolor,type);
             zeroIndex=find(Pslice.PSTHrange==0);
             [row,col]=size(Rslice);  %row 10(stim) col 2(light)
             range=Pslice.PSTHrange;
             
             if zeroIndex<=3
                 error('Fewer than 2 pretouch(t<0) data points. Get new R slice with longer PSTHrange')
             elseif length(range)-zeroIndex<=2
                 error('Fewer than 2 posttouch(t>=0) data points. Get new R slice with longer PSTHrange')
             end
             
             %p check
             if isfield(p,'sample_window') && length(p.sample_window)==2
                 [I1,~]=findClosest(p.sample_window(1),range(1:end-1),[]);
                 [I2,~]=findClosest(p.sample_window(2),range(1:end-1),[]);
                 sample_windowI=[I1 I2];
             else
                 sample_windowI=[1 length(range)-1];
             end
             newrange=range(sample_windowI(1):sample_windowI(2));
             
              if isfield(p,'alpha') && isnumeric(p.alpha)
                  alpha=p.alpha;
              else
                  alpha=0.05;
              end
             if isfield(p,'tail') && (isstring(p.tail) || ischar(p.tail))
                  tail=p.tail;
              else
                  tail='both';
             end
             
             if isfield(p,'method') 
                 method=p.method;
             else
                 method='Wilcoxon'; %'ttest'
             end
             

 
 

             
             outputL=sum(Rbool(:));
             
             
             %1: BLF(b4 piston) vs TAF(after piston)
             % already done in mastercode stored in obj.PW_H
             if ismember(1,stat2do)
                 error('Check obj.PW_H. If not there, runin  mastercode')
             end
             %-------------------------
             
             %data slicing
             preT=cell(1,outputL);
             postT=cell(1,outputL);
             BaselineFr=zeros(1,outputL);
             %              ind = sub2ind(size(Rbool),i,j);
             
             count=0;
             savei=zeros(1,outputL);
             savej=zeros(1,outputL);
             preT_L=length(sample_windowI(1):zeroIndex-1);
             postT_L=length(zeroIndex:sample_windowI(2));
             
             for j=1:col
                 for i=1:row
                     if ~Rbool(i,j)
                         continue;
                     else
                         count=count+1;
                         savei(count)=i;
                         savej(count)=j;
                     end
                     Rij=Rslice(i,j);
                     mat=Rij.PSTHmat{1};
                     preT{count}=mat(:,sample_windowI(1):zeroIndex-1);
                     postT{count}=mat(:,zeroIndex:sample_windowI(2));
                     BaselineFr(count)=Rij.BLF;    %note BLF used here would be different for Ctrl and Lighton
                 end
             end
             
             %--------------------
             %2: pretouch vs posttouch
             % unfiltered touch samples used and stored in PW_H_Touch
             %using filtered samples from obj.R hereon
             %ttest or wilcoxon
             if ismember(2,stat2do)
                 
                 H=zeros(1,outputL);
                 P=zeros(1,outputL);
                 for i=1:outputL
                     x=sum(preT{i},2)/preT_L;  %using rate instead of spike count since pre and post can have different window sizes
                     y=sum(postT{i},2)/postT_L;
                     if isempty(x) || isempty(y)
                         H(i)=nan;P(i)=nan;
                         continue
                     end
                     
                     if strcmp(method,'ttest')
                        [H(i),P(i)] = ttest2(x,y,'Alpha',alpha,'Tail',tail);
                     elseif strcmp(method,'Wilcoxon')
                     	[P(i),H(i)] = ranksum(x,y,'Alpha',alpha,'Tail',tail);        
                     else
                        [P(i),H(i)] = ranksum(x,y,'Alpha',alpha,'Tail',tail);        
                     end

                 end
                 %save into obj.R.stat
                 for k=1:outputL
                     obj.R(savei(k),savej(k),fp_id).stat.TouchSig=H(k);
                     obj.R(savei(k),savej(k),fp_id).stat.TouchP=P(k);
                 end
             end
              %--------------------
             %3: posttouch anova across all stims(n=outputL) given in Rbool
             if ismember(3,stat2do)
                 Responses=[];
                 group=[];
                 for i=1:outputL
                     x=sum(postT{i},2)/postT_L;
                     Responses=[Responses;x(:)];
                     g=i*ones(length(x),1);
                     group=[group;g(:)];
                 end
                 [p,tbl,stats] = anova1(Responses,group,'off');    
                 try
                 [c,m,h,nms] = multcompare(stats,'Alpha',alpha,'Display','off');
                 catch
                     c=nan;
                     nms=nan;
                 end
                 
                 
                 %save into obj.Rstat
                 Rstats.p=p;
                 Rstats.c=c;
                 Rstats.savei=savei;
                 Rstats.savej=savej;
                 Rstats.nms=nms;
                 obj.Rstat=Rstats;
             end
              %--------------------
             %4: pretouch vs constant BLF ttest/wilcoxon
             if ismember(4,stat2do)
                  H=zeros(1,outputL);
                  P=zeros(1,outputL);
                 for i=1:outputL
                     x=sum(preT{i},2)/preT_L;
                     y=BaselineFr(i);
                     if isempty(x) || isempty(y)
                         H(i)=nan;P(i)=nan;
                         continue
                     end
                    if strcmp(method,'ttest')
                        [H(i),P(i)] = ttest2(x,y,'Alpha',alpha,'Tail',tail);
                     elseif strcmp(method,'Wilcoxon')
                     	[P(i),H(i)] = ranksum(x,y,'Alpha',alpha,'Tail',tail);        
                     else
                        [P(i),H(i)] = ranksum(x,y,'Alpha',alpha,'Tail',tail);        
                     end    

                 end
                 %save into obj.R.stat
                 for k=1:outputL
                     obj.R(savei(k),savej(k),fp_id).stat.PistonSig=H(k);
                 end  
             end
             
             %5: wilcoxon rank sum test between two StimC (e.g. for ctrl vs light)
             
             if ismember(5,stat2do)
                if outputL~=2
                    error('Rbool should have length=2')
                end          
                x=sum(postT{1},2)/postT_L;  %using rate instead of spike count since pre and post can have different window sizes
                y=sum(postT{2},2)/postT_L;
                
                if isempty(x) || isempty(y)
                    h=nan;p=nan;
                elseif strcmp(method,'ttest')
                    [h,p] = ttest2(x,y,'Alpha',alpha,'Tail',tail);
                elseif strcmp(method,'Wilcoxon')
                    [p,h] = ranksum(x,y,'Alpha',alpha,'Tail',tail);
                else
                    [p,h] = ranksum(x,y,'Alpha',alpha,'Tail',tail);
                end
   
                if isstruct(obj.Rstat_wilcoxon)
                    Rstat_w=obj.Rstat_wilcoxon;
                end
                Rstat_w(savei(1),savej(2)-1).p=p;
                Rstat_w(savei(1),savej(2)-1).h=h;
                Rstat_w(savei(1),savej(2)-1).Rtype=fp_id;
                
                obj.Rstat_wilcoxon=Rstat_w;
             end
             
             
             
             
%              gs=cellfun(@str2num,nms);
             % AnovaP{i}=c;
%              for j=1:size(c,1)
%                  Pmat(gs(c(j,1)),gs(c(j,2)))=c(j,end);
%              end
             
%              sigpairs=find(c(:,6)<=alpha);
%              notsigpairs=find(~(c(:,6)<=alpha));

%              figure
%              [si,sj]=subplotDim(length(sigpairs));
%              for i=1:length(sigpairs)
%                  subplot(si,sj,i)
%                  arr=zeros(1,10);
%                  arr(c(sigpairs(i),1))=1;arr(c(sigpairs(i),2))=1;
%                  rbool=logical([arr;zeros(1,10)]);
%                  obj.plotPSTH(rbool,[],1);
%              end
             
         end
         function [H,P]=stats_2slices(obj,Rbool,type,p)
             if ~iscell(Rbool) || length(Rbool)~=2
                 error('rbool should be a cell array of length 2')
             end
             if length(type)~=2
                 error('type should be an array of 2')
             end
             for i=1:length(type)
                 [Rslice,Pslice,rbool,~,~]=obj.checkinput(Rbool{i},[],type(i));       
                 Rslices(i)=Rslice(rbool);
             end
             zeroIndex=find(Pslice.PSTHrange==0);
%              [row,col]=size(Rslice);  %row 10(stim) col 2(light)
             range=Pslice.PSTHrange;
             %p check
             if isfield(p,'sample_window') && isvector(p.sample_window) && length(p.sample_window)==2
                 [I1,~]=findClosest(p.sample_window(1),range(1:end-1),[]);
                 [I2,~]=findClosest(p.sample_window(2),range(1:end-1),[]);
                 sample_windowI=[I1 I2;I1 I2];
             elseif isfield(p,'sample_window') && ismatrix(p.sample_window)
                 [I1,~]=findClosest(p.sample_window(1,1),range(1:end-1),[]);
                 [I2,~]=findClosest(p.sample_window(1,2),range(1:end-1),[]);
                 [I3,~]=findClosest(p.sample_window(2,1),range(1:end-1),[]);
                 [I4,~]=findClosest(p.sample_window(2,2),range(1:end-1),[]);
                 sample_windowI=[I1 I2;I3 I4];
             else
                 sample_windowI=[1 length(range)-1;1 length(range)-1];
             end
%              newrange=range(sample_windowI(1):sample_windowI(2));
             
              if isfield(p,'alpha') && isnumeric(p.alpha)
                  alpha=p.alpha;
              else
                  alpha=0.05;
              end
             if isfield(p,'tail') && (isstring(p.tail) || ischar(p.tail))
                  tail=p.tail;
              else
                  tail='both';
             end
             
             if isfield(p,'method') 
                 method=p.method;
             else
                 method='ttest';%'Wilcoxon'; %
             end

             %data slicing

             %              ind = sub2ind(size(Rbool),i,j);
             mat_type=cell(1,length(type));
             for i=1:length(type)
                 Ri=Rslices(i);
                 mat=Ri.PSTHmat{1};
                 mat_type{i}=mat(:,sample_windowI(i,1):sample_windowI(i,2));
             end
             
             x=mean(mat_type{1},2);
             y=mean(mat_type{2},2);
             
             if strcmp(method,'ttest')
                 [H,P] = ttest2(x,y,'Alpha',alpha,'Tail',tail);
             elseif strcmp(method,'Wilcoxon')
                 [P,H] = ranksum(x,y,'Alpha',alpha,'Tail',tail);
             else
                 [H,P] = ttest2(x,y,'Alpha',alpha,'Tail',tail);
             end
         end
         function [TouchH,PistonH,anovaH,LightH]=stat_summary(obj,type)
             temp=obj.R(:,:,type);
             TouchH=zeros(size(temp));
             PistonH=zeros(size(temp));
             for i=1:size(temp,1)
                 for j=1:size(temp,2)
                     try
                        TouchH(i,j)=temp(i,j).stat.TouchSig;
                     catch
                        TouchH(i,j)=nan;
                     end
                     try
                         PistonH(i,j)=temp(i,j).stat.PistonSig;
                     catch
                         PistonH(i,j)=nan;
                     end
                 end
             end
             anovaH=(obj.Rstat.p<0.05);
             try
                LightH=[obj.Rstat_wilcoxon.h];
             catch
                 LightH=nan;
             end
%              touchsig(j)=(4+j,1).stat.TouchSig;
         end
         function [MUs,STDs,Ns]=PSTH_summary(obj,type)
             temp=obj.R(:,:,type);
             MUs=zeros(size(temp));
             STDs=zeros(size(temp));
             Ns=zeros(size(temp));
             for i=1:size(temp,1)
                 for j=1:size(temp,2)
                     MUs(i,j)=temp(i,j).MU;
                     STDs(i,j)=temp(i,j).SD;
                     Ns(i,j)=temp(i,j).N;
                 end
             end 
         end
         
         function [mu,sd]=norm_factor(obj,Rbool,type)%,Rbool,type)
             %get the mu and std for all trials (full window used) of given
             %Rbool slice
             %should use ctrl trials aligned to pistonstart or trialstart
             %to cover more data
             %check with res
             
%              Rbool=false(size(obj.R,2),size(obj.R,1));
%              Rbool(1,end)=1;

%              Rbool=true(size(obj.R,2),size(obj.R,1));
%              type=find(strcmp({obj.PARAMS.type},'PistonStart'));

             %end:free running;1:Ctrl;3:PistonStart;
             mat=[];
             for i=1:length(type)
                 [Rslice,Pslice,Rbool,Rcolor,fp_id]=obj.checkinput(Rbool,[],type(i));
                 Rs=Rslice(Rbool);
                 for j=1:length(Rs)
                     mat=[mat;Rs(j).PSTHmat{:}];
                 end
%                  t=t(1:length(Rs(1).PSTHs{1}));
             end
             dt=Pslice.PSTHrange(2)-Pslice.PSTHrange(1);

             [TrN,winN]=size(mat);
             
             x=sum(mat,2)/(winN*dt);
             
             mu=mean(x);
             sd=std(x);
             
%              mu=Rs.MU;
%              std=Rs.SD;
             
         end
         
         
         function [TRType,sig_arr]=touchResponseType(obj,optic)
%              type   BLF->TAF   Pre->post
%              1:     0           0
%              2:     0           1   (facilitated)
%              3:     0           -1  (suppressed)
%              4:     1           0
%              5:     1           1
%              6:     1           -1
%                   Rs().stat.PistonSig  Rs().stat.TouchSig
             StimN=size(obj.R,1);
             if ~exist('optic','var') || ~optic
                Rbool=logical([ones(1,StimN);zeros(1,StimN)]);
             else
                 Rbool=logical([zeros(1,StimN);ones(1,StimN)]);
             end
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,[],1);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             
             L=length(Rs);
             sig_arr=nan(2,L);
             for i=1:L
                 sig_arr(1,i)=Rs(i).stat.PistonSig;
                 sig_arr(2,i)=Rs(i).stat.TouchSig;
             end  
             sig_arr(isnan(sig_arr))=0;
             
             for i=1:L
                 if sig_arr(2,i)~=0  %not zero, determine faciliatated or suppreseed                 
                     temp1=abs(Rs(i).PREZERO-Rs(i).PEAKS)>abs(Rs(i).PREZERO-Rs(i).TROUGHS);
                     temp2=Rs(i).PREZERO<Rs(i).MU;
                     temp3=Rs(i).sign;                     
%                      if temp1~=temp2 && temp2~=temp3
%                          figure
%                          plot(Pslice.PSTHrange(1:end-1),Rs(i).PSTHs{1})
%                          123;
%                      end
                     if temp2
                         sig_arr(2,i)=1;
                     else
                         sig_arr(2,i)=-1;
                     end                     
                 end
             end
             
             TRType=nan(1,L);
             for i=1:L
                 switch sig_arr(1,i)
                     case 0
                         switch sig_arr(2,i)
                             case 1
                                 TRType(i)=2;
                             case 0
                                 TRType(i)=1;
                             case -1
                                 TRType(i)=3;
                         end   
                     case 1
                         switch sig_arr(2,i)
                             case 1
                                 TRType(i)=5;
                             case 0
                                 TRType(i)=4;
                             case -1
                                 TRType(i)=6;
                         end  
                 end
             end
             
             if sum(isnan(TRType))>0
                 123;
             end
             
%              obj.Rstat.c
             
         end
         function plotRaster(obj,Rbool,Rcolor,type,p)
             
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));
             
              if exist('p','var') && isfield(p,'sortI'), sortI=p.sortI;
             else sortI='none';end

             gcf,hold on
             data=[];corrTr=[];corrR=[];
             if strcmp(sortI,'chrono')    %all Rs in one single plot in chrono order
                 for i=1:length(Rs)
                     r=Rs(i);
                     data=[data;r.PSTHrelT{1}]; 
                     corrTr{i,1}=r.CorrTrial{1};
                     corrR{i,1}=ones(length(r.PSTHrelT{1}),1)*i;
                 end
                 corrTr=cell2mat(corrTr);
                 corrR=cell2mat(corrR);
                 [~,I]=sortrows(corrTr);
                 data=data(I);
                 corrR=corrR(I);
                 py=[1 length(data)];
                 rasterplot(data,Rcolor(corrR,:),py,[t(1) t(end)]);

             else %all Rs in one single plot stacked and sort by r.TouchDeltaT within each R 
                 count=0;
                 for i=1:length(Rs)
                     r=Rs(i);
                     data=r.PSTHrelT{1}; 
                     TouchN=length(data);
                     dt0=zeros(TouchN,1);
                     if isfield(r,'WID_FW') && ~isempty(r.WID_FW) && ~strcmp(sortI,'none')
                         [~,I]=sortrows([r.WID_FW r.TouchDeltaT]);
                         data=data(I);
                         dt=[dt0 r.TouchDeltaT(I,:)];
                     elseif ~strcmp(sortI,'none') && ~isempty(r.TouchDeltaT)
                         try
                         [~,I]=sortrows([r.TouchDeltaT]);
                         data=data(I);
                         dt=[dt0 r.TouchDeltaT(I,:)];
                         catch
                             dt=dt0;  %prob not using touch related alignment
                         end
                     else
                         dt=dt0;
                     end
                     py=[count+1 count+length(data)];
                     count=count+length(data);
                     rasterplot(data,Rcolor(i,:),py,[t(1) t(end)]);
                     hold on
                     try
                     for k=1:size(dt,2)
                        h=plot(dt(:,k),py(1):py(end));
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                     end
                     catch
                         disp('Cannot plot raster')
                     end
    %                  scatter(dt,,'MarkerFaceColor',c(j,:))
                 end
             end
            
             
             
             
         end
         function [y,t]=plotTrialTouchType(obj,Rbool,Rcolor,type,p)
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));
             
              if exist('p','var') && isfield(p,'countType'), countType=p.countType;
             else countType='TouchN';end
                if exist('p','var') && isfield(p,'isplot'), isplot=p.isplot;
             else isplot=1;end
             
             gcf,hold on
             data=[];corrTr=[];corrR=[];
             
                 for i=1:length(Rs)
                     r=Rs(i);
                     data=[data;r.PSTHrelT{1}]; 
                     corrTr{i,1}=r.CorrTrial{1};
                     corrR{i,1}=ones(length(r.PSTHrelT{1}),1)*i;
                 end
                 corrTr=cell2mat(corrTr);
                 corrR=cell2mat(corrR);
                 [~,I]=sortrows(corrTr);
                 
                 data=data(I);
                 corrR=corrR(I);
                 py=[1 length(data)];
%                  rasterplot(data,Rcolor(corrR,:),py,[t(1) t(end)]);
                 G=findgroups(corrTr(I));
                 for i=1:length(unique(G))
                     count(i)=sum(G==i);
                     type(i)=corrR(find(G==i,1,'first'));
                     spks=cell2mat(data(G==i));
                     Spikecounts(i)=length(spks);
                     SpikecountsPost(i)=sum(spks>=0);
                 end
                 t=1:length(unique(G));
                 switch countType
                     case 'TouchN'
                         y=count;
                     case 'SpikeN'
                         y=Spikecounts;
                     case 'SpikePerTouch'
                         y=Spikecounts./count;
                     case 'SpikePerPostTouch'
                         y=SpikecountsPost./count;
                 end
                 if isplot
                     patch([t nan],[y nan],[type nan],'EdgeColor','interp','Marker','o','MarkerFaceColor','flat');
                 end
             
         end
         
         function plotLinear(obj,Rbool,Rcolor,type) 
%              simply plot R(s).psth_linear with error R(s).confidenceI_linear if
%              found. Otherwise user should use linearR_bootstrapped first



             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));
             strs=Pslice.titleStr(Rbool);
             for i=1:length(Rs)
                 try
                     y=Rs.psth_linear;
                     CI_u=Rs.confidenceI_linear(1,:);
                     CI_l=Rs.confidenceI_linear(2,:);
                 catch
                     disp('Missing Rs.psth_linear or Rs.confidenceI_linear, run N.linearR_bootstrapped first')
                 end
                 C1=Rcolor(i,:);
                 C2=C1*0.8;
                ALPHA=0.3;
                h=plotWFilledError(t,y,C1,CI_l,CI_u,C2,ALPHA);
                h.DisplayName=strs;
                hold on
             end
         end
         
         function [psth_linear,confidenceI,R_temp,t] =linearR_deblf_bootstrapped(obj,Rbool,type,otherR,p)
             %Give singular value pretouch rate linear adjusted by blf before pistonIn
             %i.e. looks at overall response change from spontaneous
             if exist('p','var') && isfield(p,'isplot')
                 isplot=p.isplot;
             else
                 isplot=1;
             end
             if exist('p','var') && isfield(p,'color')
                 c=p.color;
             else
                 c=[0 0 0];
             end
             %with bootstrapping to get sample distribution
             [Rslice,Pslice,Rbool,~]=checkinput(obj,Rbool,[],type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));    
             %add other R slices if provided
             if exist('otherR','var')
                 if length(otherR)==1 && isfield(otherR,PSTHs) && isfield(otherR,PREZERO) && isfield(otherR,PSTHmat)
                    Rs=[Rs;otherR];
                 end
             end
             %calculate prezero avg
             prezeroMU=cell(length(Rs),1);
             blfiring=zeros(length(Rs),1);
             for i=1:length(Rs)
                prezeroMU{i}=sum(Rs(i).PSTHmat{1}(:,t<0),2)/(-1*t(1));
                blfiring(i)=Rs(i).BLF;
             end
             
             if length(unique(blfiring))==1
                 blf=blfiring(1);
             else
                 blf=mean(blfiring);
             end
             prezeroMU_deblf=cell(length(Rs),1);
             for i=1:length(Rs)
                prezeroMU_deblf{i}=prezeroMU{i}-blf;
             end
             
             % sum of all slices by bootstrap
             nRep=1000;
             psth_linear=zeros(1,1);
             confidenceI=zeros(2,1);
             obsMat=zeros(nRep,1);
             for i=1:length(psth_linear)
                 if length(Rs)==2
                     [H,sampStat,CI,ti_obs]=bootstrapping(prezeroMU_deblf{1}(:,i),prezeroMU_deblf{2}(:,i),@(x1,x2) mean(x1)+mean(x2),nRep);
                 elseif length(Rs)==3
                     [H,sampStat,CI,ti_obs]=bootstrapping3(prezeroMU_deblf{1}(:,i),prezeroMU_deblf{2}(:,i),prezeroMU_deblf{3}(:,i),@(x1,x2,x3) mean(x1)+mean(x2)+mean(x3),nRep);
%                      [bootstat,bootsam] = bootstrp(1000,@sum,Rmu_dePrezero{1}(:,i),Rmu_dePrezero{2}(:,i),Rmu_dePrezero{3}(:,i));
                 else
                     error('Need 2 or 3 to get sum')
                 end
                 psth_linear(i)=sampStat+blf;
                 confidenceI(1,i)=CI(1)+blf;
                 confidenceI(2,i)=CI(2)+blf;
                 obsMat(:,i)=ti_obs+blf;
             end
             %create fake struct R for the bootstrapped linearly combined
             %sample distribution
             R_temp.PSTHmat{1}=obsMat;
             R_temp.blf=blf;
             R_temp.PSTHs{1}=psth_linear; 
             if isplot
                plotWFilledError(0,psth_linear,c,confidenceI(1,:),confidenceI(2,:),c*0.8,0.3);
             end
         end

         function [psth_linear,confidenceI,R_temp,t] =linearR_bootstrapped(obj,Rbool,type,otherR,p)
             % Give posttouch linear adjusted by mean of (all single whisker
             % pretouch)
             %i.e. looks at only fine temporal response change from touch
             if exist('p','var') && isfield(p,'isplot')
                 isplot=p.isplot;
             else
                 isplot=1;
             end
             if exist('p','var') && isfield(p,'color')
                 c=p.color;
             else
                 c=[0 0 0];
             end
             if exist('p','var') && isfield(p,'blfType')
                 blfType=p.blfType;
             else
                 blfType='prepiston';  %pretouch
             end
             
             %with bootstrapping to get sample distribution
%              [Rslice,Pslice,Rbool,Rcolor,fp_id,Rbool_flipped]
             [Rslice,Pslice,Rbool,~,~,~]=checkinput(obj,Rbool,[],type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));       
             
             %add other R slices if provided
             if exist('otherR','var')
                 if length(otherR)==1 && isfield(otherR,PSTHs) && isfield(otherR,PREZERO) && isfield(otherR,PSTHmat)
                    Rs=[Rs;otherR];
                 end
             end
            
             %calculate prezero avg
             prezeroMU=zeros(length(Rs),1);
             if strcmp(blfType,'pretouch')       
                 for i=1:length(Rs)
                     prezeroMU(i)=Rs(i).PREZERO;
                 end
                 prezeromean=mean(prezeroMU);
             elseif strcmp(blfType,'prepiston')
                 for i=1:length(Rs)
                     prezeroMU(i)=Rs(i).BLF;
                 end
                 prezeromean=mean(prezeroMU);
             elseif strcmp(blfType,'runningCtrl')
                 for i=1:length(Rs)
                     prezeroMU(i)=Rs(i).BLF;
                 end
                 prezeromean=mean(prezeroMU);
             else
                 error('unknown blfType')
             end
%              prezero_binT=-t(1);
             
             Rmu_dePrezero=cell(length(Rs),1);
             for i=1:length(Rs)
                 Rmu_dePrezero{i}=Rs(i).PSTHmat{1}/Pslice.bin-prezeroMU(i);
             end
             % sum of all slices by bootstrap
             nRep=1000;
             psth_linear=zeros(1,length(Rs(1).PSTHs{1}));
             confidenceI=zeros(2,length(Rs(1).PSTHs{1}));
             obsMat=zeros(nRep,length(Rs(1).PSTHs{1}));
             for i=1:length(psth_linear)
                 if length(Rs)==2
                     [H,sampStat,CI,ti_obs]=bootstrapping(Rmu_dePrezero{1}(:,i),Rmu_dePrezero{2}(:,i),@(x1,x2) mean(x1)+mean(x2),nRep);
                 elseif length(Rs)==3
                     [H,sampStat,CI,ti_obs]=bootstrapping3(Rmu_dePrezero{1}(:,i),Rmu_dePrezero{2}(:,i),Rmu_dePrezero{3}(:,i),@(x1,x2,x3) mean(x1)+mean(x2)+mean(x3),nRep);
%                      [bootstat,bootsam] = bootstrp(1000,@sum,Rmu_dePrezero{1}(:,i),Rmu_dePrezero{2}(:,i),Rmu_dePrezero{3}(:,i));
                 else
                     error('Need 2 or 3 to get sum')
                 end
                 psth_linear(i)=sampStat+prezeromean;
                 confidenceI(1,i)=CI(1)+prezeromean;
                 confidenceI(2,i)=CI(2)+prezeromean;
                 obsMat(:,i)=ti_obs+prezeromean;
             end
             %create fake struct R for the bootstrapped linearly combined
             %sample distribution
             R_temp.PSTHmat{1}=obsMat;
             R_temp.PREZERO=prezeromean;
             R_temp.PSTHs{1}=psth_linear; 
             if isplot
                h=plotWFilledError(t,psth_linear,c,confidenceI(1,:),confidenceI(2,:),c*0.8,0.3);
                h.DisplayName='Lin';
             end
         end
         function [Y,Er,t]=plotPSTH(obj,Rbool,Rcolor,type,p)
             % Rbool determines which stim to plot
             % Rbool=[1 0 0....;1 0 0....] Nrow=2 (usually second row is light); Ncolumn=Sstim+MStim+Free (~10)
             % Rcolor=[0 0 0;0 0 1;....]   Nrow=N of Stim to plot; Ncolumn=RGB
             % p params:
             % plotType='linear'(default),'heatmap'
             % isplot=1 (default)
             % norm='none'(default),'maxabs','spontaneous'
             % plotError=1(default),0
             % interpolate=0(default),1 **only apply to heatmap for now
             
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             if isempty(Rs)
                 Y=[];Er=[];t=[];
                 return
             end
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));
             strs=Pslice.titleStr(Rbool);
             
             try
                [Rslice_ret,~,~]=checkinput(obj,Rbool,Rcolor,'Retraction_Tonly');   
                Rslice_ret=Rslice_ret(Rbool);
                Rslice_ret=Rslice_ret(:);
             catch
             end

             % check p
             if exist('p','var') && isfield(p,'plotType'), plotType=p.plotType;
             else plotType='linear';end
                
             if exist('p','var') && isfield(p,'isplot'), isplot=p.isplot;
             else isplot=1;end
             
             if exist('p','var') && isfield(p,'norm'), Norm=p.norm;
             else Norm='none';end
             switch Norm
                 case 'maxabs'
                     h = @(x) x/max(abs(x));
                     he = @(x1,x) x1/max(abs(x));
                 case 'spontaneous'
                     [mu_norm,std_norm]=obj.norm_factor;
                     h = @(x) (x-mu_norm)./std_norm;
                     he = @(x1,x) x1./std_norm;
                     
                 otherwise
                     h = @(x) x;
                     he = @(x1,x) x1;
             end
             if exist('p','var') && isfield(p,'plotError'), plotError=p.plotError;
             else plotError=1;end
             
             if exist('p','var') && isfield(p,'interpolate'),interpolate=p.interpolate;
             else interpolate=0;end
             if interpolate
                 MultipleFactor=3;
                h_int = @(x) smoothdata(resample(x,MultipleFactor,1),'gaussian',3); 
%                 C = smoothdata(A,'gaussian',20);
             else
                 MultipleFactor=1;
                 h_int = @(x) x;
             end
             
                        
             if strcmp(plotType,'heatmap')
                 
                 Hmap=nan(length(Rs),length(t));
                 for i=1:length(Rs)
                      if Rs(i).N==0 && (strcmp('Touch',type) || type==1)  %no touch, i.e. free whisking stim
                         if exist('Rslice_ret','var')
                             Hmap(i,:)=h(Rslice_ret(i).PSTHs{1});
                         else
                             Hmap(i,:)=h(Rs(i).TAF*ones(length(t)));
                         end
                     else
                         Hmap(i,:)=h(Rs(i).PSTHs{1});                         
                     end
                 end
                 Hmap=h_int(Hmap')';
                 t_backup=t;
                 t=h_int(t);
                 ystr=Pslice.titleStr(Rbool);
                 for i=1:length(Rs)
                    if ~isnan(Rs(i).stat.TouchSig) && Rs(i).stat.TouchSig
                        ystr{i}=['*' ystr{i}];
                    end
                 end
                 if isplot
                     hold off
                    heatmap(t,ystr,Hmap,'Colormap', hot);
                    xlabel('Time from onset (s)')
                 end      
                 xstr = num2cell(t);
                 count=0;
                 for i=1:length(xstr)
                     if rem(i,MultipleFactor)==1
                         count=count+1;
                         xstr{i}=t_backup(count);
                     else
                         xstr{i}=' ';
                     end
                 end
                 if isplot
                    ax=gca;
                    ax.XDisplayLabels=xstr;
                 end
                 Y=Hmap;
                 Er=[];

             else  %linear(default)
                 Y=zeros(length(Rs),length(t));
                 Er=zeros(length(Rs),length(t));
                 for i=1:length(Rs)
                     str=sprintf('%s n:%d',strs(i),Rs(i).N);
                     if Rs(i).N==0 && (strcmp('Touch',type) || type==1)  %no touch, i.e. free whisking stim
                         if exist('Rslice_ret','var')
                             %                          plot(t,Rslice_ret(i).PSTHs{1},'Color',Rcolor(i,:));
                             x=Rslice_ret(i).PSTHs{1};
                             x_sem=Rslice_ret(i).SEMs{1};
                             if isplot
                                 gcf;hold on  
                                if plotError
                                 plotWE(t,h(x),he(x_sem,x),Rcolor(i,:),str)
                                else
                                 plot(t,h(x),'Color',Rcolor(i,:),'DisplayName',str);
                                end
                             end
                         else
                             x=Rs(i).TAF*ones(length(t));
                             if isplot
                                 gcf;hold on  
                                plot(t,h(x),'Color',Rcolor(i,:),'DisplayName',str);
                             end
                         end
                     else
                         %                       plot(t,Rs(i).PSTHs{1},'Color',Rcolor(i,:));
                         x=Rs(i).PSTHs{1};
                         x_sem=Rs(i).SEMs{1};
                         if isplot
                             gcf;hold on  
                             if plotError
                                 plotWE(t,h(x),he(x_sem,x),Rcolor(i,:),str);
                             else
                                 plot(t,h(x),'Color',Rcolor(i,:),'DisplayName',str);
                             end
                         end
                     end
                     Y(i,:)=h(x);  
                     Er(i,:)=he(x_sem,x);  
                 end   
                 if isplot
                    xlabel('Time from onset (s)')
                    ylabel('Spike rate s-1')
                    H=legend;
                    set(H,'Location','northwest','NumColumns',1,'FontSize',6);
                    H.ItemTokenSize = [5,10]; %reduce line size
                 end
             end
   
         end
         
         function plotTriangleDiagram(obj,Exp,Rbool,type,p)
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,[],type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             if isempty(Rs)
%                  Y=[];t=[];
                 return
             end
             t=Pslice.PSTHrange;
             t=t(1:length(Rs(1).PSTHs{1}));
             MUs=[Rs.MU];   % MUs=cellfun(@(x) x.MU,Rs);
             blf=mean([Rs.BLF]);  %should be the same for all Rs anyway
             
             
             z=MUs;  %z determines the color value
             z_zero=0;  %value of z that represents zero on colorbar  %blf
             scalingFactor=1/max(abs(z-z_zero));
%              scalingFactor=1/std(z);
             z_scaled=(z-z_zero)*scalingFactor;
             
             C=createColorMap(501,[1 0 0;1 1 0;0 1 0;0 1 1;0 0 1]);
             cvalue=[-250:250]/250; %[0:500]/500;
             
             [Iidx,~]=findClosest(z_scaled,cvalue,[]);
             z_c=C(Iidx,:);
             
             blf_scaled=(blf-z_zero)*scalingFactor;
             [Iidx,~]=findClosest(blf_scaled,cvalue,[]);
             blf_c=C(Iidx,:);
             
             recipe=Exp.MStim.recipe;
             
             C1_vertex=[0 0];
             B1_vertex=[-1/3 (1/3)^0.5]*3/2;
             D1_vertex=[-1/3 -(1/3)^0.5]*3/2;
             C2_vertex=[2/3 0]*3/2;
             
             vertices=[C1_vertex;
                 B1_vertex;
                 D1_vertex;
                 C2_vertex];
             
             W_order=[2 1 3 4]'; %C1 B1 D1 C2
             
             vertices_sorted=vertices(W_order,:);
            
             %scatter singles
             hold on
             for i=5:length(Rs)
                Ws_used=find(recipe(:,i));
                c=z_c(i,:);
                line(vertices_sorted(Ws_used,1),vertices_sorted(Ws_used,2),'Color',c,'LineWidth',3);
             end
              %plot lines for pairs
             for i=1:4
                Ws_used=i;
                c=z_c(i,:);
                scatter(vertices_sorted(Ws_used,1),vertices_sorted(Ws_used,2),100,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',c);
             end
             axis equal
             colormap(C)
             cbh = colorbar ; %Create Colorbar
             cbh.Ticks = linspace(0, 1, 3) ; %Create n ticks from zero to 1
             cbh.TickLabels = num2cell([cvalue(1) cvalue(ceil(length(cvalue)/2)) cvalue(end)]) ;
%              set(gca,'Color',blf_c)
         end
         
         function plotFirstSpike(obj,Exp,Rbool,Rcolor,type)
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);

             temp=squeeze([Rs.firstSpkLats]);
             for j=1:length(Rs)
                 [h,edges]=histcounts(temp{j},'BinWidth',0.002,'Normalization','probability'); t=edges(1:end-1)+(edges(2)/2-edges(1)/2);
                 try
                    plot(t,h,'LineStyle','-','Color',Rcolor(j,:));    hold on,
                 catch
                 end
             end
             xlabel('First Spike from touch (s)')
             ylabel('Spk Prob.')
             legend(Pslice.titleStr{Rbool})
             title(sprintf('IDCID:%d,%d %s pw:%d shankpw:%d depth:%d\n',obj.ID,obj.CID,obj.Shape,obj.pw,Exp.KSort.shankPW(obj.Shank),obj.Depth))
         end
         
         function plotSpkTendencyJitter(obj,Rbool,Rcolor,type)
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);

             temp1=squeeze([Rs.firstSpkTendency]);
             temp2=squeeze([Rs.firstSpkJitter]);  
             
%              temp=squeeze([Rs.firstSpkLats]);temp2=nan(1,length(Rs));  %jitter but norm by number of touches
%              for j=1:length(Rs)
%                 [h,edges]=histcounts(temp{j},'BinWidth',0.002,'Normalization','probability');  
%                 temp2(j)=std(h,'omitnan');
%              end
             gcf;hold on;
             for i=1:length(Rs)
                scatter(temp1(i),temp2(i),[],Rcolor(i,:),'filled')
             end
             xlabel('Spk Tendency')
             ylabel('Spk Jitter')
             legend(Pslice.titleStr{Rbool})
         end
         function plotSpkCount(obj,Rbool,Rcolor,type)
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,type);   
             Rs=Rslice(Rbool);
             Rs=Rs(:);
             temp=squeeze([Rs.countpertrial]);
             for j=1:length(Rs)
                 [h,edges]=histcounts(temp{j}); t=edges(1:end-1)+(edges(2)/2-edges(1)/2);
                 try
                 plot(t,h/Rs(j).N,'LineStyle','-','Color',Rcolor(j,:));    hold on,
                 catch
                 end
             end
             xlabel('Spike Counts')
             ylabel('Touch Prob.')
             legend(Pslice.titleStr{Rbool})
         end
         function [VectorL,SpkProb,binCentre]=plotTuning(obj,Exp,Whisker,StimCs,Rbool,Rcolor,typeX,typeY,p)
             
             
             [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,1);  
             titleStr=Pslice.titleStr(Rbool);
             
             if size(StimCs)==size(Rslice)
                 StimCi=StimCs(Rbool);
                 %                 titleStr=titleStr(Rbool);
             elseif size(StimCs)==flip(size(Rslice))
                 StimCs=StimCs';
                 StimCi=StimCs(Rbool);
                 %                 titleStr=titleStr(Rbool');
             else
                 error('Check StimC input dimension')
             end
%              titleStr=reshape(Pslice.titleStr,flip(size(StimCs)))';  %quick fix only
             
             if ~exist('p','var')
                 p=[];
             end
             if ~isfield(p,'WhiskerRef'), p.WhiskerRef='Stim';end
             if ~isfield(p,'samples_from_touchingCyclesOnly'), p.samples_from_touchingCyclesOnly=1;end
             if ~isfield(p,'MaxSampleN'), p.MaxSampleN=[];end%2500;
             if ~isfield(p,'SampleMethod'), p.SampleMethod='all';end%2500;      
             if ~isfield(p,'Discretize_Edges'), p.Discretize_Edges=[];end%-pi:pi*2/20:pi;
             if ~isfield(p,'Discretize_binN'), p.Discretize_binN=[];end
             if ~isfield(p,'Color'), p.Color=Rcolor;end
             if ~isfield(p,'plotFig'), p.plotFig=1;end   %plot tuning curve    
             if ~isfield(p,'plottype'), p.plottype='Pol';end %'Cart' 'Pol'
             if ~isfield(p,'useRate'), p.useRate=0;end %'Cart' 'Pol'
             if ~isfield(p,'plotvectorL'), p.plotvectorL=0;end   %plot vector length
             
             if p.plotFig && p.plotvectorL
                 error('Pick one plot only')
             end
             %p.titleStr

             
             StimCii=StimCi;
             for i=1:length(StimCi)
                 Si=StimCi{i}.PWadjusted(Exp,Whisker,obj);
                 if ~isempty(Si.WID)
                     Si=updateStimC(Si,Exp,Pslice(1));
                 end
                 StimCii{i}=Si;
             end
             [SpkProb,Edges,N]=tuningcurve(obj,Exp,Whisker,StimCii,typeX,typeY,p);
             d=Edges(2)-Edges(1);
             binCentre=Edges(1:end-1)+d/2;
             titleStr=cellfun(@(x,y) sprintf('%s N=%d',x,y),toColumn(titleStr),toColumn(num2cell(N)),'UniformOutput',false);
             
%              if p.plotFig
%                  hold off  %avoid cart/pol icompatibilty with current figure
%              end
            if p.plotFig
             for i=1:length(StimCi)
                 if strcmp(p.plottype,'Cart') 
                     plot(binCentre,SpkProb(i,:),'Color',Rcolor(i,:),'LineWidth',1);
                     hold on
                 else%if p.plotFig
                     polarplot([binCentre binCentre(1)],[SpkProb(i,:) SpkProb(i,1)],'Color',Rcolor(i,:),'LineWidth',1);
                     hold on
                 end
                 
                 %[x,y] = pol2cart(angle(L),abs(L));
                 %                 quiver(0,0,x,y);
             end
             %              cellfun(@plot,titleStr,Y)
             legend(titleStr)
             title([typeX ' tuned ' typeY])
            end
             
             VectorL=zeros(length(StimCi),1);
             for i=1:length(StimCi)
                 VectorL(i)=vector_length(binCentre,SpkProb(i,:));
             end
             
             if p.plotvectorL
                 for i=1:length(StimCi)    
                     L=VectorL(i);
                     polarplot([angle(L) angle(L)],[0 abs(L)],'Color',Rcolor(i,:),'LineWidth',3)  %,'HandleVisibility','off'
                     hold on
                 end
                 legend(titleStr)
             end
         end
         
         function plotSTA_Wtouch(obj,Exp,Whisker,StimCs,Rbool,Rcolor,STArange,touchType)
             
             if ~exist('touchType','var') || isempty(touchType)
                touchType=1:4;
                multiplePlots=1;
             elseif length(touchType)==1
                multiplePlots=0;
             else
                 multiplePlots=1;
             end
             [ai,aj]=subplotDim(length(touchType));
%              figure
             for m=1:length(touchType)
                 if multiplePlots
                    subplot(ai,aj,touchType(m))
                 end
                 switch touchType(m)
                     case 1
                         touchbinaryType='onset_only';
                         smooth_it=1;
                     case 2
                         touchbinaryType='offset_only';
                         smooth_it=1;
                     case 3
                         touchbinaryType='continuous';
                         smooth_it=0;
                     case 4
                         touchbinaryType='onoffset_only';
                         smooth_it=1;
                     otherwise
                         error('touchType should range from 1 to 4')
                 end
                 
                 
                 
                 spkType='all';
                 %              spkType='StartofTrain_only';
%                  StimLineStyle={'-','--',':','-.'};
                 WhiskerColor=[1 0 0;0 1 0;0 0 1;1 0 1];
                 try
                    [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,'Touch');
                 catch
                     [Rslice,Pslice,Rbool,Rcolor]=checkinput(obj,Rbool,Rcolor,1);
                 end
                 titleStr=Pslice.titleStr(Rbool);
                 if size(StimCs)==size(Rslice)
                     StimCi=StimCs(Rbool);
                     %                 titleStr=titleStr(Rbool);
                 elseif size(StimCs)==flip(size(Rslice))
                     StimCs=StimCs';
                     StimCi=StimCs(Rbool);
                     %                 titleStr=titleStr(Rbool');
                 else
                     error('Check StimC input dimension')
                 end
                 %Get binary touch for each WT
                 WhkBhvr=cell(length(Whisker),Exp.TrN);
                 for i=1:length(Whisker)
                     W=Whisker(i);
                     for j=1:Exp.TrN
                         WT=W.trial(j);
                         switch touchbinaryType
                             case 'onset_only'
                                 WhkBhvr{i,j}=createBinary(WT.FrameN,WT.TouchFrame,1);
                             case 'offset_only'
                                 WhkBhvr{i,j}=createBinary(WT.FrameN,WT.ReleaseFrame,1);
                             case 'continuous'
                                 WhkBhvr{i,j}=createBinary(WT.FrameN,WT.TouchFrame,WT.ReleaseFrame);
                             case 'onoffset_only'
                                 WhkBhvr{i,j}=createBinary(WT.FrameN,sort([WT.TouchFrame;WT.ReleaseFrame]),1);
                         end
                         
                     end
                 end
                 
                 for i=1:length(StimCi)
                     
                     S=StimCi{i};
                     STA_Bhvr=cell(1,length(Whisker));
                     for j=1:length(Whisker)
                         if ~ismember(j,S.WID)% && ~isempty(S.WID)
                             continue
                         end
                         [STA_Bhvr{j},t]=obj.STA(Exp,WhkBhvr(j,:),S,STArange);
                     end
                     
                     for j=1:length(Whisker)
                         BrightnessTerm=0.7^(j-1);
                         if ~isempty(STA_Bhvr{j})
                             if smooth_it
                                plot(t,smooth(STA_Bhvr{j}),'Color',min([Rcolor(i,:)*BrightnessTerm;1 1 1],[],1),'DisplayName',sprintf('%s W%d',S.type,j));
                             else
                                 plot(t,STA_Bhvr{j},'Color',min([Rcolor(i,:)*BrightnessTerm;1 1 1],[],1),'DisplayName',sprintf('%s W%d',S.type,j));                                
                             end
                             hold on
                         end
                     end
                     xlabel('Time from Spike(s)')
                     ylabel('Touch prob.')
                     title(touchbinaryType)
                     legend
                 end
             end
         end
         
         function [Rslice,Pslice,Rbool,Rcolor,fp_id,Rbool_flipped]=checkinput(obj,Rbool,Rcolor,type)
             if  isempty(obj.R) || ~isstruct(obj.R)
                error('Obtain and save R from RES_script first')
             end
             if ~islogical(Rbool)
                 error('Rbool should be logical');
             end
             
             if isstring(type) || ischar(type)
                fp_id=find(strcmp(strtrim(string(char(obj.PARAMS.type))),type));
             elseif isnumeric(type) && length(type)==1
                 fp_id=type;
             else
                 fp_id=[];
             end
             
             if isempty(fp_id)
                 error('No %s Res found. Run res_script first',type)
             elseif length(fp_id)>1
                 error('One type per time plz')
             end
             
             Rslice=obj.R(:,:,fp_id);
             Pslice=obj.PARAMS(fp_id);
             if size(Rbool)~=size(Rslice)
                Rbool=Rbool';
                if size(Rbool)~=size(Rslice)
                    error('Rbool dimension should match StimC * w/wo light (2)')
                end
                Rbool_flipped=1;
             else
                Rbool_flipped=0;
             end
               
             if sum(Rbool,'all')>size(Rcolor,1)            
%                  if ~isempty(Rcolor)
%                      fprintf('Not enough colors for %d plots, using defaults\n',sum(Rbool,'all'))
%                  end
                 Rcolor=createColorMap(sum(Rbool,'all'),[0 0 0;0 0 1;0 1 0;1 0 0]);
             else
                 Rcolor=Rcolor(1:sum(Rbool,'all'),:);
             end
         end
         
         function mwbest=bestMW(obj)
            %output significant control MW id sorted by descending MU response
             
              mu=obj.touchMU(1,5:9);
              I_sig=logical(obj.PW_H_Touch(1,5:9));
              MWid=5:9;
              
             [mu_best,I]=sort(mu(I_sig),'descend');
             mwbest=MWid(I_sig);
             I_sort=I(~isnan(mu_best));
             mwbest=mwbest(I_sort);
             
         end
         

         function emptyTouch=raster(obj,timepoints,PSTHrange,color,py)
            px=zeros(1,2);
            if isempty(py)
                py=[]; %keep it empty
            end
            spkT=toColumn(obj.SpikeTime);
            px(1)=PSTHrange(1);px(2)=PSTHrange(end);
            data=cell(1,length(timepoints));
            for i=1:length(timepoints)
                data{i}=spkT(timepoints(i)+px(1)<spkT & timepoints(i)+px(2)>spkT)-timepoints(i);
            end
            emptyTouch=rasterplot(data,color,py,[PSTHrange(1) PSTHrange(end)]);
         end
         function [PSTHavg, SEM,SS,HaveSpikes]=PSTH_F(obj,timepoints,PSTHrange,CorrTrial)
%              Same as PSTH but only consider touch with at least 1 spike when plotting PSTH
%              Also gives filtered bool indexes
%              And removed trials histogram
                spkT=toColumn(obj.SpikeTime);
                px=[PSTHrange(1) PSTHrange(end)];
                HaveSpikes=true(1,length(timepoints));
                for i=1:length(timepoints)
                    data=spkT(timepoints(i)+px(1)<spkT & timepoints(i)+px(2)>spkT)-timepoints(i);
                    if  isempty(data)
                        HaveSpikes(i)=0;
                    end
                end
               if ~isempty(CorrTrial)
                   prevFig = gcf;
                   figure
                   histogram(CorrTrial(HaveSpikes))
                   set(0, 'currentfigure', prevFig);
               end
               [PSTHavg, SEM,SS]=PSTH(obj,timepoints(HaveSpikes),PSTHrange);
             
             
         end
         function [blf,sd,n]=baselineFiring(obj,Exp)
             bin = 0.01;
             startTime=0;
             winSize=Exp.Stim.Piston.Delay;
             I=~Exp.Stim.Optic.On;
             PSTHrange = (startTime:bin:startTime+winSize);
             [~, ~,SS]=PSTH(obj,Exp.TrialStartT(I),PSTHrange,[]);
             ratepertrial=mean(SS.PSTHmat,2)/bin;
             sd=std(ratepertrial);
             blf=mean(ratepertrial);
             n=length(ratepertrial);
         end
         function [blf,sd,n,ratepertrial]=baselineFiringLightOn(obj,Exp)
             bin = 0.01;
             startTime=0;
             winSize=Exp.Stim.Piston.Delay;
             I=Exp.Stim.Optic.On;
             PSTHrange = (startTime:bin:startTime+winSize);
             [~, ~,SS]=PSTH(obj,Exp.TrialStartT(I),PSTHrange,[]);
             ratepertrial=mean(SS.PSTHmat,2)/bin;
             sd=std(ratepertrial);
             blf=mean(ratepertrial);
             n=length(ratepertrial);
         end
         function obj=baselineDrift(obj,Exp,range)
             bin = 0.01;
             startTime=range(1);
             endTime=range(2);
%              startTime=-0.5;
%              endTime=0.4;%Exp.Stim.Piston.Delay;
             PSTHrange = (startTime:bin:endTime);
             [~, ~,SS]=PSTH(obj,Exp.TrialStartT,PSTHrange,[]);
             if Exp.TrN~=size(SS.PSTHmat,1)
                 error('Checkpoint, TrialNumber unequal')
             end
             Drift=zeros(1,Exp.TrN);
             for i=1:size(SS.PSTHmat,1)
                Drift(i)=mean(SS.PSTHmat(i,:))/bin;
             end
             
             plot(Drift)
             xlabel('TrialN')
             ylabel('Blf')
             
             c=-1/(sqrt(2)*erfcinv(3/2));
             th=c*median(abs(Drift-median(Drift)));         %th is too high, change to sth else when free
             hold on
             addline('y',th,'r','-');
%              obj.ActiveTrials=Drift>th;
%              obj.ActiveTh=th;

%              h=histogram(Drift);
%              y=h.Values;
%              x=h.BinEdges(1:end-1)+h.BinWidth/2;
%              f = fit(x.',y.','gauss1');
%              mu=f.b1;
%              std=f.c1/(2)^0.5;
             
         end
                      
             
         
         function taf=TrialAvgFiring(obj,Exp,StimC)
             %Trial structure may be different for different experiment so
             %check the window used. can use TrialTouchRaster.m to
             %vivsually check
             bin = 0.01;    
             if length(Exp.videoFps)==1
                startTime=min(min(Exp.PistonBuffer))/Exp.videoFps;   
             else
                startTime=Exp.Stim.Piston.Delay+0.1;
             end
%              startTime=min(min(Exp.PistonBuffer))/500;
             winSize=0.5;
             PSTHrange = (startTime:bin:startTime+winSize);
             [PSTHavg, ~,~]=PSTH(obj,Exp.TrialStartT(StimC.validTrials),PSTHrange,[]);
             taf=mean(PSTHavg);
         end
         function [Y,Edges,N]=tuningcurve(obj,Exp,Whisker,StimCs,Feat2TuneX,Feat2TuneY,p)
             
             %p
%              p.WhiskerRef
%              p.samples_from_touchingCyclesOnly
%              p.MaxSampleN
%              p.SampleMethod    "wReplacement" otherwise w/o
%              p.Discretize_Edges
%              p.Discretize_binN
%              p.Color
%              p.titleStr
%              p.useRate  %default=0:use probability

%           Input 
%              Feat2TuneX:continuous whisker feature found in Whisker.Trial, e.g. angle,
%           phase...

% Feat2TuneY: probability of dependent binary variable, 
% e.g. spike (otherwise default)
% Touch_Onset
% Touch_Offset
% Touch_Continuous



             %params
             if isfield(p,'WhiskerRef')
                 WhiskerRef=p.WhiskerRef;
             else
                WhiskerRef='Stim';
             end
             if isfield(p,'samples_from_touchingCyclesOnly')
                %==1, only touching cycles;==0, only non touching cycles;==[], all 
                samples_from_touchingCyclesOnly=p.samples_from_touchingCyclesOnly;              
                %note: if ==1 & WhiskerRef is 'PW', may filter out all samples if StimC does not have said whisker piston
             else
                 samples_from_touchingCyclesOnly=[];
             end
             if isfield(p,'MaxSampleN')
                MaxSampleN=p.MaxSampleN;
             else
                MaxSampleN=[];
             end
             if isfield(p,'SampleMethod')
                SampleMethod=p.SampleMethod;
             else
                SampleMethod=[];
             end
             if isfield(p,'Discretize_Edges')   %edges precedes binN
                Discretize_Edges=p.Discretize_Edges;
             else
                 Discretize_Edges=[];
             end
             if isfield(p,'Discretize_binN') && ~isempty(p.Discretize_binN)
                Discretize_binN=p.Discretize_binN;
             else
                 Discretize_binN=15;
             end
             if isfield(p,'Color')
                 c=p.Color;
             else
                 c=[];
             end
             if isfield(p,'titleStr')
                 titleStrs=p.titleStr;
             else
                 titleStrs=[];
             end
             if isfield(p,'useRate')
                 useRate=p.useRate;
             else
                 useRate=0;
             end
%              if isfield(p,'plotFig')
%                  isplot=p.plotFig;
%              else
%                  isplot=0;
%              end
             isplot=0;
               
             
             %get FeatX
             if ~iscell(StimCs)
                 error('Restructure StimCs into a cell array')
             end
             %get binned features
             if ~isempty(Discretize_Edges)
                 Y=zeros(length(StimCs),length(Discretize_Edges)-1);
             else
                 Y=zeros(length(StimCs),Discretize_binN);
             end
             
             N=zeros(1,length(StimCs));
             
             for i=1:length(StimCs)
                 S=StimCs{i};
                 TOI=S.validTrials;
                 
                 WhkBInT=cell(length(TOI),1);
                 WhkTr=cell(length(TOI),1);
                 WhkF=cell(length(TOI),1);
                 
                 switch WhiskerRef
                     case 'Stim'
                         if ~isempty(S.WID)
                            WOI=S.WID;
                            WOI=WOI(1);
                         else
                             WOI=2;
                         end
                     case 'PW'          
                         WOI=obj.pw;
                     case 'C1'
                         WOI=2;
                     case 'Avg'
                         WOI=1:4;      
                 end
                 if ~isempty(titleStrs)
                    titleStr=titleStrs(i);
                 else
                     titleStr=S.type;
                 end
                 
                 for j=1:length(TOI)
                     tr=TOI(j);
                     Feat=cell(1,length(WOI));
                     for k=1:length(WOI)
                        WT=Whisker(WOI(k)).trial(tr);   
                        Feat{k}=toColumn(WT.(Feat2TuneX));
                        if length(Feat{k})~=WT.FrameN
                            error('Invalid feature %s, dimension should match FrameN',Feat2TuneX)
                        end
                     end                     
                     FEAT=mean(cell2mat(Feat),2);  %average feats from all WOI                    
                     WhkBInT{j}=FEAT;     
                     WhkTr{j}=ones(Exp.FrameN(tr),1)*tr;
                     WhkF{j}=toColumn(1:Exp.FrameN(tr));
                 end
                 
                 %Get Filter bools
                 WhkFeatCat=cell2mat(WhkBInT);
                 WhkTrCat=cell2mat(WhkTr);
                 WhkFCat=cell2mat(WhkF);
                 
                 filter_bool=true(size(WhkFeatCat));
                 
                 for j=1:length(TOI)
                     tr=TOI(j);
                     trSampleBool=WhkTrCat==tr;
                     for k=1  %using first element in WOI
                        WT=Whisker(WOI(k)).trial(tr); 
                     end
                     if ~isempty(samples_from_touchingCyclesOnly)
                         PstartF=WT.P_Frame(1:end-1);  
                         PendF=WT.P_Frame(2:end);
                         if samples_from_touchingCyclesOnly && ~isempty(S.WID)
                             gt=find(WT.GoodTouch(1:end-1));
                             goodtouchbool=cell2mat(S.isGoodPairs(tr));  %good touch within trial tr after filtering with StimC
                             touching_bool=createBinary(WT.FrameN,PstartF(gt(goodtouchbool)),PendF(gt(goodtouchbool)));
                         else
                              %no touching, StimC independent
                             touching_bool=createBinary(WT.FrameN,PstartF(~WT.GoodTouch(1:end-1)),PendF(~WT.GoodTouch(1:end-1)));
                         end
                         filter_bool(trSampleBool)= (filter_bool(trSampleBool) & touching_bool);
                     end
                     
                     %vvvvv  add additional filters here vvvvv
                     
                     %------------------
                 end
                 
                 % Obtain SampleN considering MaxSampleN
                 if ~isempty(MaxSampleN)
                    if sum(filter_bool)>=MaxSampleN
                        SampleN=MaxSampleN;
                    else
                        SampleN=sum(filter_bool);
                    end
                 else
                     SampleN=sum(filter_bool);
                 end
                 % Sample by SampleMethod
                 switch SampleMethod
                     case 'wReplacement'
                         filter_id = datasample(find(filter_bool),SampleN);      
                     otherwise  %sample all without replacement
                         filter_id = datasample(find(filter_bool),SampleN,'Replace',false);
                 end
                         
                 N(i)=SampleN;
      
                 %filter and discretise FEATURE                             
                 if ~isempty(Discretize_Edges)
                    [Groups,Edges] = discretize(WhkFeatCat(filter_id),Discretize_Edges);
                 else
                     [Groups,Edges] = discretize(WhkFeatCat(filter_id),Discretize_binN);
                 end
                 
                 Edge_D=Edges(2)-Edges(1);
                 BinCentres=Edges(1:end-1)+Edge_D/2;
                 WhkFeatCat_discretized=BinCentres(Groups);  %..........Whisker Features in binned samples
                 
                 %FeatY
                 FeatY=cell(Exp.TrN,1);
                 switch Feat2TuneY
                     case 'Touch_Onset'
                         for k=1:Exp.TrN
                             for m=1:length(WOI)
                                WT=Whisker(WOI(m)).trial(k);
                                FeatY{k}=createBinary(WT.FrameN,WT.TouchFrame,1);
                             end
                         end
                     case 'Touch_Offset'
                         for k=1:Exp.TrN
                             for m=1:length(WOI)
                                WT=Whisker(WOI(m)).trial(k);
                                FeatY{k}=createBinary(WT.FrameN,WT.ReleaseFrame,1);
                             end
                         end    
                     case 'Touch_Continuous'
                         for k=1:Exp.TrN
                             for m=1:length(WOI)
                                WT=Whisker(WOI(m)).trial(k);
                                FeatY{k}=createBinary(WT.FrameN,WT.TouchFrame,WT.ReleaseFrame);
                             end
                         end 
                     otherwise
                         % get binned spikes
                         [SpikeBInT_all,~] = obj.SpikeBinInTrials(Exp,[]);
                         FeatY=SpikeBInT_all;
                 end
                 
                 %filter FeatY
                 FeatY_f=cell2mat(FeatY(TOI));
                 FeatY_f=logical(FeatY_f(filter_id));   %..........Spike binary in binned samples
                 
                 %construct tuning curve
                 h=histcounts(WhkFeatCat_discretized(FeatY_f),Edges);
                 h_all=histcounts(WhkFeatCat_discretized,Edges);
                 
                 YProb=h./(h_all);
                 if useRate
                    YProb=YProb*mean(Exp.videoFps);  %mean in case non consistent videoFps per trial
                    y_type='_Rate s-1';
                 else
                     y_type='_Prob';
                 end
                 %bar(BinCentres,SpkProb)
                 if isplot
                     hold on,
                     if ~isempty(c)
                         plot(BinCentres,YProb,'Color',c(i,:),'DisplayName',titleStr)
                     else
                         plot(BinCentres,YProb,'DisplayName',titleStr)
                     end
                     xlabel(Feat2TuneX)
                     ylabel([Feat2TuneY y_type])
                 end
                 Y(i,:)=YProb;
             end
             
         end
                  
         function [percentSpikes,x,Groups,Edges]=tuningcurve_old(obj,WhkBhvr,corrTime,binN,SpkTimeBin)
             %this function checks for spikes at time resolution so it is
             %very slow and computational heavy, use binned spikes instead,
             %i.e. tuningcurve()
             [Groups,Edges,corrTime] =DiscretizeBehvr(WhkBhvr,corrTime,binN);
             [percentSpikes,x]=obj.tuningcurveDiscretized(Groups,Edges,corrTime,SpkTimeBin);            
         end
         
         function [percentSpikes,x]=tuningcurveDiscretized_old(obj,Groups,Edges,corrTime,SpkTimeBin)
             x=Edges(1:end-1)+(Edges(2)-Edges(1))/2;
             percentSpikes=zeros(1,length(x));
%              for i=1:length(x)
%                  fprintf('Getting spikecounts for N%d bin%d\n',obj.CID,i)
%                  T1s=corrTime(Groups==i);
%                  [Iidx,~]=findClosest(T1s,obj.SpikeTime,SpkTimeBin);
%                  percentSpikes(i)=sum(~isnan(Iidx))/length(T1s);      
%              end
tic
            [Iidx,~]=findClosest(obj.SpikeTime,corrTime,SpkTimeBin);
            toc
            spkphase=Groups(Iidx(~isnan(Iidx)));
            
             for i=1:length(x)
                 fprintf('Getting spikecounts for N%d bin%d\n',obj.CID,i)
                 F_wSpks=sum(spkphase==i);
                 percentSpikes(i)=F_wSpks/sum(Groups==i);             
             end
         end
         function [STA_Bhvr,STA_range]=STA(obj,Exp,WhkBhvr,StimC,STArange)
             SpkT=obj.SpikeTime;
             bin=1/Exp.videoFps;
             if length(WhkBhvr)~=Exp.TrN
                 error('Wrong WhkBhvr format, need to be cell array of length Exp.TrN')
             end
              if isempty(STArange)
                STArange=[-0.05 0.1];
              end
             STArange=round(STArange*Exp.videoFps);
             windowL=STArange(end)-STArange(1)+1;
             STA_range=[STArange(1):STArange(end)]*bin;
             
             
             if ~isempty(StimC)
                 tr=StimC.validTrials;
             else
                 tr=1:Exp.TrN;
             end
             stackedB=[];
             for i=1:length(tr)
                 t=tr(i);
                 TDelay=max(Exp.PistonBuffer(t,:))/Exp.videoFps;
                 if isnan(TDelay)  %no piston out, i.e. free whisking, so all data points are used
                     TDelay=0;
                 end
                 alignTime=SpkT(SpkT>=Exp.FrameT{t}(1)+TDelay & SpkT<=Exp.FrameT{t}(end));
                 whkrb=WhkBhvr{t};
                 [Iidx,~]=findClosest(alignTime,Exp.FrameT{t},bin);
                 if sum(isnan(Iidx))>=1
                     123;
                 end
                 Iidx=Iidx(~isnan(Iidx));
                 offset=-1*windowL;
                whkrb_padded=NANpad(whkrb,length(whkrb)+windowL*2,offset);
                stackedB_tr=nan(length(Iidx),windowL);
                 for j=1:length(Iidx)
                     I=Iidx(j);
                     win=I+[STArange(1):STArange(end)]-offset;
                     stackedB_tr(j,:)=whkrb_padded(win);
                 end
                 stackedB=[stackedB;stackedB_tr];
             end
             STA_Bhvr=mean(stackedB,1,'omitnan');
             
         end
     end
end
classdef StimCondition
     properties( Access = public )
         % All your arguments here with their default values
           %to link to the corresponding exp parmas, imagine an array of neurons combined from all experiments with subset of them linking to a particular exp array
         WID;   %only these whiskers are touching pistons (vector)
         WID_FW; %sorted by touch irder, matrix TouchN x WhiskerN
         
         type;  %description of StimC (string)
         Wn;    %total no. of touching whisker
         WhiskerBool;  % whiskers touching (logical0
         
         Light=[];
         Hz=nan;
         Run=[];  %isempty=>not considered; (<=0) =>not run, th set as 2.0; (0< <1),is run, th set as obj.Run, else th=2.5
         ITIth;
         firstFew=[];  %only take the first few number of touches per trial, last few if a negative value is given
         FrameRange=[];
         
         
         % for now only PSTH are using behavior hits, if used more often
         % consider setting a parameter bool HITS with length=TrN
         
         validTrials;
         % touch info
         TouchTimes;  %double array
         ReleaseTimes;
         TouchFrames;   %cell array
         ReleaseFrames;
         isGoodPairs; %%cell array (trialN) of bool vector (TouchN) of WID(1)
         
         
         TouchDeltaF;
         TouchDeltaT; % (T2-T1)%matrix of length(TouchTimes) x NofPaired whiskers, empty if single whisker only
         ReleaseDeltaF;
         ReleaseDeltaT;% (R2-R1)
         CorrTrials; %array of length(TouchTimes)
         RunBool; %to be implemented, obtain runbool based off of Exp.RunSpeed for given StimC Run condition
         
         % touch features
         Toverlap;
%        TouchDur=toColumn(StimC.ReleaseTimes-StimC.TouchTimes);

         
         
         
     end
     methods
         function obj = StimCondition(Exp,Whisker,type,WID,Light,Hz,Run,ITIth,firstFew,FrameRange)   %WID is a vector of whiskers ID that are touching, and the post-analysis is in reference to the first whisker given
            obj.Wn=length(WID);
            obj.type=type;            
            if isvector(WID)
                obj.WID = WID;
                arr=zeros(1,Exp.Stim.Piston.tot);
                for i=1:obj.Wn
                    arr(WID(i))=1;
                end
                obj.WhiskerBool=arr;           
            elseif isempty(WID)
                obj.WID = WID;
                obj.WhiskerBool=zeros(1,Exp.Stim.Piston.tot);
            end
            if ~isempty(Light)
                obj.Light    = Light;
            end
            if ~isempty(Hz)
                obj.Hz    = Hz;
            end
            if ~isempty(Run)
                obj.Run    = Run;
            end
            obj.ITIth = ITIth;
            obj.firstFew=firstFew;
            obj.FrameRange = FrameRange;            
            
            obj=obj.findTouchTimes(Exp,Whisker);
         end
         
         function obj=findTouchTimes(obj,Exp,Whisker)
             if obj.Wn>=2  %there are more than two whiskers and their touchtimes need to be matched
                 obj=CloseTouchTimes(obj,Exp,Whisker);
             elseif obj.Wn==1  %only 1 whisker
                 if length(Whisker)==1
                    PW=Whisker;
                 else
                    PW=Whisker(obj.WID(1));
                 end
                 if isempty(obj.Light)
                    LightC=true(Exp.TrN,1);
                 elseif obj.Light
                    LightC=logical(Exp.Stim.Optic.f==obj.Hz);
                 else
                    LightC=isnan(Exp.Stim.Optic.f);
                 end
                 CorrespondingTrials= find(PistonComb(obj.WhiskerBool,Exp.Stim.Piston.Mat) & LightC);
                 CorrespondingTrials=CorrespondingTrials(~ismember(CorrespondingTrials,Exp.BadTrials));
                 obj.validTrials=CorrespondingTrials;
                 
                 TouchTimeByTrials=cell(length(CorrespondingTrials),1);
                 ReleaseTimeByTrials=cell(length(CorrespondingTrials),1);
                 TouchFs=cell(Exp.TrN,1);
                 ReleaseFs=cell(Exp.TrN,1);
                 GoodIDs=cell(Exp.TrN,1);
                 CorrTrialsByTrials=cell(length(CorrespondingTrials),1);
                 for j=1:length(CorrespondingTrials)
                     tr=CorrespondingTrials(j);
                     T1s=PW.trial(tr).TouchFrame;
                     R1s=PW.trial(tr).ReleaseFrame;
                     if ~isempty(obj.FrameRange) && length(obj.FrameRange)==2
                         FrameRangeBool=(T1s>obj.FrameRange(1) & R1s<=obj.FrameRange(2));
                     else
                         FrameRangeBool=true(length(T1s),1);
                     end
                     
                     % add run bool with length(T1s)
                     if ~isempty(obj.Run)
                         if obj.Run>0 && obj.Run<1  % is running, th based on given obj.Run if >0 & <1
                            runBool=(Exp.RunSpeed{tr}(T1s)>obj.Run);
                         elseif obj.Run==0  %is not running
                            runBool=(Exp.RunSpeed{tr}(T1s)<=0.25);
                         else % is running, th based on default 0.25 if obj.Run is out-of-bounds
                             runBool=(Exp.RunSpeed{tr}(T1s)>0.25);
                         end
                     else % do not consider running
                         runBool=true(length(T1s),1);
                     end
                     
                     %only get the first few touches (optional)
                     if ~isempty(obj.firstFew)
                         FirstFew=obj.firstFew;
                         %goodpairsFrame=goodpairsFrame(1:min([length(goodpairsFrame) FirstFew]));
                         if FirstFew>0
                             goodID=find(logical(runBool),FirstFew,'first');
                         elseif FirstFew<0 
                             goodID=find(logical(runBool),abs(FirstFew),'last');
                         end                  
                     else
                         goodID=logical(runBool);                      
                     end
                     
                     goodpairsTFrame=T1s(goodID & FrameRangeBool);
                     goodpairsRFrame=R1s(goodID & FrameRangeBool);
                     
                     TouchTimeByTrials{j,1}=toColumn(Exp.FrameT{tr}(goodpairsTFrame));
                     ReleaseTimeByTrials{j,1}=toColumn(Exp.FrameT{tr}(goodpairsRFrame));                     
                     TouchFs{tr,1}=toColumn(goodpairsTFrame);
                     ReleaseFs{tr,1}=toColumn(goodpairsRFrame);
                     CorrTrialsByTrials{j,1}=ones(length(goodpairsTFrame),1)*tr;
                     GoodIDs{tr,1}=goodID;
                     
%                      g=sprintf('%d ', obj.WID);
%                      fprintf('Pair touches for W %strial %d: %d/%d\n',g,tr,length(TouchTimeByTrials{j,1}),length(T1s))

                     
                 end
                 obj.TouchTimes=cell2mat(TouchTimeByTrials);
                 obj.ReleaseTimes=cell2mat(ReleaseTimeByTrials);
                 obj.TouchFrames=TouchFs;
                 obj.ReleaseFrames=ReleaseFs;
                 obj.isGoodPairs=GoodIDs;
                 obj.CorrTrials=cell2mat(CorrTrialsByTrials);
             else 
                 if isempty(obj.Light)
                    LightC=true(Exp.TrN,1);
                 elseif obj.Light
                    LightC=logical(Exp.Stim.Optic.f==obj.Hz);
                 else
                    LightC=isnan(Exp.Stim.Optic.f);
                 end
                 CorrespondingTrials= find(PistonComb(obj.WhiskerBool,Exp.Stim.Piston.Mat) & LightC);
                 CorrespondingTrials=CorrespondingTrials(~ismember(CorrespondingTrials,Exp.BadTrials));
                 obj.validTrials=CorrespondingTrials;
             end
         end
         
         
         function obj=CloseTouchTimes(obj,Exp,Whisker)
             if obj.Wn<2 || length(Whisker)<2
                 error('At least 2 whiskers are needed to find CloseTouchTimes')
             end
             if isempty(obj.ITIth) %no ITI
                 th=[];
             else
                 th=obj.ITIth.*Exp.videoFps;
             end
             
             wid=obj.WID;
             PW=Whisker(wid(1));
             WhiskersLeft=obj.Wn-1;
             
             PairTouchTimes=cell(Exp.TrN,1);
             PairReleaseTimes=cell(Exp.TrN,1);
             TouchFs=cell(Exp.TrN,1);
             ReleaseFs=cell(Exp.TrN,1);             
             PairdeltaF=cell(Exp.TrN,WhiskersLeft);
             PairRdeltaF=cell(Exp.TrN,WhiskersLeft);
             CorrTrialsByTrials=cell(Exp.TrN,1);
             GoodIDs=cell(Exp.TrN,1);
             %determine light bool
             if isempty(obj.Light)
                 LightC=true(Exp.TrN,1);
             elseif obj.Light
                 LightC=logical(Exp.Stim.Optic.f==obj.Hz);
             else
                 LightC=isnan(Exp.Stim.Optic.f);
             end
             CorrespondingTrials= find(PistonComb(obj.WhiskerBool,Exp.Stim.Piston.Mat) & LightC);
             CorrespondingTrials=CorrespondingTrials(~ismember(CorrespondingTrials,Exp.BadTrials));
             obj.validTrials=CorrespondingTrials;
             
             
             
             for i=1:Exp.TrN 
                 if ~ismember(i,CorrespondingTrials)
                     continue;
                 end           
                 T1s=PW.trial(i).TouchFrame;
                 R1s=PW.trial(i).ReleaseFrame;   
                 if ~isempty(obj.FrameRange) && length(obj.FrameRange)==2
                     FrameRangeBool=(T1s>obj.FrameRange(1) & R1s<=obj.FrameRange(2));
                 else
                     FrameRangeBool=true(length(T1s),1);
                 end
                 
                 
                 goodpairsID=false(length(T1s),WhiskersLeft);
                 Iidx=cell(1,WhiskersLeft);

                 deltaF=cell(1,WhiskersLeft);
                 deltaFR=cell(1,WhiskersLeft);
                 
                 for PPWid=1:WhiskersLeft
                     PPW=Whisker(wid(PPWid+1));
                     T2s=PPW.trial(i).TouchFrame;
                     R2s=PPW.trial(i).ReleaseFrame;
                     if ~isempty(obj.FrameRange) && length(obj.FrameRange)==2
                         I=(T2s>obj.FrameRange(1) & R2s<=obj.FrameRange(2));
                         T2s=T2s(I);
                         R2s=R2s(I);
                     end
                     
                     
                     
                     if isempty(th) %no ITI
                         threshold=[];
                     elseif length(th)==1  %all trials have fixed videoFps, i.e. length(Exp.videoFps)==1
                        threshold=th;
                     else   %(rare unless there is unstable frame drops) variable fps,i.e. length(Exp.videoFps)==Exp.TrN
                        threshold=th(i);
                     end
%                      if ~isempty(T2s)
%                          123
%                      end
                     [Iidx{PPWid},deltaF{PPWid}]=findClosest(T1s,T2s,threshold);  %isCloseT2 also gives you the corresponding pair touch on the other whisker if needed
                     A=Iidx{PPWid};
                     [a,b]=findgroups(A);
                     notRepeatingTouch=true(length(A),1);
                     for m=1:length(b)
                         groupm=find(a==m);
                         if length(groupm)>1
                             [~,minI]=min(abs((deltaF{PPWid}(groupm))));
                             notRepeatingTouch(groupm)=0;
                             notRepeatingTouch(groupm(minI))=1;          
                         end
                     end
                     
                     goodpairsID(:,PPWid)=(~isnan(Iidx{PPWid}) & notRepeatingTouch);   %bool good pairs out of all touches                  
                     deltaFR{PPWid}=ones(length(T1s),1);
                     try
                     deltaFR{PPWid}(~goodpairsID(:,PPWid))=nan;
                     deltaFR{PPWid}(goodpairsID(:,PPWid))=R2s(Iidx{PPWid}(~isnan(Iidx{PPWid})))-R1s(goodpairsID(:,PPWid));
                     catch
%                          123
                     end
                 
                 end  
                 
                 goodMultipairsID=logical(sum(goodpairsID,2)==WhiskersLeft);  %all columns(i.e. all paired whiskers) should be one
                 
                 % add run bool with length(T1s)
                 if ~isempty(obj.Run)
                     if obj.Run>0 && obj.Run<1  % is running, th based on given obj.Run if >0 & <1
                        runBool=(Exp.RunSpeed{i}(T1s)>obj.Run);
                     elseif obj.Run==0  %is not running
                        runBool=(Exp.RunSpeed{i}(T1s)<=0.25);
                     else % is running, th based on default 0.25 if obj.Run is out-of-bounds
                         runBool=(Exp.RunSpeed{i}(T1s)>0.25);
                     end
                 else % do not consider running
                     runBool=true(length(T1s),1);
                 end

                 %only get the first few touches (optional)
                 if ~isempty(obj.firstFew)
                     FirstFew=obj.firstFew;
                     if FirstFew>0
                         goodID=find(logical(goodMultipairsID & runBool & FrameRangeBool),FirstFew,'first');
                     elseif FirstFew<0
                         goodID=find(logical(goodMultipairsID & runBool & FrameRangeBool),abs(FirstFew),'last');
                     end
                     goodpairsTFrame=T1s(goodID);
                     goodpairsRFrame=R1s(goodID);
                     
                 else
                     goodID=logical(goodMultipairsID & runBool & FrameRangeBool);
                     goodpairsTFrame=T1s(goodID);
                     goodpairsRFrame=R1s(goodID);                     
                 end
                 
                 PairTouchTimes{i,1}=toColumn(Exp.FrameT{i}(goodpairsTFrame));
                 PairReleaseTimes{i,1}=toColumn(Exp.FrameT{i}(goodpairsRFrame));
                 
                 TouchFs{i,1}=goodpairsTFrame;
                 ReleaseFs{i,1}=goodpairsRFrame;

                 CorrTrialsByTrials{i,1}=ones(length(goodpairsTFrame),1)*i;
                 
                 for PPWid=1:WhiskersLeft   
                     PairdeltaF{i,PPWid}=toColumn(deltaF{PPWid}(goodID));
                     PairRdeltaF{i,PPWid}=toColumn(deltaFR{PPWid}(goodID));
                 end  
                 GoodIDs{i,1}=goodID;
%                     g=sprintf('%d ', obj.WID);
%                      fprintf('Pair touches for W %strial %d: %d/%d\n',g,i,length(PairTouchTimes{i,1}),length(T1s))
 
             end
             
             obj.TouchTimes=cell2mat(PairTouchTimes);
             obj.TouchFrames=TouchFs;
             obj.ReleaseTimes=cell2mat(PairReleaseTimes);
             obj.ReleaseFrames=ReleaseFs;
             obj.isGoodPairs=GoodIDs;
             obj.TouchDeltaF=cell2mat(PairdeltaF);
             obj.ReleaseDeltaF=cell2mat(PairRdeltaF);
             obj.CorrTrials=cell2mat(CorrTrialsByTrials);
             if length(Exp.videoFps)==1
                obj.TouchDeltaT=obj.TouchDeltaF/Exp.videoFps;
                obj.ReleaseDeltaT=obj.ReleaseDeltaF/Exp.videoFps;
             else
                obj.TouchDeltaT=obj.TouchDeltaF./Exp.videoFps(obj.CorrTrials);
                obj.ReleaseDeltaT=obj.ReleaseDeltaF./Exp.videoFps(obj.CorrTrials);
             end
         end
         
         function obj=updateTouches(obj,I)
             % update touches based on given logical or array I
             if islogical(I)
                 I=find(I);
             end
             I=sort(I);

             TF=obj.TouchFrames;      
             RF=obj.ReleaseFrames;
             goodID=obj.isGoodPairs;
             
             count=0;
             for i=1:length(TF)
                 if ~isempty(TF{i})
                     TN=length(TF{i});
                     goodTinTr=ones(TN,1);
                     previousGoodTouchID=find(goodID{i});
                     for j=1:TN
                         if ~ismember(count+j,I)
                             goodTinTr(j)=0;
                             goodID{i}(previousGoodTouchID(j))=0;
                         end
                     end
                     TF{i}=TF{i}(logical(goodTinTr));
                     RF{i}=RF{i}(logical(goodTinTr));
                     count=count+TN;
                 end
             end
             %check data integrity:
             if count~=length(obj.TouchTimes)
                 error('Original TouchFrame Number is not equal to Original TouchTimes Number, check code or redetermine TouchFrame with updated TouchTimes & Exp params')
             end
             obj.TouchTimes=obj.TouchTimes(I);
             obj.ReleaseTimes=obj.ReleaseTimes(I);
             obj.CorrTrials=obj.CorrTrials(I);
             if ~isempty(obj.TouchDeltaT)
                 obj.TouchDeltaT=obj.TouchDeltaT(I,:);
                 obj.ReleaseDeltaT=obj.ReleaseDeltaT(I,:);
                 obj.TouchDeltaF=obj.TouchDeltaF(I,:);
                 obj.ReleaseDeltaF=obj.ReleaseDeltaF(I,:);
             end
             obj.TouchFrames=TF;
             obj.ReleaseFrames=RF;
             obj.isGoodPairs=goodID;
             obj.validTrials=unique(obj.CorrTrials);
             if any(strcmp(properties(obj), 'WID_FW')) && ~isempty(obj.WID_FW) && ~all(isnan(obj.WID_FW),'all')
                obj.WID_FW=obj.WID_FW(I,:);
             end
             if any(strcmp(properties(obj), 'Toverlap')) && ~isempty(obj.Toverlap) && ~all(isnan(obj.Toverlap),'all')
                obj.Toverlap=obj.Toverlap(I,:);
             end
         end
         
         function StimC_adjusted=PWadjusted(obj,Exp,Whisker,Neuron)
%              Recalculate all StimC properties based on given neurons' strongest whisker response
%              and gives StimC_adjusted
 
             wid=obj.WID;
            
             NeuN=length(Neuron);
%              TN=length(obj.TouchTimes);
%              TT=cell(1,NeuN);
%              DT=cell(1,NeuN);
             StimC_adjusted=StimCondition.empty(NeuN,0);
               for i=1:NeuN
                   N_PWs=Neuron(i).PW;
%                     N_PWs=N_PWs(ismember(N_PWs,wid));
                   if isempty(wid)  %free whisking
                       StimC_adjusted(i)=obj;
                       continue
                   end
                   [allRmembers,id]=ismember(wid,N_PWs);
                   %check integrity
                   if ~allRmembers
                       error('Unidentified StimC wid from neuron principal whiskers')
                   end
                   [~,I]=sort(id);
                   obj_adj = StimCondition(Exp,Whisker,strcat(obj.type,"_PW"),wid(I),obj.Light,obj.Hz,obj.Run,obj.ITIth,obj.firstFew,obj.FrameRange);
                                      
                   
                   
                   %                  obj = StimCondition(Exp,Whisker,type,WID,Light,Hz,Run,ITIth,firstFew,FrameRange)   %WID is a vector of whiskers ID that are touching, and the post-analysis is in reference to the first whisker given

%                    dt=[zeros(TN,1) obj.TouchDeltaT];
%                    dt=dt(:,I);
%                    dtPW=dt(:,1);
%                    dt=dt-dtPW;
%                    dt=dt(:,2:end);
%                    TT{i}=obj.TouchTimes+dtPW;
%                    DT{i}=dt;

                   %create new adjusted StimC (not completed, still need releasetimes and Frames)
%                    obj_adj=obj;
%                    obj_adj.WID=wid(I);
%                    obj_adj.TouchTimes=TT{i};
%                    obj_adj.TouchDeltaT=DT{i};
                   
                   StimC_adjusted(i)=obj_adj;
                       
               end
             

         end
         
         function StimC_adjusted=FWadjusted(obj,Exp)
            
             if length(Exp.videoFps)~=1
                 error('Incompatible code with variable videoFps')
             end
             
             wid=obj.WID;
             
             if length(wid)==1
                 StimC_adjusted=obj;
                 return
             end
             
             dF=obj.TouchDeltaF;
             dFr=obj.ReleaseDeltaF;
             TouchN=size(dF,1);
             
             Trials=obj.CorrTrials;
             uniqueTrials=unique(Trials);
             
             [rearranged_dF,I_dF]=sort([zeros(TouchN,1) dF],2,'ascend');
             newWID=wid(I_dF);
             newdF=rearranged_dF-rearranged_dF(:,1);
             newTouchDeltaF=newdF(:,2:end);
             newTouchDeltaT=newTouchDeltaF*1/Exp.videoFps;
             
             rearranged_newdFr=[zeros(TouchN,1) dFr];
             for i=1:TouchN
                rearranged_newdFr(i,:)=rearranged_newdFr(i,I_dF(i,:));
             end 
             newdFr=rearranged_newdFr-rearranged_newdFr(:,1);
             newReleaseDeltaF=newdFr(:,2:end);
             newReleaseDeltaT=newReleaseDeltaF*1/Exp.videoFps;
             
             
             
             newTouchFrames=cell(Exp.TrN,1);
             newTouchTimes=zeros(TouchN,1);
             newReleaseFrames=cell(Exp.TrN,1);
             newReleaseTimes=zeros(TouchN,1);
             
             count=0;
             for i=1:length(uniqueTrials)
                 
                 tr=uniqueTrials(i);
                 TNinTr=sum(Trials==tr);                
                 TF=obj.TouchFrames{tr};
                 RF=obj.ReleaseFrames{tr};
                 
                 if TNinTr~=length(TF)
                     error('Inconsistent TouchN, corrupted')
                 end
                 
                 for j=1:TNinTr
                     count=count+1;
                     FrameN_adjustment=rearranged_dF(count,1);
                     TF(j)=TF(j)+FrameN_adjustment;
                     TT=Exp.FrameT{tr}(TF(j));
                     newTouchTimes(count)=TT;
                     
                     FrameN_adjustment=rearranged_newdFr(count,1);
                     RF(j)=RF(j)+FrameN_adjustment;
                     RT=Exp.FrameT{tr}(RF(j));
                     newReleaseTimes(count)=RT;
                 end              
                 newTouchFrames{tr}=TF;
                 newReleaseFrames{tr}=RF;    
             end
             if count~=TouchN
                     error('Inconsistent TouchN, corrupted')
             end
                 
             StimC_adjusted=obj;
             StimC_adjusted.TouchDeltaT=newTouchDeltaT;
             StimC_adjusted.ReleaseDeltaT=newReleaseDeltaT;
             StimC_adjusted.TouchDeltaF=newTouchDeltaF;
             StimC_adjusted.ReleaseDeltaF=newReleaseDeltaF;
             StimC_adjusted.TouchTimes=newTouchTimes;
             StimC_adjusted.ReleaseTimes=newReleaseTimes;
             StimC_adjusted.TouchFrames=newTouchFrames;
             StimC_adjusted.ReleaseFrames=newReleaseFrames;
             StimC_adjusted.WID_FW=newWID;
             StimC_adjusted.type=strcat(obj.type,"_FW");
             
         end
         
         function touchraster(obj,Exp,color)
%              skipEmpty=0; %only affect data output format
             
             if isempty(color)
                 color='k';
             end
             TT=obj.TouchTimes;
             CorrT=obj.CorrTrials;
             UniqTrials=unique(obj.validTrials);

             data=cell(1,length(UniqTrials));
             gap=0.1;firstSegDur=1;lastSegDur=1;
             for i=1:length(UniqTrials)
                 TTi=TT(CorrT==UniqTrials(i));
                 firstSeg=TTi(TTi<=Exp.FrameT{UniqTrials(i)}(1)+firstSegDur)-Exp.FrameT{UniqTrials(i)}(1);
                 lastSeg=TTi(TTi>=Exp.FrameT{UniqTrials(i)}(end)-lastSegDur)-(Exp.FrameT{UniqTrials(i)}(end)-lastSegDur)+firstSegDur;
                 data{i}=[toColumn(firstSeg);toColumn(lastSeg+gap)];
             end
             rasterplot(data,color,[],[0 firstSegDur+lastSegDur+gap])
             addline('x',firstSegDur,'k','--');
             addline('x',firstSegDur+gap,'k','--');
             addline('x',Exp.Stim.Piston.Delay-Exp.Cam.Delay,'r','--');
             addline('x',firstSegDur+lastSegDur+gap-(Exp.Stim.Piston.OffWin-Exp.Cam.OffWin),'r','--');
         end
         function [x,y,STD,n,mat]=PSTH(obj,Exp,Whisker,type,alignmentType,hits,PSTHrange,yrange)
             %output:
%              x: Frame vector
%              y: mean type value at each frame
%              STD: std
%              n: sample size
%              mat: sample population data matrix

%              hits: filters out trials
%             TrialLimit=1:13;
%             TrialLimit=13:25;
            if length(Exp.videoFps)~=1
                fps=Exp.OldFps;
                disp('Variable videoFps, using Exp.OldFps')
            else
                fps=Exp.videoFps;
            end
             CameraStart=0;
             
             Hz30Delay=(Exp.Stim.Optic.Delay-Exp.Cam.Delay+0.5)*fps;
             
             OpticDelay=(Exp.Stim.Optic.Delay-Exp.Cam.Delay)*fps;
             PistonDelay=(Exp.Stim.Piston.Delay-Exp.Cam.Delay)*fps;
             WaterDelay=(2.15-Exp.Cam.Delay)*fps;  % 0.65 PistonDelay + 1.5s - 0.15 CamDelay wrt Trialstart
             touchFs=obj.TouchFrames;
             try
             pulseOnset=Exp.Stim.Optic.OnsetF;
             pulseOffset=Exp.Stim.Optic.OffsetF;
             catch
             end
             %              zeroIndex=find(PSTHrange==0);
             FramesBeforeTp=round(PSTHrange(1)*fps);
             FramesAfterTp=round(PSTHrange(end)*fps);
             if length(PSTHrange)>2
                 binSize=PSTHrange(2)-PSTHrange(1);
             else
                 binSize=1/fps;
             end
             WkrN=length(Whisker);
             WhkrPlots=cell(WkrN,1);
             trials=unique(obj.validTrials);
             if ~isempty(hits)
                 if length(hits)==Exp.TrN  %full Exp.TrN array length      %&& sum(hits<0 & hits>1)==0
                    trials=trials(hits(trials)>0);
                 else 
                     trials=hits;
                 end
             end
%              disp(trials)
%              timepoints=obj.TouchTimes;

            %limit trials to a certain number
            if exist('TrialLimit','var') && ~isempty(TrialLimit)
                trials=trials(TrialLimit);
            end
            
             for i=1:WkrN
%                  activeW=find(Exp.Stim.Piston.active);
                 W=Whisker(i);
                 WhkrPlot=[];
                 for j=1:length(trials)
                     WT=W.trial(trials(j));
                     switch type
                         case 'Angle'
                             WkrBhvr=WT.CenAngle.sig;
                             israte=0;
                         case 'Phase'
                             WkrBhvr=WT.Phase;
                             israte=0;
                         case 'Curvature'
                             WkrBhvr=WT.Curvature;
                             israte=0;
                         case 'CurvDebase'
                             WkrBhvr=WT.CurvDebase;
                             israte=0;
                         case 'Amplitude'
                             WkrBhvr=WT.Amplitude;
                             israte=0;
                         case 'Setpoint'
                             WkrBhvr=WT.Setpoint;
                             israte=0;
                         case 'SetpointNorm'
                             WkrBhvr=WT.Setpoint;
                             WkrBhvr=minmaxnorm(WkrBhvr);
                             israte=0;
                         case 'Sync'
                             WkrBhvr=Exp.Sync{trials(j)};
                             israte=0;
                         case 'Directionality'
                             WkrBhvr=Exp.Directionality{trials(j)};
                             israte=0;
                         case 'Spatiality'
                             WkrBhvr=Exp.Spatiality{trials(j)};
                             israte=0;
                         case 'RunSpeed'
                             WkrBhvr=Exp.RunSpeed{trials(j)};
                             israte=0;
                         case 'RunSpeedcm'
                             WkrBhvr=Exp.Runcmpers{trials(j)};
                             israte=0;
                         case 'Touch'
                             WkrBhvr=createBinary(WT.FrameN,obj.TouchFrames{trials(j)},1);
                             israte=1;
                         otherwise
                             error('Type %s not recognised.',type);
                     end

                    switch alignmentType  %in frames
                        case 'Touch'
                            tp=touchFs{trials(j)};
                        case 'FirstTouch'
                            try
                                tp=touchFs{trials(j)}(1);
                            catch
                                tp=[];
                            end
                        case 'TouchBeforeWater'
                            tp=touchFs{trials(j)};
                            tp=tp(tp<=WaterDelay);
                        case 'Piston'
                            tp=round(PistonDelay);
                        case 'Light'
                            tp=round(OpticDelay);
                        case 'PulseOnset'
                            OpticOnsetFsize=fps/Exp.Stim.Optic.f(trials(j));
                            tp=round(OpticDelay):OpticOnsetFsize:WT.FrameN;%round(pulseOnset{trials(j)});  
                        case 'Pseudo10HzPulseOnsetWhenRun'
                            OpticOnsetFsize=fps/10;
                            tp=round(OpticDelay):OpticOnsetFsize:WT.FrameN;%round(pulseOnset{trials(j)});  
                            runBool=(Exp.RunSpeed{trials(j)}>=0.25);
                            tp=tp(logical(runBool(tp)));
                        case 'Pseudo20HzPulseOnsetWhenRun'
                            OpticOnsetFsize=fps/20;
%                             tp=round(round(OpticDelay):OpticOnsetFsize:Hz30Delay);%round(pulseOnset{trials(j)});  
                            tp=round(round(OpticDelay):OpticOnsetFsize:WT.FrameN);%round(pulseOnset{trials(j)});  
                            runBool=(Exp.RunSpeed{trials(j)}>=0.25);
                            tp=tp(logical(runBool(tp)));
                        case 'Pseudo30HzPulseOnsetWhenRun'
                            OpticOnsetFsize=fps/30;
                            tp=round(round(Hz30Delay):OpticOnsetFsize:WT.FrameN);%round(pulseOnset{trials(j)});  
                            runBool=(Exp.RunSpeed{trials(j)}>=0.5);
                            tp=tp(logical(runBool(tp)));
                        case 'Pseudo30HzPulseOnsetNoRun'
                            OpticOnsetFsize=fps/30;
                            tp=round(round(Hz30Delay):OpticOnsetFsize:WT.FrameN);%round(pulseOnset{trials(j)});  
                            runBool=(Exp.RunSpeed{trials(j)}<0.25);
                            tp=tp(logical(runBool(tp)));
                        case 'PulseOnsetWhenNoRun'
                            OpticOnsetFsize=fps/Exp.Stim.Optic.f(trials(j));
                            tp=round(OpticDelay):OpticOnsetFsize:WT.FrameN;%round(pulseOnset{trials(j)});  
                            runBool=(Exp.RunSpeed{trials(j)}<0.25);
                            tp=tp(logical(runBool(tp))); 
                            
                        case 'PulseOffset'
                            tp=round(pulseOffset{trials(j)});     
                        case 'Camera'
                            tp=round(CameraStart);
                        case 'Water'
                            tp=round(WaterDelay);
                        case 'SyncReset'
                            tp=Exp.ResetTime{trials(j)};
                        case 'SyncDrift'
                            tp=Exp.DriftTime{trials(j)};
                            

                        otherwise
                            error('Alignment Type %s not recognised.',alignmentType);
                    end
                    if isempty(tp)
                        fprintf('Aligntype is empty for StimC: %s;Trial %d\n',obj.type,trials(j))
                    end
                    
                     WP=[];  
                     for k=1:length(tp)
                         try
                            WP_videoFps=WkrBhvr(tp(k)+FramesBeforeTp:tp(k)+FramesAfterTp);         
                         catch
                             frontNAN=nan(abs(min([0 tp(k)+FramesBeforeTp-1])),1);
                             backNAN=nan(max([0 tp(k)+FramesAfterTp-length(WkrBhvr)]),1);
                             WP_videoFps=[frontNAN;WkrBhvr(max([1 tp(k)+FramesBeforeTp]):min([tp(k)+FramesAfterTp length(WkrBhvr)]));backNAN];
                         end
                        if binSize~=1/fps
                            [temp,~] =  histcounts(find(WP_videoFps~=0 & ~isnan(WP_videoFps))/fps,[toColumn(PSTHrange);PSTHrange(end)+binSize]);
                            if israte
                                WP(k,:)=temp/binSize;%(2:end);
                            else
                                error('Need code writing to resample non-binary data. If possible just use videoFps')
%                                 WP(k,:)=temp;
                            end
%                             WP(k,:)=resample(WP_videoFps,1/binSize,Exp.videoFps);
                        else
                            WP(k,:)=WP_videoFps;
                        end 
                         
                     end
                     
                     
                     
                     
                     WhkrPlot=[WhkrPlot;WP];
                 end
                 WhkrPlots{i}=WhkrPlot;
                 sizePerColumn{i}=sum(~isnan(WhkrPlots{i}));
                 STD{i}=std(WhkrPlots{i},'omitnan');
                 SEM{i}=STD{i}./(sizePerColumn{i}.^0.5);
             end
            %only look at one given whisker since touchtimes are
            %based on it, defined as whisker K
            
%             if ~exist('K','var')
%                 if ~isempty(obj.WID)
%                     K=obj.WID(1);  %look at the main touching whisker
%                 else
%                     K=1;  %no whisker touching, look at W1
%                 end
%             end
            x=[FramesBeforeTp*binSize:binSize:FramesAfterTp*binSize];
            for k=1:WkrN   
%                 hold on
%                 plot(FramesBeforeTp:FramesAfterTp,WhkrPlots{i}','color',[0,0,0]+0.3)
%                 x=FramesBeforeTp:FramesAfterTp;
                
                try
                    y{k}=mean(WhkrPlots{k},1,'omitnan');
                catch
                    error('fu')
                end
%                 STD{k}=STD{k};
                y1{k}=y{k}-SEM{k};
                y2{k}=y{k}+SEM{k};
%                 plotWFilledError(x,y,[0 0 0],y1,y2,[0.3 0.3 0.3],0.3)
%                 plotWFilledError(x,y,[1 0 0],y1,y2,[0.8 0.3 0.3],0.3)
                n{k}=sizePerColumn{k};
                mat{k}=WhkrPlots{k};
%                 fprintf('n = %d',n)
%                 plot(FramesBeforeTp:FramesAfterTp,mean(WhkrPlots{i},1,'omitnan'),'r')
%                 if ~isempty(yrange)
%                     ylim(yrange)
%                 end
            end
            
         end
         
         function [PPP,ff,t]=TrialAveragedPowerSpectrum(obj,Exp,Whisker,type)
             if length(Exp.videoFps)~=1
                 error('Incompatible code with variable videoFps')
             end
             Tres=0.2;
             f=5:1:40;
             win=[];
             windowEnd=round(1.3*Exp.videoFps);
             
             WkrN=length(Whisker);
             PPP=cell(WkrN,1);
             tr=obj.validTrials;
             
%              roughly determine runbool, >90% running during trial
runWin=round(Exp.Stim.Piston.Delay*Exp.videoFps);
             runbool=zeros(length(tr),1);
             for i=1:length(tr)
                runbool(i)=(sum(Exp.RunSpeed{tr(i)}(1:runWin)>=0.25)/runWin)>=0.9;
             end
             tr=tr(logical(runbool));
             fprintf('Trials running: %d/%d\n',sum(runbool),length(runbool))
             
             PistonDelay=round((Exp.Stim.Piston.Delay-Exp.Cam.Delay)*Exp.videoFps);
             OpticDelay=round((Exp.Stim.Optic.Delay-Exp.Cam.Delay)*Exp.videoFps);
             Hz30Delay=round((Exp.Stim.Optic.Delay-Exp.Cam.Delay+0.5)*Exp.videoFps);
             
             
             
             PXX_BL=cell(WkrN,1);
             PXX_Stim=cell(WkrN,1);
             
             
             for i=1:WkrN
                 BaselineX=[];
                 StimX=[];
                 AllX=[];
                 RunX=[];
                 chkpt_tr=zeros(length(tr),2);
                 chkpt_bl=zeros(length(tr),1);
                 
                 W=Whisker(i);
                 Pxx_BL=zeros(length(tr),length(f));
                 Pxx_Stim=zeros(length(tr),length(f));
                 
                 PP=cell(length(tr),1);
                 for j=1:length(tr)
                     WT=W.trial(tr(j));
                     x=WT.CenAngle.sig;
                     switch type
                         case 'Piston'
                             BaselineWin=1:PistonDelay;
                             StimWin=PistonDelay:WT.FrameN-round(Exp.Stim.Piston.OffWin*Exp.videoFps);
                         case 'Light'
                              BaselineWin=1:OpticDelay;
                             StimWin=OpticDelay:WT.FrameN-round(Exp.Stim.Optic.OffWin*Exp.videoFps);
                         case '30Hz'
                              BaselineWin=1:OpticDelay;
                             StimWin=Hz30Delay:WT.FrameN-round(Exp.Stim.Piston.OffWin*Exp.videoFps);
                         case '20Hz'
                              BaselineWin=1:OpticDelay;
                             StimWin=OpticDelay:Hz30Delay;
                     end
                     
%                      %run filter
%                      runBool=(Exp.RunSpeed{tr(j)}>=0.25);
%                      BaselineWin=BaselineWin(logical(runBool(BaselineWin)));
%                      StimWin=StimWin(logical(runBool(StimWin)));
                     try
                     XB=toColumn(x(BaselineWin));
                     XS=toColumn(x(StimWin));
                     catch
                         error('Window based on type %s are out of camera bounds',type)
                     end
%                      XB=bandpass(XB,[5 50],Exp.videoFps);
%                      XS=bandpass(XS,[5 50],Exp.videoFps);

%                      try
%                         [pxx_bl,~] = pwelch(XB,[],[],f,Exp.videoFps);
%                      catch  %signal too short
%                          pxx_bl=nan(1,length(f));
%                          disp('yikes')
%                      end
%                      try
%                      [pxx_Stim,~] = pwelch(XS,[],[],f,Exp.videoFps);
%                      catch
%                          pxx_Stim=nan(1,length(f));
%                          disp('eww')
%                      end
%                      Pxx_BL(j,:)=pxx_bl;
%                      Pxx_Stim(j,:)=pxx_Stim;
                    BaselineX=[BaselineX;XB];
                    StimX=[StimX;XS];
                    chkpt_bl(j,1)=length(AllX)+length(XB);
                    chkpt_tr(j,1)=length(AllX);
                    chkpt_tr(j,2)=tr(j);
                    AllX=[AllX;XB;XS];
                    trialX=[XB;XS];
                    try
                    trialX=trialX(1:windowEnd);
                    catch
                        if windowEnd>=min(Exp.FrameN)
                            error('WindowEnd(%d) is out of bounds(%d) for trial %d, likely due to short trials',windowEnd,length(trialX),tr(j))
                        else
                            error('WindowEnd(%d) is out of bounds(%d) for trial %d, likely due to early stimWinOff',windowEnd,length(trialX),tr(j))
                        end
                    end
                    
                    trialX=bandpass(trialX,[4 40],Exp.videoFps);
                    [p,ff,t]=pspectrum(trialX,Exp.videoFps,'spectrogram', ...
    'FrequencyLimits',[f(1) f(end)],'TimeResolution',Tres,'OverlapPercent',95);
                    PP{j}=p;
                    RunX{j}=toColumn(Exp.RunSpeed{tr(j)}(1:windowEnd));
                    
                 end
%                  BaselineX=bandpass(BaselineX,[5 50],Exp.videoFps);
%                  StimX=bandpass(StimX,[5 50],Exp.videoFps);
                 
%                  [pxx_bl,~] = pwelch(BaselineX,[],[],f,Exp.videoFps);
%                  [pxx_Stim,~] = pwelch(StimX,win,[],f,Exp.videoFps);
% figure(i)
% close
% figure(i)
% pspectrum(AllX,Exp.videoFps,'spectrogram', ...
%     'FrequencyLimits',[f(1) f(end)],'TimeResolution',0.5);

PPP_i=zeros(size(PP{1},1),size(PP{1},2));
RunAvg=zeros(length(RunX{1}),1);
for k=1:length(PP)
    PPP_i=PPP_i+PP{k};
    RunAvg=RunAvg+RunX{k};
end
PPP{i}=PPP_i/length(PP);
RunAvg=RunAvg/length(PP);

[X,Y] = meshgrid(t,ff);
% s=surf(X,Y,10*log(PPP));
figure
subplot(4,2,1)

for k=1:8
subplot(4,2,k)
switch k
    case 1
        PSD=PPP{i};
        Runplot=RunAvg;
        str=sprintf('Trial Average n=%d',length(PP));
    otherwise
        try
        PSD=PP{k-1};
        Runplot=RunX{k-1};
        str=sprintf('Trial: %d',tr(k-1));
        
        catch
            break
        end
end


s=surf(X,Y,10*log(PSD));

s.EdgeColor = 'none';
view(0,90)
caxis([-80 -20])
addline('x',PistonDelay/Exp.videoFps,'r','-');
if k>1
    TTw=Whisker(i).trial(tr(k-1)).TouchFrame;
    TTw=TTw(TTw<windowEnd);
    for m=1:length(TTw)
        addline('x',TTw(m)/Exp.videoFps,'k','-');
    end
    if ~isempty(obj.TouchFrames)
        TT=obj.TouchFrames{tr(k-1)};
        TT=TT(TT<windowEnd);
        for m=1:length(TT)
            addline('x',TT(m)/Exp.videoFps,'r',':');
        end
    end
end

yyaxis right
ylim([0 1])
plot(0:1/Exp.videoFps:(windowEnd-1)/Exp.videoFps,Runplot)
title(str)
end
suptitle(sprintf('%s Whisker %d',obj.type,Whisker(i).W_id))
% 
% xx=xlim;
% if xx(2)<30
%     scale=60;
% else
%     scale=1;
% end
% 
% yyaxis right
% plot(0:1/Exp.videoFps/scale:(length(RunX)-1)/Exp.videoFps/scale,RunX)
% ylim([0 1])
% for k=1:size(chkpt_tr,1)
%     hold on
%     addline('x',chkpt_tr(k,1)/Exp.videoFps/scale,'r',[])
% end
% for k=1:size(chkpt_bl,1)
%     hold on
%     addline('x',chkpt_bl(k,1)/Exp.videoFps/scale,'r','--')
% end
%                  pspectrum(BaselineX,Exp.videoFps,'spectrogram', ...
%     'FrequencyLimits',[f(1) f(end)],'TimeResolution',0.1)
%                     pspectrum(StimX,Exp.videoFps,'spectrogram', ...
%     'FrequencyLimits',[f(1) f(end)],'TimeResolution',0.5)


%                  %concatenate
%                  PXX_BL{i}=[];
%                  PXX_Stim{i}=[];
                 %trial average
%                  PXX_BL{i}=Pxx_BL;
%                  PXX_Stim{i}=Pxx_Stim;               
             end
             
         end
         
         
     end
end
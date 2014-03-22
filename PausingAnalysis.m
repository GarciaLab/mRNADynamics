function PausingAnalysis

%Runs PlotTraces for different datasets. 


%% Set folders, load the data sets and define parameters

%Parameters

%Averaging:
MinParticles=3;         %Minimum number of particles to use to obtain averages
MinEmbryos=3;           %Minimum number of embryos to use to obtain averages

%Rate of elongation
ElongationRate=1.54;    %In kb/minutes.
SEElongationRate=0.14;
GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
Delay=GeneLength/ElongationRate;    %Minutes for PolII to fall off after reaching
                                    %the first MS2 site.

%I'm giving just one data set so it knows whic Dropbox folder to use
[SourcePath,FISHPath,DropboxFolder,MS2CodePath,SchnitzcellsFolder]=...
    DetermineLocalFolders('2013-10-01-SnaESnaP');                             
                                    


% Load the data sets


%Some labels to use for plots
Labels='k.r.g.b.y.c.m.ksrsgsbsyscsmskorogoboyocomok^r^g^b^y^c^m^';


%Some parameters
IntArea=109;        %Area of integration



%SnaE-SnaP
[StatusNum,StatusTxt,StatusRaw]=xlsread([DropboxFolder,filesep,'DataStatusPausing.xlsx'],'SnaSna');

CompileRow=find(strcmp(StatusRaw(:,1),'AnalyzeLiveData Compile Particles'));
CompiledSets=find(strcmp(StatusRaw(CompileRow,:),'READY')|strcmp(StatusRaw(CompileRow,:),'ApproveAll'));

clear SetNames
clear Schnitzcells
clear Ellipses

for i=1:length(CompiledSets)
    SetName=StatusRaw{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    DataSna(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    Schnitzcells(i)=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
    SetNames{i}=SetName;
    %Load Ellipses
    Ellipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);
end

%Now add the extra information
for i=1:length(DataSna)
    DataSna(i).SetName=SetNames{i};
    DataSna(i).schnitzcells=Schnitzcells(i).schnitzcells;
    DataSna(i).Ellipses=Ellipses(i).Ellipses;
end


%SnaE-SogP
[StatusNum,StatusTxt,StatusRaw]=xlsread([DropboxFolder,filesep,'DataStatusPausing.xlsx'],'SnaSog');

CompileRow=find(strcmp(StatusRaw(:,1),'AnalyzeLiveData Compile Particles'));
CompiledSets=find(strcmp(StatusRaw(CompileRow,:),'READY')|strcmp(StatusRaw(CompileRow,:),'ApproveAll'));

clear SetNames
clear Schnitzcells
clear Ellipses

for i=1:length(CompiledSets)
    SetName=StatusRaw{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    DataSog(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    Schnitzcells(i)=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
    SetNames{i}=SetName;
    %Load Ellipses
    Ellipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);
end

%Now add the extra information
for i=1:length(DataSog)
    DataSog(i).SetName=SetNames{i};
    DataSog(i).schnitzcells=Schnitzcells(i).schnitzcells;
    DataSog(i).Ellipses=Ellipses(i).Ellipses;
end


%SnaE-ThsP
[StatusNum,StatusTxt,StatusRaw]=xlsread([DropboxFolder,filesep,'DataStatusPausing.xlsx'],'SnaThs');

CompileRow=find(strcmp(StatusRaw(:,1),'AnalyzeLiveData Compile Particles'));
CompiledSets=find(strcmp(StatusRaw(CompileRow,:),'READY')|strcmp(StatusRaw(CompileRow,:),'ApproveAll'));

clear SetNames
clear Schnitzcells
clear Ellipses

for i=1:length(CompiledSets)
    SetName=StatusRaw{6,CompiledSets(i)};
    Quotes=strfind(SetName,'''');
    Prefix=SetName((Quotes(1)+1):(Quotes(end)-1));
    DataThs(i)=load([DropboxFolder,filesep,Prefix,filesep,'CompiledParticles.mat']);
    Schnitzcells(i)=load([DropboxFolder,filesep,Prefix(1:end),filesep,Prefix(1:end),'_lin.mat']);
    SetNames{i}=SetName;
    %Load Ellipses
    Ellipses(i)=load([DropboxFolder,filesep,Prefix,filesep,'Ellipses.mat']);
end

%Now add the extra information

for i=1:length(DataThs)
    DataThs(i).SetName=SetNames{i};
    DataThs(i).schnitzcells=Schnitzcells(i).schnitzcells;
    DataThs(i).Ellipses=Ellipses(i).Ellipses;
end



%% Compare traces

% %Single traces
% figure(1)
% clf
% hold all
% for i=1:length(DataSna)
%     for j=1:length(DataSna(i).CompiledParticles)
%         if DataSna(i).CompiledParticles(j).nc==14
%             plot(DataSna(i).ElapsedTime(DataSna(i).CompiledParticles(j).Frame)-...
%                 DataSna(i).ElapsedTime(DataSna(i).nc14),...
%                 DataSna(i).CompiledParticles(j).Fluo,'.-')
%         end
%     end
% end
% for i=1:length(DataSog)
%     for j=1:length(DataSog(i).CompiledParticles)
%         if DataSog(i).CompiledParticles(j).nc==14
%             plot(DataSog(i).ElapsedTime(DataSog(i).CompiledParticles(j).Frame)-...
%                 DataSog(i).ElapsedTime(DataSog(i).nc14),...
%                 DataSog(i).CompiledParticles(j).Fluo,'o-','MarkerFaceColor','w')
%         end
%     end
% end
% 
% figure(2)
% clf
% hold all
% for i=1:length(DataThs)
%     for j=1:length(DataThs(i).CompiledParticles)
%         if DataThs(i).CompiledParticles(j).nc==14
%             plot(DataThs(i).ElapsedTime(DataThs(i).CompiledParticles(j).Frame)-...
%                 DataThs(i).ElapsedTime(DataThs(i).nc14),...
%                 DataThs(i).CompiledParticles(j).Fluo,'s-','MarkerFaceColor','w')
%             
%         end
%     end
% end
% hold off      
% 
% 



%Averages

figure(2)
clf
hold all
PlotHandle=[];
PlotHandleForLegend=[];
for i=1:length(DataSna)
        PlotHandle(end+1)=errorbar(DataSna(i).ElapsedTime(DataSna(i).NParticlesAll>=MinParticles)-...
            DataSna(i).ElapsedTime(DataSna(i).nc14),...
            DataSna(i).MeanVectorAll(DataSna(i).NParticlesAll>=MinParticles),...
            DataSna(i).SDVectorAll(DataSna(i).NParticlesAll>=MinParticles)./...
            sqrt(DataSna(i).NParticlesAll(DataSna(i).NParticlesAll>=MinParticles)),'.-');
end
PlotHandleForLegend(end+1)=PlotHandle(end);
for i=1:length(DataSog)
        PlotHandle(end+1)=errorbar(DataSog(i).ElapsedTime(DataSog(i).NParticlesAll>=MinParticles)-...
            DataSog(i).ElapsedTime(DataSog(i).nc14),...
            DataSog(i).MeanVectorAll(DataSog(i).NParticlesAll>=MinParticles),...
            DataSog(i).SDVectorAll(DataSog(i).NParticlesAll>=MinParticles)./...
            sqrt(DataSog(i).NParticlesAll(DataSog(i).NParticlesAll>=MinParticles)),'o-');
end
PlotHandleForLegend(end+1)=PlotHandle(end);
for i=1:length(DataThs)
        PlotHandle(end+1)=errorbar(DataThs(i).ElapsedTime(DataThs(i).NParticlesAll>=MinParticles)-...
            DataThs(i).ElapsedTime(DataThs(i).nc14),...
            DataThs(i).MeanVectorAll(DataThs(i).NParticlesAll>=MinParticles),...
            DataThs(i).SDVectorAll(DataThs(i).NParticlesAll>=MinParticles)./...
            sqrt(DataThs(i).NParticlesAll(DataThs(i).NParticlesAll>=MinParticles)),'s-');
end
hold off
PlotHandleForLegend(end+1)=PlotHandle(end);
xlim([0,60])
legend(PlotHandleForLegend,'Sna','Sog','Ths','Location','SouthEast')
box on
StandardFigure(PlotHandle,gca)

title('Average intensity')

%CV
figure(3)
clf
hold all
for i=1:length(DataSna)
        plot(DataSna(i).ElapsedTime(DataSna(i).NParticlesAll>=MinParticles)-...
            DataSna(i).ElapsedTime(DataSna(i).nc14),...
            DataSna(i).SDVectorAll(DataSna(i).NParticlesAll>=MinParticles)./...
            DataSna(i).MeanVectorAll(DataSna(i).NParticlesAll>=MinParticles),'.-')
end
for i=1:length(DataSog)
        plot(DataSog(i).ElapsedTime(DataSog(i).NParticlesAll>=MinParticles)-...
            DataSog(i).ElapsedTime(DataSog(i).nc14),...
            DataSog(i).SDVectorAll(DataSog(i).NParticlesAll>=MinParticles)./...
            DataSog(i).MeanVectorAll(DataSog(i).NParticlesAll>=MinParticles),'o-')
end
for i=1:length(DataThs)
        plot(DataThs(i).ElapsedTime(DataThs(i).NParticlesAll>=MinParticles)-...
            DataThs(i).ElapsedTime(DataThs(i).nc14),...
            DataThs(i).SDVectorAll(DataThs(i).NParticlesAll>=MinParticles)./...
            DataThs(i).MeanVectorAll(DataThs(i).NParticlesAll>=MinParticles),'s-')
end
hold off
xlim([0,60])

title('CV')


%Fano
figure(4)
clf
hold all
for i=1:length(DataSna)
        plot(DataSna(i).ElapsedTime(DataSna(i).NParticlesAll>=MinParticles)-...
            DataSna(i).ElapsedTime(DataSna(i).nc14),...
            DataSna(i).SDVectorAll(DataSna(i).NParticlesAll>=MinParticles).^2./...
            DataSna(i).MeanVectorAll(DataSna(i).NParticlesAll>=MinParticles),'.-')
end
for i=1:length(DataSog)
        plot(DataSog(i).ElapsedTime(DataSog(i).NParticlesAll>=MinParticles)-...
            DataSog(i).ElapsedTime(DataSog(i).nc14),...
            DataSog(i).SDVectorAll(DataSog(i).NParticlesAll>=MinParticles).^2./...
            DataSog(i).MeanVectorAll(DataSog(i).NParticlesAll>=MinParticles),'o-')
end
for i=1:length(DataThs)
        plot(DataThs(i).ElapsedTime(DataThs(i).NParticlesAll>=MinParticles)-...
            DataThs(i).ElapsedTime(DataThs(i).nc14),...
            DataThs(i).SDVectorAll(DataThs(i).NParticlesAll>=MinParticles).^2./...
            DataThs(i).MeanVectorAll(DataThs(i).NParticlesAll>=MinParticles),'s-')
end
hold off
xlim([0,60])

title('Fanno Factor')

%% First frame
% 
% 
% for i=1:length(DataSna)
%     TimeStartSna{i}=[];
%     for j=1:length(DataSna(i).CompiledParticles)
%         if DataSna(i).CompiledParticles(j).nc==14
%             TimeStartSna{i}=[TimeStartSna{i},...
%                 DataSna(i).ElapsedTime(DataSna(i).CompiledParticles(j).FirstFrame)-...
%                 DataSna(i).ElapsedTime(DataSna(i).nc14)];
%         end
%     end
% end
% for i=1:length(DataSog)
%     TimeStartSog{i}=[];
%     for j=1:length(DataSog(i).CompiledParticles)
%         if DataSog(i).CompiledParticles(j).nc==14
%             TimeStartSog{i}=[TimeStartSog{i},...
%                 DataSog(i).ElapsedTime(DataSog(i).CompiledParticles(j).FirstFrame)-...
%                 DataSog(i).ElapsedTime(DataSog(i).nc14)];
%         end
%     end
% end
% for i=1:length(DataThs)
%     TimeStartThs{i}=[];
%     for j=1:length(DataThs(i).CompiledParticles)
%         if DataThs(i).CompiledParticles(j).nc==14
%             TimeStartThs{i}=[TimeStartThs{i},...
%                 DataThs(i).ElapsedTime(DataThs(i).CompiledParticles(j).FirstFrame)-...
%                 DataThs(i).ElapsedTime(DataThs(i).nc14)];
%         end
%     end
% end
% 
% [nSna{1},xOut]=hist(TimeStartSna{1},[0:5:45]);
% for i=2:length(DataSna)
%     nSna{i}=hist(TimeStartSna{i},xOut);
% end
% for i=1:length(DataThs)
%     nThs{i}=hist(TimeStartThs{i},xOut);
% end
% for i=1:length(DataSog)
%     nSog{i}=hist(TimeStartSog{i},xOut);
% end
% 
% figure(1)
% clf
% hold all
% for i=1:length(DataSna)
%     plot(xOut,nSna{i}/sum(nSna{i}),'.-')
% end
% for i=1:length(DataSog)
%     plot(xOut,nSog{i}/sum(nSog{i}),'o-')
% end
% for i=1:length(DataThs)
%     plot(xOut,nThs{i}/sum(nThs{i}),'s-')
% end
% hold off        
%         



%% Fraction of active nuclei


% 
% figure(5)
% clf
% hold all
% for i=1:length(DataSna)
%         plot(DataSna(i).ElapsedTime-...
%             DataSna(i).ElapsedTime(DataSna(i).nc14),...
%             DataSna(i).NParticlesAll/sum(DataSna(i).ncFilter(:,end)),'.-')
% end
% for i=1:length(DataSog)
%         plot(DataSog(i).ElapsedTime-...
%             DataSog(i).ElapsedTime(DataSog(i).nc14),...
%             DataSog(i).NParticlesAll/sum(DataSog(i).ncFilter(:,end)),'o-')
% end
% for i=1:length(DataThs)
%         plot(DataThs(i).ElapsedTime-...
%             DataThs(i).ElapsedTime(DataThs(i).nc14),...
%             DataThs(i).NParticlesAll/sum(DataThs(i).ncFilter(:,end)),'s-')
% end
% hold off
% xlim([0,45])
% 

% Threshold=0.5E4;
% 
% figure(2)
% clf
% hold all
% for i=1:length(DataSna)
%     NumberOnSna{i}=zeros(length(DataSna(i).ElapsedTime),1);
%     for Frame=DataSna(i).nc14:length(DataSna(i).ElapsedTime)
%         for j=1:length(DataSna(i).CompiledParticles)
%             if sum(DataSna(i).CompiledParticles(j).Frame==Frame)
%                 if DataSna(i).CompiledParticles(j).Fluo(DataSna(i).CompiledParticles(j).Frame==Frame)>=...
%                         Threshold
%                     NumberOnSna{i}(Frame)=NumberOnSna{i}(Frame)+1;
%                 end
%             end
%         end
%     end
%     plot(DataSna(i).ElapsedTime-...
%         DataSna(i).ElapsedTime(DataSna(i).nc14),...
%         NumberOnSna{i}/max(NumberOnSna{i}),'.-')
% end
% for i=1:length(DataSog)
%     NumberOnSog{i}=zeros(length(DataSog(i).ElapsedTime),1);
%     for Frame=DataSog(i).nc14:length(DataSog(i).ElapsedTime)
%         for j=1:length(DataSog(i).CompiledParticles)
%             if sum(DataSog(i).CompiledParticles(j).Frame==Frame)
%                 if DataSog(i).CompiledParticles(j).Fluo(DataSog(i).CompiledParticles(j).Frame==Frame)>=...
%                         Threshold
%                     NumberOnSog{i}(Frame)=NumberOnSog{i}(Frame)+1;
%                 end
%             end
%         end
%     end
%     plot(DataSog(i).ElapsedTime-...
%         DataSog(i).ElapsedTime(DataSog(i).nc14),...
%         NumberOnSog{i}/max(NumberOnSog{i}),'o-')
% end
% for i=1:length(DataThs)
%     NumberOnThs{i}=zeros(length(DataThs(i).ElapsedTime),1);
%     for Frame=DataThs(i).nc14:length(DataThs(i).ElapsedTime)
%         for j=1:length(DataThs(i).CompiledParticles)
%             if sum(DataThs(i).CompiledParticles(j).Frame==Frame)
%                 if DataThs(i).CompiledParticles(j).Fluo(DataThs(i).CompiledParticles(j).Frame==Frame)>=...
%                         Threshold
%                     NumberOnThs{i}(Frame)=NumberOnThs{i}(Frame)+1;
%                 end
%             end
%         end
%     end
%     plot(DataThs(i).ElapsedTime-...
%         DataThs(i).ElapsedTime(DataThs(i).nc14),...
%         NumberOnThs{i}/max(NumberOnThs{i}),'s-')
% end
% hold off
%                 
% %%%%%%%% Cumulative Plot of Data
% 
% ThsSort=[];
% 
% for i=1:length(DataThs);
%  
%     
% Ths(i).Time=DataThs(i).ElapsedTime-...
% DataThs(i).ElapsedTime(DataThs(i).nc14);
% 
% Ths(i).Frac = DataThs(i).NParticlesAll/sum(DataThs(i).ncFilter(:,end));
% 
% Ths(i).Sort=sort(Ths(i).Frac(find(Ths(i).Time>=0&Ths(i).Time<=40)));
% 
% ThsSort=sort([ThsSort,Ths(i).Sort]);
% 
% end
% 
% ThsX = linspace(0,1,length(ThsSort));
% 
% 
% SnaSort=[];
% 
% for i=1:length(DataSna);
%     
% Sna(i).Time=DataSna(i).ElapsedTime-...
% DataSna(i).ElapsedTime(DataSna(i).nc14);
% 
% Sna(i).Frac = DataSna(i).NParticlesAll/sum(DataSna(i).ncFilter(:,end));
% 
% Sna(i).Sort=sort(Sna(i).Frac(find(Sna(i).Time>=0&Sna(i).Time<=40)));
% 
% SnaSort=sort([SnaSort,Sna(i).Sort]);
% 
% end
% 
% SnaX = linspace(0,1,length(SnaSort));
% 
% SogSort=[];
% 
% for i=1:length(DataSog);
%     
% Sog(i).Time=DataSog(i).ElapsedTime-...
% DataSog(i).ElapsedTime(DataSog(i).nc14);
% 
% Sog(i).Frac = DataSog(i).NParticlesAll/sum(DataSog(i).ncFilter(:,end));
% 
% Sog(i).Sort=sort(Sog(i).Frac(find(Sog(i).Time>=0&Sog(i).Time<=40)));
% 
% SogSort=sort([SogSort,Sog(i).Sort]);
% 
% end
% 
% SogX = linspace(0,1,length(SogSort));
% 
% 
% figure(6)
% hold all
% 
% plot(SnaX,SnaSort,'r.')
% plot(SogX,SogSort,'g.')
% plot(ThsX,ThsSort,'b.')
% 
% hold off


%%%%%%%%%

% Single traces handful
% 
% figure(7)
% clf
% hold all
% Counter=1;
% CounterMax=5;
% 
% for i=1:length(DataSna)
%     for j=1:length(DataSna(i).CompiledParticles)
%         if DataSna(i).CompiledParticles(j).nc==14 & Counter<=CounterMax
%             plot(DataSna(i).ElapsedTime(DataSna(i).CompiledParticles(j).Frame)-...
%                 DataSna(i).ElapsedTime(DataSna(i).nc14),...
%                 DataSna(i).CompiledParticles(j).Fluo,'.-')
%             Counter=Counter+1;
%         end
%     end
% end
% 
% Counter=1;
% for i=1:length(DataSog)
%     for j=1:length(DataSog(i).CompiledParticles)
%         if DataSog(i).CompiledParticles(j).nc==14 & Counter<=CounterMax
%             plot(DataSog(i).ElapsedTime(DataSog(i).CompiledParticles(j).Frame)-...
%                 DataSog(i).ElapsedTime(DataSog(i).nc14),...
%                 DataSog(i).CompiledParticles(j).Fluo,'o-','MarkerFaceColor','w')
%                    Counter=Counter+1;
%         end
%     end
% end
% 
% Counter=1;
% for i=1:length(DataThs)
%     for j=1:length(DataThs(i).CompiledParticles)
%         if DataThs(i).CompiledParticles(j).nc==14 & Counter<=CounterMax
%             plot(DataThs(i).ElapsedTime(DataThs(i).CompiledParticles(j).Frame)-...
%                 DataThs(i).ElapsedTime(DataThs(i).nc14),...
%                 DataThs(i).CompiledParticles(j).Fluo,'s-','MarkerFaceColor','w')
%                 Counter=Counter+1;
%         end
%     end
% end
% hold off
% xlim([0,45])

% 
% %% First Detection Times
% 
% figure(8)
% clf
% hold all
% 
% FirstPasageTimesSna=[];
% 
% for i=1:length(DataSna)
%     for j=1:length(DataSna(i).CompiledParticles)
%         if DataSna(i).CompiledParticles(j).nc==14
%            
%             FirstPasageTimesSna=[FirstPasageTimesSna; DataSna(i).ElapsedTime(DataSna(i).CompiledParticles(j).FirstFrame)-...
%                 DataSna(i).ElapsedTime(DataSna(i).nc14)];
%             
%         end
%     end
% end
% 
% FirstPasageTimesSog=[];
% 
% for i=1:length(DataSog)
%     for j=1:length(DataSog(i).CompiledParticles)
%         if DataSog(i).CompiledParticles(j).nc==14
% 
%              FirstPasageTimesSog=[FirstPasageTimesSog; DataSog(i).ElapsedTime(DataSog(i).CompiledParticles(j).FirstFrame)-...
%                 DataSog(i).ElapsedTime(DataSog(i).nc14)];
%             
%         end
%     end
% end
% 
% 
% FirstPasageTimesThs=[];
% 
%  for i=1:length(DataThs)
%      for j=1:length(DataThs(i).CompiledParticles)
%          if DataThs(i).CompiledParticles(j).nc==14
%            
%           FirstPasageTimesThs=[FirstPasageTimesThs; DataThs(i).ElapsedTime(DataThs(i).CompiledParticles(j).FirstFrame)-...
%           DataThs(i).ElapsedTime(DataThs(i).nc14)];
%             
%              
%          end
%      end
%  end
%  
%  xHist = [0:2:30];
%  
%  NSna = hist(FirstPasageTimesSna,[0:2:30]);
%  NSog = hist(FirstPasageTimesSog,[0:2:30]);
%  NThs = hist(FirstPasageTimesThs,[0:2:30]);
%  
%  plot(xHist(1:end-1),NSna(1:end-1)./sum(NSna(1:end-1)),'r');
%  plot(xHist(1:end-1),NSog(1:end-1)./sum(NSog(1:end-1)),'g');
%  plot(xHist(1:end-1),NThs(1:end-1)./sum(NThs(1:end-1)),'b');
%  
%  ylabel('Scaled Histogram number of cells')
%  xlabel('Time at first expression in cc14')
%  hold off      

save_to_base(1)

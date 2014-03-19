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

%Single traces
figure
clf
hold all

NumberSnail(1).Num=[3,4,7,8,11,12,13,17,18,19,20,33]; 
NumberSnail(2).Num=[18,20,45,60,68];
NumberSnail(3).Num=[41,43,46,49,50,51,52,59,62,63];

for i=1:length(DataSna)
    for j=1:length(DataSna(i).CompiledParticles)
        if DataSna(i).CompiledParticles(j).nc==14 & ismember(j,NumberSnail(i).Num)
            plot(DataSna(i).ElapsedTime(DataSna(i).CompiledParticles(j).Frame)-...
                DataSna(i).ElapsedTime(DataSna(i).nc14),...
                DataSna(i).CompiledParticles(j).Fluo,'.-')
        end
    end
end

hold off

figure
clf
hold all

NumberSog(1).Num=[42,43,49,50,51,52,54,57,58,59,60,68,69]; 
NumberSog(2).Num=[22,23,25,28,32,33,34,35,36,45,47,54,59];
NumberSog(3).Num=[27,29,32,37,38,40,46,47,48,51,57];

 for i=1:length(DataSog)
     for j=1:length(DataSog(i).CompiledParticles)
         if DataSog(i).CompiledParticles(j).nc==14 & ismember(j,NumberSog(i).Num)
             plot(DataSog(i).ElapsedTime(DataSog(i).CompiledParticles(j).Frame)-...
                 DataSog(i).ElapsedTime(DataSog(i).nc14),...
                 DataSog(i).CompiledParticles(j).Fluo,'o-','MarkerFaceColor','w')
         end
     end
 end
 
 figure
 clf
 hold all
 
NumberThs(1).Num=[1,2,3,5,7,8,9,11,13,24,26]; 
NumberThs(2).Num=[2,3,7,20,23,24,31];
NumberThs(3).Num=[22,28,30,34,50,52,64,79];
 
 
for i=1:length(DataThs)
    for j=1:length(DataThs(i).CompiledParticles)
        if DataThs(i).CompiledParticles(j).nc==14 & ismember(j,NumberThs(i).Num)
            plot(DataThs(i).ElapsedTime(DataThs(i).CompiledParticles(j).Frame)-...
                DataThs(i).ElapsedTime(DataThs(i).nc14),...
                DataThs(i).CompiledParticles(j).Fluo,'s-','MarkerFaceColor','w')
        end
    end
end
hold off      





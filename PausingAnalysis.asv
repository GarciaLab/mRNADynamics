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
    DetermineLocalFolders('2013-10-01-SnaESnaP')     ;                             
                                    


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
figure(1)
clf
hold all
for i=1:length(DataSna)
    for j=1:length(DataSna(i).CompiledParticles)
        if DataSna(i).CompiledParticles(j).nc==14
            plot(DataSna(i).ElapsedTime(DataSna(i).CompiledParticles(j).Frame)-...
                DataSna(i).ElapsedTime(DataSna(i).nc14),...
                DataSna(i).CompiledParticles(j).Fluo,'.-')
        end
    end
end
for i=1:length(DataSog)
    for j=1:length(DataSog(i).CompiledParticles)
        if DataSog(i).CompiledParticles(j).nc==14
            plot(DataSog(i).ElapsedTime(DataSog(i).CompiledParticles(j).Frame)-...
                DataSog(i).ElapsedTime(DataSog(i).nc14),...
                DataSog(i).CompiledParticles(j).Fluo,'o-','MarkerFaceColor','w')
        end
    end
end
for i=1:length(DataThs)
    for j=1:length(DataThs(i).CompiledParticles)
        if DataThs(i).CompiledParticles(j).nc==14
            plot(DataThs(i).ElapsedTime(DataThs(i).CompiledParticles(j).Frame)-...
                DataThs(i).ElapsedTime(DataThs(i).nc14),...
                DataThs(i).CompiledParticles(j).Fluo,'s-','MarkerFaceColor','w')
        end
    end
end
hold off      





%Averages
figure(2)
clf
hold all
for i=1:length(DataSna)
        errorbar(DataSna(i).ElapsedTime(DataSna(i).NParticlesAll>=MinParticles)-...
            DataSna(i).ElapsedTime(DataSna(i).nc14),...
            DataSna(i).MeanVectorAll(DataSna(i).NParticlesAll>=MinParticles),...
            DataSna(i).SDVectorAll(DataSna(i).NParticlesAll>=MinParticles)./...
            sqrt(DataSna(i).NParticlesAll(DataSna(i).NParticlesAll>=MinParticles)),'.-')
end
for i=1:length(DataSog)
        errorbar(DataSog(i).ElapsedTime(DataSog(i).NParticlesAll>=MinParticles)-...
            DataSog(i).ElapsedTime(DataSog(i).nc14),...
            DataSog(i).MeanVectorAll(DataSog(i).NParticlesAll>=MinParticles),...
            DataSog(i).SDVectorAll(DataSog(i).NParticlesAll>=MinParticles)./...
            sqrt(DataSog(i).NParticlesAll(DataSog(i).NParticlesAll>=MinParticles)),'o-')
end
for i=1:length(DataThs)
        errorbar(DataThs(i).ElapsedTime(DataThs(i).NParticlesAll>=MinParticles)-...
            DataThs(i).ElapsedTime(DataThs(i).nc14),...
            DataThs(i).MeanVectorAll(DataThs(i).NParticlesAll>=MinParticles),...
            DataThs(i).SDVectorAll(DataThs(i).NParticlesAll>=MinParticles)./...
            sqrt(DataThs(i).NParticlesAll(DataThs(i).NParticlesAll>=MinParticles)),'s-')
end
hold off
xlim([0,60])



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



%% Fraction of active nuclei

close all

figure(1)
clf
hold all
for i=1:length(DataSna)
        plot(DataSna(i).ElapsedTime-...
            DataSna(i).ElapsedTime(DataSna(i).nc14),...
            DataSna(i).NParticlesAll/sum(DataSna(i).ncFilter(:,end)),'.-')
end
for i=1:length(DataSog)
        plot(DataSog(i).ElapsedTime-...
            DataSog(i).ElapsedTime(DataSog(i).nc14),...
            DataSog(i).NParticlesAll/sum(DataSog(i).ncFilter(:,end)),'o-')
end
for i=1:length(DataThs)
        plot(DataThs(i).ElapsedTime-...
            DataThs(i).ElapsedTime(DataThs(i).nc14),...
            DataThs(i).NParticlesAll/sum(DataThs(i).ncFilter(:,end)),'s-')
end
hold off
xlim([0,45])





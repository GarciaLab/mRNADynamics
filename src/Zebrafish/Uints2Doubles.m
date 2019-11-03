function [Spots] = Uints2Doubles(varargin)
%Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;

[Prefix, ForceAP, SkipTraces, SkipFluctuations, SkipFits, SkipMovie, ...
    SkipAll, ApproveAll, MinParticles, minTime, ROI, intArea, noHist, ...
    ROI1, ROI2, slimVersion, manualSingleFits, optionalResults, yToManualAlignmentPrompt, minBinSize, edgeWidth] = determineCompileParticlesOptions(varargin);

FilePrefix=[Prefix,'_'];

[rawDataPath,ProcPath,DropboxFolder,MS2CodePath, PreProcPath,...
    rawDataFolder, Prefix, ExperimentType,Channel1,Channel2,OutputFolder,...
    Channel3, spotChannels, MovieDataBaseFolder, movieDatabase]...
    = readMovieDatabase(Prefix, optionalResults);

disp('Loading Spots.mat...');
load([DropboxFolder,filesep,Prefix,filesep,'Spots.mat']);

for i=1:length(Spots);
    if isempty(Spots(i).Fits);
    else  
    Fields=fieldnames(Spots(i).Fits);
    break
    end
end

    Field16s=zeros(length(Fields));
    Field8s=zeros(length(Fields));

for i=1:length(Spots);
    for j=1:length(Spots(i).Fits);
        for k=1:length(Fields);
        if isa(Spots(i).Fits(j).(Fields{k}),'uint16');
            Spots(i).Fits(j).(Fields{k})=double(Spots(i).Fits(j).(Fields{k}));
            Field16s(k)=1;
            
        end
     end
    end
end

for i=1:length(Spots);
    for j=1:length(Spots(i).Fits);
        for k=1:length(Fields);
        if isa(Spots(i).Fits(j).(Fields{k}),'uint8');
            Spots(i).Fits(j).(Fields{k})=double(Spots(i).Fits(j).(Fields{k}));
             Field8s(k)=1;
        end
     end
    end
end
 disp('Saving Spots.mat...');
 save([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'],'Spots');
  save([DropboxFolder,filesep,Prefix,filesep,'SpotFields.mat'],'Field16s','Field8s');
 
 disp('Loading Schnitzcells.mat...');
 load([DropboxFolder,filesep,Prefix,filesep,Prefix '_lin.mat']);
 

 FieldSchs=fieldnames(schnitzcells);
 Field16Schs=zeros(length(FieldSchs));

for i=1:length(schnitzcells);
    for j=1:length(FieldSchs);
        if isa(schnitzcells(i).(FieldSchs{j}),'uint16')
            schnitzcells(i).(FieldSchs{j})=double(schnitzcells(i).(FieldSchs{j}));
            Field16Schs(j)=1;
        end
     end
end

 disp('Saving Schnitzcells.mat...'); 
 save([DropboxFolder,filesep,Prefix,filesep,Prefix '_lin.mat'],'schnitzcells');
  save([DropboxFolder,filesep,Prefix,filesep,'SchFields'],'Field16Schs')
end

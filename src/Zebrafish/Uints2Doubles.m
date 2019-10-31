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


Fields=fieldnames(Spots(1).Fits);

for i=1:length(Spots);
    for j=1:length(Spots(i).Fits);
        for k=1:length(Fields);
        if isa(Spots(i).Fits(j).(Fields{k}),'uint16')|| isa(Spots(i).Fits(j).(Fields{k}),'uint8');
            Spots(i).Fits(j).(Fields{k})=double(Spots(i).Fits(j).(Fields{k}));
        end
     end
    end
end
 disp('Saving Spots.mat...');
 save([DropboxFolder,filesep,Prefix,filesep,'Spots.mat'],'Spots');

 disp('Loading Schnitzcells.mat...');
 load([DropboxFolder,filesep,Prefix,filesep,Prefix '_lin.mat']);
 
 clear Fields;
 Fields=fieldnames(schnitzcells);

for i=1:length(schnitzcells);
    for j=1:length(Fields);
        if isa(schnitzcells(i).(Fields{j}),'uint16')
            schnitzcells(i).(Fields{j})=double(schnitzcells(i).(Fields{j}));
        end
     end
end

 disp('Saving Schnitzcells.mat...'); 
 save([DropboxFolder,filesep,Prefix,filesep,Prefix '_lin.mat'],'schnitzcells')
end

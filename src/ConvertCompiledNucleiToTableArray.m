function [CompiledNucleiTable] = ConvertCompiledNucleiToTableArray(Prefix, varargin)
% ConvertCompiledNucleiToTableArray(Prefix, varargin)
%
% DESCRIPTION
% This function converts the struct array output by the fucntion 
% CompileNuclearProtein to a table array format. 
%
% ARGUMENTS
% Prefix: prefix string of the data set to analyze. 
% varargin: A cell in which the first element nis the prefix string of the data set
%           to analyze. Subsequent elements can be the options below.
% 
% Author (contact): Gabriella Martini (martini@berkeley.edu)
% Created: 11/29/19
% Last Updated: 12/04/19
%
% Documented by: Gabriella Martini (martini@berkeley.edu)


%% Get Information about about folders
[~,~,DefaultDropboxFolder,~,~]=...
    DetermineLocalFolders;


%Load the folder information
[SourcePath, FISHPath, DefaultDropboxFolder, DropboxFolder, MS2CodePath, PreProcPath,...
configValues, movieDatabasePath] = DetermineAllLocalFolders(Prefix);

[Date, ExperimentType, ExperimentAxis, CoatProtein, StemLoop, APResolution,...
Channel1, Channel2, Objective, Power, DataFolder, DropboxFolderName, Comments,...
nc9, nc10, nc11, nc12, nc13, nc14, CF] = getExperimentDataFromMovieDatabase(Prefix, DefaultDropboxFolder);



%Load all the information
disp('Loading CompiledNuclei.mat...');
load([DropboxFolder,filesep,Prefix,filesep,'CompiledNuclei.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'FrameInfo.mat']);
load([DropboxFolder,filesep,Prefix,filesep,'APDivision.mat']);

APDivisionTimes = zeros(size(APDivision));
for i=1:size(APDivision, 1)
    for j = 1:size(APDivision, 2)
        if APDivision(i,j) ~= 0
            APDivisionTimes(i,j) = FrameInfo(APDivision(i,j)).Time;
        end
    end
end


% 
%% 
numNuclei = size(CompiledNuclei, 2);
varnames =  { 'Prefix','NucleusID', 'Frame', 'Time', 'TimeNC', 'FrameCount', 'xPos',...
    'yPos', 'Radius', 'cellno', 'nc', 'schnitz', 'MeanAP', 'MedianAP', 'APbin',...
    'MeanDV', 'MedianDV', 'Fluo', 'FluoMax', 'ncStart'};
vartypes = { 'string','int64', 'int64', 'double', 'double', 'int64', 'int64',...
    'int64', 'double', 'int64', 'int64', 'int64', 'double', 'double','double',...
    'double', 'double', 'double', 'double', 'int64'};


nvars = length(varnames);

compnucfields = fieldnames(CompiledNuclei(1));
n = 1;
totalSamples = 0;
while n <= numNuclei
    totalSamples = totalSamples + size(CompiledNuclei(n).Frames, 1);
    n = n + 1;
end

CompiledNucleiTable = table('Size', [totalSamples nvars],...
    'VariableTypes', vartypes,'VariableNames', varnames);

sample = 1;
n = 1;
nc_frames = [nc9, nc10, nc11, nc12, nc13, nc14];
while n <= numNuclei
     FrameCount = size(CompiledNuclei(n).Frames, 1);
     nc = CompiledNuclei(n).nc;
     schnitz = CompiledNuclei(n).schnitz;
     MeanAP = CompiledNuclei(n).MeanAP;
     MedianAP = CompiledNuclei(n).MedianAP;
     APbin = round(MeanAP/APResolution)*APResolution;
     MeanDV = CompiledNuclei(n).MeanDV;
     MedianDV = CompiledNuclei(n).MedianDV;
     ncStart = nc_frames(nc-8);
     ncStart = max(ncStart, 1);
     ncOffset = FrameInfo(ncStart).Time/60;
     for f=1:size(CompiledNuclei(n).Frames, 1)
        CompiledNucleiTable.Prefix(sample) = Prefix;
        CompiledNucleiTable.NucleusID(sample) = n; 
        CompiledNucleiTable.MeanAP(sample) = MeanAP;
        CompiledNucleiTable.MedianAP(sample) = MedianAP;
        CompiledNucleiTable.APbin(sample) = APbin;
        CompiledNucleiTable.Frame(sample) = CompiledNuclei(n).Frames(f);
        CompiledNucleiTable.FrameNC(sample) = CompiledNuclei(n).Frames(f)-APDivision(uint16(nc),uint16(APbin/APResolution)+1);
        
        CompiledNucleiTable.Time(sample) = FrameInfo(CompiledNucleiTable.Frame(sample)).Time/60;
        %CompiledNucleiTable.TimeNC(sample) = CompiledNucleiTable.Time(sample)-ncOffset;
        CompiledNucleiTable.TimeNC(sample) = CompiledNucleiTable.Time(sample)-APDivisionTimes(uint16(nc), uint16(APbin/APResolution))/60;
        CompiledNucleiTable.FrameCount(sample) = FrameCount;
        CompiledNucleiTable.xPos(sample) = CompiledNuclei(n).xPos(f);
        CompiledNucleiTable.yPos(sample) = CompiledNuclei(n).yPos(f);
        CompiledNucleiTable.Radius(sample) = CompiledNuclei(n).Radius(f);
        CompiledNucleiTable.cellno(sample) = CompiledNuclei(n).cellno(f);
        CompiledNucleiTable.nc(sample) = nc;
        CompiledNucleiTable.schnitz(sample) = schnitz;
        
        CompiledNucleiTable.MeanDV(sample) = MeanDV;
        CompiledNucleiTable.MedianDV(sample) = MedianDV;
        CompiledNucleiTable.Fluo(sample) = CompiledNuclei(n).FluoTimeTrace(f);
        CompiledNucleiTable.FluoMax(sample) = CompiledNuclei(n).FluoMax(f);            
        CompiledNucleiTable.ncStart(sample) = ncStart;
        sample = sample + 1;
     end
     n = n + 1;
end

save([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat'],...
        'CompiledNucleiTable');

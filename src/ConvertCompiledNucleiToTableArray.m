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
    'MeanDV', 'MedianDV', 'ncStart', 'Fluo', 'FluoMax', 'Fluo2', 'MaxSlice2', 'Fluo2Delta1','Time2Delta1',  'Fluo2Delta2','Time2Delta2',...
    'Fluo2Delta3', 'Time2Delta3','Fluo2Delta4','Time2Delta4', 'Fluo2Delta5','Time2Delta5', 'Fluo2Delta6','Time2Delta6', 'Fluo2Delta7','Time2Delta7', 'Fluo2Delta8','Time2Delta8',...
    'Fluo2Delta9', 'Time2Delta9', 'Fluo2Delta10', 'Time2Delta10', 'Flag1', 'Flag2', 'Flag3', 'Flag4'};
vartypes = { 'string','int64', 'int64', 'double', 'double', 'int64', 'int64',...
    'int64', 'double', 'int64', 'int64', 'int64', 'double', 'double','double',...
    'double', 'double','int64', 'double', 'double', 'double', 'int64', 'double', 'double','double', 'double',...
    'double', 'double','double', 'double','double', 'double','double', 'double','double', 'double','double', 'double',...
    'double', 'double','double', 'double','logical', 'logical', 'logical', 'logical'};


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
NCDeltas = {[], [], [], [], [], []};

h=waitbar(0,'Converting Compiled Nuclei');
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
     FrameTimes = vertcat(FrameInfo.Time)/60;
     NCFrameTimes = FrameTimes(CompiledNuclei(n).Frames)-ncOffset;
     
     if (nc > 12) & (length(NCFrameTimes) < 20)
         Flag1 = true;
     elseif (nc <= 12) & (length(NCFrameTimes) < 10)
         Flag1 = true;
     else
         Flag1 = false;
     end
     if min(NCFrameTimes) > 10
         Flag3 = true;
     else
         Flag3 = false;
     end
     NucleusFluoDeltas = diff(CompiledNuclei(n).FluoTimeTrace2);
     TimeDeltas = diff(NCFrameTimes);
     NucleusDeltas = NucleusFluoDeltas./(TimeDeltas.');
     NCDeltas{nc-8} = [NCDeltas{nc-8} NucleusDeltas];

     for f=1:FrameCount
        waitbar(sample/totalSamples,h)
        CompiledNucleiTable.Prefix(sample) = Prefix;
        CompiledNucleiTable.NucleusID(sample) = n; 
        CompiledNucleiTable.MeanAP(sample) = MeanAP;
        CompiledNucleiTable.MedianAP(sample) = MedianAP;
        CompiledNucleiTable.APbin(sample) = APbin;
        CompiledNucleiTable.Frame(sample) = CompiledNuclei(n).Frames(f);
        CompiledNucleiTable.FrameNC(sample) = CompiledNuclei(n).Frames(f)-APDivision(uint16(nc),uint16(APbin/APResolution)+1);
        
        CompiledNucleiTable.Time(sample) = FrameInfo(CompiledNucleiTable.Frame(sample)).Time/60;
        %CompiledNucleiTable.TimeNC(sample) = CompiledNucleiTable.Time(sample)-ncOffset;
        CompiledNucleiTable.TimeNC(sample) = CompiledNucleiTable.Time(sample)-APDivisionTimes(uint16(nc), uint16(APbin/APResolution)+1)/60;
        CompiledNucleiTable.FrameCount(sample) = FrameCount;
        CompiledNucleiTable.xPos(sample) = CompiledNuclei(n).xPos(f);
        CompiledNucleiTable.yPos(sample) = CompiledNuclei(n).yPos(f);
        CompiledNucleiTable.Radius(sample) = CompiledNuclei(n).Radius(f);
        CompiledNucleiTable.cellno(sample) = CompiledNuclei(n).cellno(f);
        CompiledNucleiTable.nc(sample) = nc;
        CompiledNucleiTable.schnitz(sample) = schnitz;
        
        CompiledNucleiTable.MeanDV(sample) = MeanDV;
        CompiledNucleiTable.MedianDV(sample) = MedianDV;
        CompiledNucleiTable.ncStart(sample) = ncStart;
        CompiledNucleiTable.Fluo(sample) = CompiledNuclei(n).FluoTimeTrace(f);
        CompiledNucleiTable.FluoMax(sample) = CompiledNuclei(n).FluoMax(f);
        CompiledNucleiTable.Fluo2(sample) = CompiledNuclei(n).FluoTimeTrace2(f);
        CompiledNucleiTable.MaxSlice2(sample) = CompiledNuclei(n).MaxSlice2(f);
        
        CompiledNucleiTable.Flag1(sample) = Flag1;
        CompiledNucleiTable.Flag2(sample) = CompiledNuclei(n).Flag2(f);
        CompiledNucleiTable.Flag3(sample) = Flag3;
        
        sample = sample + 1;
     end
     n = n + 1;
end

for nc=9:14
   if (size(CompiledNucleiTable(CompiledNucleiTable.nc == nc,:), 1) > 0) &...
      (length(NCDeltas{nc-8}) > 100)
        deltas = abs(NCDeltas{nc-8});
        maxdelta = mean(deltas) + 3*std(deltas);
        uniqueNCIDs = unique(CompiledNucleiTable.NucleusID);
        for i=1:length(uniqueNCIDs)
            NCID = uniqueNCIDs(i);
            NCID_idx = find(CompiledNucleiTable.NucleusID == NCID);
            Fluos = CompiledNucleiTable.Fluo2(NCID_idx);
            ncidTimes = CompiledNucleiTable.TimeNC(NCID_idx);
            FluoTab = CompiledNucleiTable(NCID_idx, {'Fluo2Delta1','Fluo2Delta2', 'Fluo2Delta3',...
                'Fluo2Delta4', 'Fluo2Delta5', 'Fluo2Delta6', 'Fluo2Delta7', 'Fluo2Delta8', 'Fluo2Delta9', 'Fluo2Delta10'}); 
            TimeTab = CompiledNucleiTable(NCID_idx, {'Time2Delta1','Time2Delta2', 'Time2Delta3',...
                'Time2Delta4', 'Time2Delta5', 'Time2Delta6', 'Time2Delta7', 'Time2Delta8', 'Time2Delta9', 'Time2Delta10'}); 
            
            CompiledNucleiTable.Flag4(NCID_idx(1)) = false;
            if length(NCID_idx) == 1
                continue
            end
            prevbadframes = 0;
            for j = 2:length(NCID_idx)
                if ncidTimes(j) < 10
                    continue 
                end
                if prevbadframes >= 10
                    CompiledNucleiTable.Flag4(NCID_idx(j)) = true;
                    continue
                elseif abs((Fluos(j)-Fluos(j-(1+prevbadframes)))/(ncidTimes(j)-ncidTimes(j-(1+prevbadframes)))) > maxdelta
                    CompiledNucleiTable.Flag4(NCID_idx(j)) = true;
                    prevbadframes = prevbadframes+1;
                else
                    CompiledNucleiTable.Flag4(NCID_idx(j))= false;
                    prevbadframes = 0;
                end     
            end
        end
   end
end

close all


save([DropboxFolder,filesep,Prefix,filesep,'CompiledNucleiTable.mat'],...
        'CompiledNucleiTable');
save([DropboxFolder,filesep,Prefix,filesep,'NCDeltas.mat'],...
    'NCDeltas');

